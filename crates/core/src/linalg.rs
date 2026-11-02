//! Dense linear algebra on top of `faer`.
//!
//! Influence matrices are stored column-major so that each column (the
//! influence of one singularity element on every collocation point) can be
//! filled by an independent rayon task, then factored in place by faer's
//! blocked, multithreaded partial-pivoting LU without an extra copy.

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::lu::partial_pivoting::{factor, solve};
use faer::{MatMut, MatRef, Par};
use rayon::prelude::*;

#[derive(Debug, Clone, PartialEq)]
pub enum LinalgError {
    /// The matrix is (numerically) singular; `col` is the first bad pivot.
    Singular { col: usize },
    /// The matrix contains NaN or infinite entries.
    NonFinite,
    /// Eigen solver did not converge.
    NoConvergence,
}

impl std::fmt::Display for LinalgError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LinalgError::Singular { col } => write!(f, "matrix is singular (pivot {col} is zero): check for duplicated or zero-area panels"),
            LinalgError::NonFinite => write!(f, "matrix contains NaN/inf values: check the geometry for degenerate panels"),
            LinalgError::NoConvergence => write!(f, "eigenvalue solver did not converge"),
        }
    }
}
impl std::error::Error for LinalgError {}

fn par() -> Par {
    Par::rayon(0)
}

/// Column-major dense matrix.
#[derive(Clone, Debug, Default)]
pub struct DenseMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl DenseMatrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        DenseMatrix { rows, cols, data: vec![0.0; rows * cols] }
    }
    #[inline]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.rows + i]
    }
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[j * self.rows + i] = v;
    }
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.data[j * self.rows + i] += v;
    }
    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }
    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }
    /// Fill every column in parallel: `f(j, column_j)`.
    pub fn fill_columns_par(&mut self, f: impl Fn(usize, &mut [f64]) + Sync + Send) {
        let r = self.rows.max(1);
        self.data.par_chunks_mut(r).enumerate().for_each(|(j, c)| f(j, c));
    }
    /// Matrix-vector product.
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.rows];
        for j in 0..self.cols {
            let xj = x[j];
            if xj != 0.0 {
                for (yi, a) in y.iter_mut().zip(self.col(j)) {
                    *yi += a * xj;
                }
            }
        }
        y
    }
    fn as_ref(&self) -> MatRef<'_, f64> {
        MatRef::from_column_major_slice(&self.data, self.rows, self.cols)
    }
    fn as_mut(&mut self) -> MatMut<'_, f64> {
        MatMut::from_column_major_slice_mut(&mut self.data, self.rows, self.cols)
    }
    pub fn from_rows(rows: &[Vec<f64>]) -> Self {
        let r = rows.len();
        let c = if r > 0 { rows[0].len() } else { 0 };
        let mut m = DenseMatrix::zeros(r, c);
        for (i, row) in rows.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                m.set(i, j, *v);
            }
        }
        m
    }
}

/// LU factorisation with partial pivoting (in place).
pub struct Lu {
    lu: DenseMatrix,
    fwd: Vec<usize>,
    bwd: Vec<usize>,
    /// estimated reciprocal "pivot growth" (min|u_ii| / max|u_ii|)
    pub pivot_ratio: f64,
}

impl Lu {
    /// Factor a square matrix, consuming it.
    pub fn factor(mut a: DenseMatrix) -> Result<Lu, LinalgError> {
        assert_eq!(a.rows, a.cols, "LU needs a square matrix");
        let n = a.rows;
        if a.data.par_iter().any(|v| !v.is_finite()) {
            return Err(LinalgError::NonFinite);
        }
        let mut fwd = vec![0usize; n];
        let mut bwd = vec![0usize; n];
        {
            let p = par();
            let mut mem = MemBuffer::new(factor::lu_in_place_scratch::<usize, f64>(n, n, p, Default::default()));
            let stack = MemStack::new(&mut mem);
            factor::lu_in_place(a.as_mut(), &mut fwd, &mut bwd, p, stack, Default::default());
        }
        let mut dmax = 0.0f64;
        let mut dmin = f64::INFINITY;
        let mut bad = None;
        for i in 0..n {
            let d = a.at(i, i).abs();
            if !(d > 0.0) && bad.is_none() {
                bad = Some(i);
            }
            dmax = dmax.max(d);
            dmin = dmin.min(d);
        }
        if let Some(col) = bad {
            return Err(LinalgError::Singular { col });
        }
        if n > 0 && dmin <= dmax * 1e-14 {
            let col = (0..n).find(|&i| a.at(i, i).abs() <= dmax * 1e-14).unwrap_or(0);
            return Err(LinalgError::Singular { col });
        }
        Ok(Lu { lu: a, fwd, bwd, pivot_ratio: if dmax > 0.0 { dmin / dmax } else { 0.0 } })
    }

    pub fn n(&self) -> usize {
        self.lu.rows
    }

    fn perm(&self) -> faer::perm::PermRef<'_, usize> {
        faer::perm::PermRef::new_checked(&self.fwd, &self.bwd, self.lu.rows)
    }

    /// Solve `A x = b` in place for one right-hand side.
    pub fn solve(&self, b: &mut [f64]) {
        let n = self.lu.rows;
        let rhs = MatMut::from_column_major_slice_mut(b, n, 1);
        self.solve_mat(rhs);
    }

    /// Solve for several right-hand sides at once (columns of `b`).
    pub fn solve_many(&self, b: &mut DenseMatrix) {
        assert_eq!(b.rows, self.lu.rows);
        let r = b.rows;
        let c = b.cols;
        let rhs = MatMut::from_column_major_slice_mut(&mut b.data, r, c);
        self.solve_mat(rhs);
    }

    fn solve_mat(&self, rhs: MatMut<'_, f64>) {
        let n = self.lu.rows;
        let p = par();
        let mut mem = MemBuffer::new(solve::solve_in_place_scratch::<usize, f64>(n, rhs.ncols(), p));
        let stack = MemStack::new(&mut mem);
        let lu = self.lu.as_ref();
        solve::solve_in_place(lu, lu, self.perm(), rhs, p, stack);
    }
}

/// Least-squares solution of an overdetermined system `A X = B` (QR).
/// `a` is rows x cols (rows >= cols); returns cols x B.cols.
pub fn lstsq(a: &DenseMatrix, b: &DenseMatrix) -> Result<DenseMatrix, LinalgError> {
    use faer::linalg::solvers::SolveLstsq;
    assert_eq!(a.rows, b.rows);
    if a.data.iter().any(|v| !v.is_finite()) || b.data.iter().any(|v| !v.is_finite()) {
        return Err(LinalgError::NonFinite);
    }
    let qr = a.as_ref().qr();
    let x = qr.solve_lstsq(b.as_ref());
    let mut out = DenseMatrix::zeros(a.cols, b.cols);
    for j in 0..b.cols {
        for i in 0..a.cols {
            out.set(i, j, x[(i, j)]);
        }
    }
    if out.data.iter().any(|v| !v.is_finite()) {
        return Err(LinalgError::Singular { col: 0 });
    }
    Ok(out)
}

/// Solve a small dense system `A x = b` (A given row-wise).
pub fn solve_dense(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, LinalgError> {
    let lu = Lu::factor(DenseMatrix::from_rows(a))?;
    let mut x = b.to_vec();
    lu.solve(&mut x);
    Ok(x)
}

/// A complex number (eigenvalue / eigenvector component).
#[derive(Clone, Copy, Debug, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
}

/// Eigenvalues and (right) eigenvectors of a real square matrix.
/// `vectors[k]` is the eigenvector of `values[k]`.
pub fn eigen(a: &[Vec<f64>]) -> Result<(Vec<Complex>, Vec<Vec<Complex>>), LinalgError> {
    let m = DenseMatrix::from_rows(a);
    if m.data.iter().any(|v| !v.is_finite()) {
        return Err(LinalgError::NonFinite);
    }
    let n = m.rows;
    let e = m.as_ref().eigen().map_err(|_| LinalgError::NoConvergence)?;
    let s = e.S();
    let u = e.U();
    let mut vals = Vec::with_capacity(n);
    let mut vecs = Vec::with_capacity(n);
    for k in 0..n {
        let z = s.column_vector()[k];
        vals.push(Complex { re: z.re, im: z.im });
        vecs.push((0..n).map(|i| { let c = u[(i, k)]; Complex { re: c.re, im: c.im } }).collect());
    }
    Ok((vals, vecs))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lu_solves() {
        let a = vec![vec![4.0, 1.0, 0.0], vec![1.0, 3.0, 1.0], vec![0.0, 1.0, 2.0]];
        let x = solve_dense(&a, &[1.0, 2.0, 3.0]).unwrap();
        let m = DenseMatrix::from_rows(&a);
        let r = m.mul_vec(&x);
        for (ri, bi) in r.iter().zip([1.0, 2.0, 3.0]) {
            assert!((ri - bi).abs() < 1e-12);
        }
        assert!(matches!(solve_dense(&[vec![1.0, 2.0], vec![2.0, 4.0]], &[1.0, 1.0]), Err(LinalgError::Singular { .. })));
    }
    #[test]
    fn lu_large_parallel() {
        let n = 300;
        let mut m = DenseMatrix::zeros(n, n);
        m.fill_columns_par(|j, c| {
            for (i, v) in c.iter_mut().enumerate() {
                *v = if i == j { n as f64 } else { 1.0 / (1.0 + (i as f64 - j as f64).abs()) };
            }
        });
        let a = m.clone();
        let x0: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let mut b = a.mul_vec(&x0);
        Lu::factor(m).unwrap().solve(&mut b);
        for (u, v) in b.iter().zip(&x0) {
            assert!((u - v).abs() < 1e-10);
        }
    }
    #[test]
    fn eig_rotation() {
        let (v, _) = eigen(&[vec![0.0, -1.0], vec![1.0, 0.0]]).unwrap();
        let mut im: Vec<f64> = v.iter().map(|c| c.im).collect();
        im.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((im[0] + 1.0).abs() < 1e-12 && (im[1] - 1.0).abs() < 1e-12);
    }
}
