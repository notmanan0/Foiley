//! Dense linear algebra on top of `faer`.
//!
//! Influence matrices are stored column-major so that each column (the
//! influence of one singularity element on every collocation point) can be
//! filled by an independent rayon task, then factored in place by faer's
//! blocked, multithreaded partial-pivoting LU without an extra copy.

use crate::progress::Progress;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::lu::partial_pivoting::{factor, solve};
use faer::reborrow::{Reborrow, ReborrowMut};
use faer::{Accum, MatMut, MatRef, Par};
use rayon::prelude::*;

#[derive(Debug, Clone, PartialEq)]
pub enum LinalgError {
    /// The matrix is (numerically) singular; `col` is the first bad pivot.
    Singular { col: usize },
    /// The matrix contains NaN or infinite entries.
    NonFinite,
    /// Eigen solver did not converge.
    NoConvergence,
    /// The factorization was cancelled through its `Progress`.
    Cancelled,
}

impl std::fmt::Display for LinalgError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LinalgError::Singular { col } => write!(f, "matrix is singular (pivot {col} is zero): check for duplicated or zero-area panels"),
            LinalgError::NonFinite => write!(f, "matrix contains NaN/inf values: check the geometry for degenerate panels"),
            LinalgError::NoConvergence => write!(f, "eigenvalue solver did not converge"),
            LinalgError::Cancelled => write!(f, "cancelled"),
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
    /// min|u_ii| / max|u_ii| over the diagonal of U.  This is only a crude
    /// warning sign: it is neither a pivot-growth factor nor a condition
    /// number (use `cond1_estimate`).  Formerly `pivot_ratio`.
    pub diag_ratio: f64,
    /// ||A||_1 of the matrix before factorisation (max column sum).
    pub norm1: f64,
}

/// Below this size `Lu::factor_with` uses the single faer call.
pub const BLOCKED_LU_MIN: usize = 2000;

static BLOCKED_MIN: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(BLOCKED_LU_MIN);

/// Override the size from which `Lu::factor_with` uses the blocked LU
/// (process-wide; for tests and benchmarks: `usize::MAX` = never).
pub fn set_blocked_lu_min(n: usize) {
    BLOCKED_MIN.store(n, std::sync::atomic::Ordering::Relaxed);
}
/// Column panel width of the blocked LU.
pub const BLOCKED_LU_PANEL: usize = 128;
/// Work per uninterruptible Schur-update tile, in flops (sets the tile
/// width from the trailing height, so the longest step stays a few tens of
/// milliseconds on a desktop CPU).
pub const BLOCKED_LU_TILE_FLOPS: f64 = 4e9;

/// Transpositions (i, j), i < j, that turn the identity into `target` when
/// applied in order to an array (so `a[i]` becomes `a[target[i]]`).
fn transpositions(target: &[usize]) -> Vec<(usize, usize)> {
    let m = target.len();
    let mut cur: Vec<usize> = (0..m).collect();
    let mut pos: Vec<usize> = (0..m).collect();
    let mut out = Vec::new();
    for i in 0..m {
        let v = target[i];
        let p = pos[v];
        if p != i {
            let u = cur[i];
            cur.swap(i, p);
            pos[u] = p;
            pos[v] = i;
            out.push((i, p));
        }
    }
    out
}

/// Right-looking blocked LU with partial pivoting, cancellable between
/// steps.  For each column panel: faer LU of the tall panel, the panel's row
/// swaps on the other columns (column-parallel), then for each tile of
/// trailing columns a unit-lower triangular solve for U12 and the Schur
/// update A22 -= L21 U12 as one parallel faer matmul, with a cancellation
/// check between tiles.  Leaves P A = L U in `a`; returns the row order
/// (`perm[i]` = original row now at position i).
fn blocked_lu(a: &mut DenseMatrix, nb: usize, tile_flops: f64, prog: &Progress) -> Result<Vec<usize>, LinalgError> {
    let n = a.rows;
    let p = par();
    let mut perm: Vec<usize> = (0..n).collect();
    let (mut pf, mut pb) = (vec![0usize; n], vec![0usize; n]);
    let mut k = 0;
    while k < n {
        if prog.is_cancelled() {
            return Err(LinalgError::Cancelled);
        }
        let b = nb.min(n - k);
        let m = n - k;
        {
            let mut mem = MemBuffer::new(factor::lu_in_place_scratch::<usize, f64>(m, b, p, Default::default()));
            let panel = a.as_mut().submatrix_mut(k, k, m, b);
            factor::lu_in_place(panel, &mut pf[..m], &mut pb[..m], p, MemStack::new(&mut mem), Default::default());
        }
        let swaps = transpositions(&pf[..m]);
        if !swaps.is_empty() {
            a.data.par_chunks_mut(n).enumerate().filter(|(j, _)| *j < k || *j >= k + b).for_each(|(_, col)| {
                for &(i, j) in &swaps {
                    col.swap(k + i, k + j);
                }
            });
            for &(i, j) in &swaps {
                perm.swap(k + i, k + j);
            }
        }
        let rest = n - k - b;
        if rest > 0 {
            let t = ((tile_flops / (2.0 * m as f64 * b as f64)) as usize).max(64).min(rest);
            let mut c0 = k + b;
            while c0 < n {
                if prog.is_cancelled() {
                    return Err(LinalgError::Cancelled);
                }
                let w = t.min(n - c0);
                let am = a.as_mut();
                let (left, right) = am.split_at_col_mut(c0);
                let lp = left.rb().submatrix(k, k, m, b);
                let (l11, l21) = (lp.submatrix(0, 0, b, b), lp.submatrix(b, 0, m - b, b));
                let tile = right.submatrix_mut(k, 0, m, w);
                let (mut a12, a22) = tile.split_at_row_mut(b);
                faer::linalg::triangular_solve::solve_unit_lower_triangular_in_place(l11, a12.rb_mut(), p);
                faer::linalg::matmul::matmul(a22, Accum::Add, l21, a12.rb(), -1.0, p);
                c0 += w;
            }
        }
        k += b;
    }
    Ok(perm)
}

impl Lu {
    /// Factor a square matrix, consuming it.
    pub fn factor(a: DenseMatrix) -> Result<Lu, LinalgError> {
        Self::factor_with(a, &Progress::none())
    }

    /// Factor, cancellable through `prog` (blocked above `BLOCKED_LU_MIN`
    /// unknowns; below it the single faer call is fast enough).
    pub fn factor_with(a: DenseMatrix, prog: &Progress) -> Result<Lu, LinalgError> {
        let blocked = a.rows >= BLOCKED_MIN.load(std::sync::atomic::Ordering::Relaxed);
        Self::factor_impl(a, prog, blocked, BLOCKED_LU_PANEL, BLOCKED_LU_TILE_FLOPS)
    }

    /// Always use the blocked algorithm with these block sizes (tests, tuning).
    pub fn factor_blocked(a: DenseMatrix, prog: &Progress, nb: usize, tile_flops: f64) -> Result<Lu, LinalgError> {
        Self::factor_impl(a, prog, true, nb.max(1), tile_flops)
    }

    fn factor_impl(mut a: DenseMatrix, prog: &Progress, blocked: bool, nb: usize, tile_flops: f64) -> Result<Lu, LinalgError> {
        assert_eq!(a.rows, a.cols, "LU needs a square matrix");
        let n = a.rows;
        // one O(N^2) pass: finiteness and the 1-norm (max column sum)
        let colsum = |c: &[f64]| c.iter().try_fold(0.0f64, |s, v| if v.is_finite() { Some(s + v.abs()) } else { None });
        let norm1 = match a.data.par_chunks(n.max(1)).map(colsum).try_reduce(|| 0.0, |x, y| Some(x.max(y))) {
            Some(v) => v,
            None => return Err(LinalgError::NonFinite),
        };
        let mut fwd = vec![0usize; n];
        let mut bwd = vec![0usize; n];
        if blocked {
            fwd = match blocked_lu(&mut a, nb, tile_flops, prog) {
                Ok(f) => f,
                Err(e) => {
                    // Releasing a multi-GB buffer can take tens of ms (page
                    // decommit); do it off the cancelling thread.
                    let data = std::mem::take(&mut a.data);
                    if data.len() > 1 << 24 {
                        std::thread::spawn(move || drop(data));
                    }
                    return Err(e);
                }
            };
            for (i, &f) in fwd.iter().enumerate() {
                bwd[f] = i;
            }
        } else {
            if prog.is_cancelled() {
                return Err(LinalgError::Cancelled);
            }
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
        Ok(Lu { lu: a, fwd, bwd, diag_ratio: if dmax > 0.0 { dmin / dmax } else { 0.0 }, norm1 })
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

    /// Solve `A^T x = b` in place for one right-hand side.
    pub fn solve_transpose(&self, b: &mut [f64]) {
        let n = self.lu.rows;
        let rhs = MatMut::from_column_major_slice_mut(b, n, 1);
        let p = par();
        let mut mem = MemBuffer::new(solve::solve_transpose_in_place_scratch::<usize, f64>(n, 1, p));
        let stack = MemStack::new(&mut mem);
        let lu = self.lu.as_ref();
        solve::solve_transpose_in_place(lu, lu, self.perm(), rhs, p, stack);
    }

    /// Hager-Higham estimate of ||A^-1||_1 (LAPACK xLACN2 iteration plus
    /// Higham's alternative vector): at most 5 solves with A and 4 with A^T
    /// using the existing factors, O(N^2).  A lower bound, usually within a
    /// small factor of the exact value.
    pub fn inv_norm1_estimate(&self) -> f64 {
        let n = self.n();
        if n == 0 {
            return 0.0;
        }
        let norm1 = |v: &[f64]| v.iter().map(|x| x.abs()).sum::<f64>();
        let sign = |v: &[f64]| v.iter().map(|x| if *x >= 0.0 { 1.0 } else { -1.0 }).collect::<Vec<f64>>();
        let argmax = |v: &[f64]| (0..v.len()).max_by(|&a, &b| v[a].abs().total_cmp(&v[b].abs())).unwrap_or(0);
        let mut x = vec![1.0 / n as f64; n];
        self.solve(&mut x);
        let mut est = norm1(&x);
        if n > 1 {
            let mut xi = sign(&x);
            let mut z = xi.clone();
            self.solve_transpose(&mut z);
            let mut j = argmax(&z);
            for _ in 0..4 {
                let mut y = vec![0.0; n];
                y[j] = 1.0;
                self.solve(&mut y);
                let old = est;
                est = norm1(&y);
                let s = sign(&y);
                if s == xi || est <= old {
                    est = est.max(old);
                    break;
                }
                xi = s;
                z.clone_from(&xi);
                self.solve_transpose(&mut z);
                let jl = j;
                j = argmax(&z);
                if z[jl].abs() == z[j].abs() {
                    break;
                }
            }
        }
        // alternative estimate: x_i = (-1)^i (1 + i/(n-1))
        let mut alt: Vec<f64> = (0..n).map(|i| (if i % 2 == 0 { 1.0 } else { -1.0 }) * (1.0 + i as f64 / (n.max(2) - 1) as f64)).collect();
        self.solve(&mut alt);
        est.max(2.0 * norm1(&alt) / (3.0 * n as f64))
    }

    /// Estimated 1-norm condition number ||A||_1 ||A^-1||_1.
    pub fn cond1_estimate(&self) -> f64 {
        self.norm1 * self.inv_norm1_estimate()
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

/// A square matrix whose rows can be re-assembled on demand in f64, in
/// blocks of `block_size()` consecutive rows.  Used for residuals that must
/// not come from a rounded or factored copy of the operator.
pub trait RowOperator: Sync {
    fn dim(&self) -> usize;
    fn block_size(&self) -> usize {
        1
    }
    /// Fill block `b` (rows `b*bs .. b*bs+bs`) row-major: `out[r*n + j] = A[b*bs + r][j]`.
    fn fill_block(&self, b: usize, out: &mut [f64]);
}

impl RowOperator for DenseMatrix {
    fn dim(&self) -> usize {
        self.rows
    }
    fn fill_block(&self, b: usize, out: &mut [f64]) {
        for (j, o) in out.iter_mut().enumerate() {
            *o = self.at(b, j);
        }
    }
}

/// A `RowOperator` from a closure.
pub struct FnRows<F: Fn(usize, &mut [f64]) + Sync> {
    pub n: usize,
    pub block: usize,
    pub fill: F,
}

impl<F: Fn(usize, &mut [f64]) + Sync> RowOperator for FnRows<F> {
    fn dim(&self) -> usize {
        self.n
    }
    fn block_size(&self) -> usize {
        self.block
    }
    fn fill_block(&self, b: usize, out: &mut [f64]) {
        (self.fill)(b, out)
    }
}

/// How much checking of a linear solve to do (off by default).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Diagnostics {
    #[default]
    Off,
    /// condition estimate plus a residual on sampled rows: O(N^2)
    Cheap,
    /// exact residual on every row (operator re-assembled in f64) and one
    /// step of iterative refinement, kept only if it lowers the residual
    Full,
}

/// Residuals of the solutions `x` (columns) on some rows of `op`.
#[derive(Clone, Debug, Default)]
pub struct RowResiduals {
    /// rows checked, ascending
    pub rows: Vec<usize>,
    /// max abs row sum of `op` over those rows (||A||_inf when all rows)
    pub a_inf: f64,
    /// per right-hand side: b at those rows, and r = b - A x there
    pub b: Vec<Vec<f64>>,
    pub r: Vec<Vec<f64>>,
}

/// Re-assemble blocks `blocks` of `op` (in parallel) and evaluate the
/// residual of every column of `x` against the same column of `b`.
pub fn row_residuals(op: &dyn RowOperator, blocks: &[usize], x: &DenseMatrix, b: &DenseMatrix) -> RowResiduals {
    let n = op.dim();
    let bs = op.block_size();
    let k = x.cols;
    let per: Vec<(usize, f64, Vec<f64>)> = blocks
        .par_iter()
        .flat_map_iter(|&blk| {
            let mut a = vec![0.0; bs * n];
            op.fill_block(blk, &mut a);
            (0..bs)
                .filter(move |r| blk * bs + r < n)
                .map(|r| {
                    let i = blk * bs + r;
                    let row = &a[r * n..(r + 1) * n];
                    let asum: f64 = row.iter().map(|v| v.abs()).sum();
                    let res: Vec<f64> = (0..k).map(|c| b.at(i, c) - row.iter().zip(x.col(c)).map(|(a, x)| a * x).sum::<f64>()).collect();
                    (i, asum, res)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let mut per = per;
    per.sort_by_key(|p| p.0);
    let rows: Vec<usize> = per.iter().map(|p| p.0).collect();
    RowResiduals {
        a_inf: per.iter().map(|p| p.1).fold(0.0, f64::max),
        b: (0..k).map(|c| rows.iter().map(|&i| b.at(i, c)).collect()).collect(),
        r: (0..k).map(|c| per.iter().map(|p| p.2[c]).collect()).collect(),
        rows,
    }
}

/// Normwise backward error ||r||_inf / (||A||_inf ||x||_inf + ||b||_inf).
pub fn backward_error(r_inf: f64, a_inf: f64, x_inf: f64, b_inf: f64) -> f64 {
    let d = a_inf * x_inf + b_inf;
    if d > 0.0 {
        r_inf / d
    } else {
        0.0
    }
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, x| m.max(x.abs()))
}

/// What the diagnostics found for a factored system and its solutions.
#[derive(Clone, Debug, Default)]
pub struct SolveReport {
    pub level: Diagnostics,
    pub n: usize,
    pub diag_ratio: f64,
    pub cond1_estimate: Option<f64>,
    pub residuals: RowResiduals,
    /// per right-hand side: refinement step kept
    pub refined: Vec<bool>,
}

impl SolveReport {
    /// Backward error of the combination x = sum c_k x_k (b = sum c_k b_k)
    /// on the checked rows; `x_inf` is ||x||_inf of the combined solution.
    pub fn backward_error(&self, c: &[f64], x_inf: f64) -> Option<f64> {
        let res = &self.residuals;
        if res.rows.is_empty() {
            return None;
        }
        let comb = |v: &[Vec<f64>]| (0..res.rows.len()).map(|i| c.iter().zip(v).map(|(c, col)| c * col[i]).sum::<f64>()).collect::<Vec<f64>>();
        let (r, b) = (comb(&res.r), comb(&res.b));
        Some(backward_error(inf_norm(&r), res.a_inf, x_inf, inf_norm(&b)))
    }
    /// ||r||_inf of the combination on the checked rows.
    pub fn residual_inf(&self, c: &[f64]) -> Option<f64> {
        let res = &self.residuals;
        if res.rows.is_empty() {
            return None;
        }
        Some((0..res.rows.len()).map(|i| c.iter().zip(&res.r).map(|(c, col)| c * col[i]).sum::<f64>().abs()).fold(0.0, f64::max))
    }
}

/// Run the requested diagnostics on `lu` (the factors of `op`), the
/// right-hand sides `b` and their solutions `x` (refined in place for
/// `Full` when that lowers a column's residual).  `must_rows` are rows that
/// are always checked in `Cheap` mode (trailing-edge / Kutta rows).
pub fn diagnose(level: Diagnostics, lu: &Lu, op: &dyn RowOperator, b: &DenseMatrix, x: &mut DenseMatrix, must_rows: &[usize], samples: usize) -> SolveReport {
    let n = op.dim();
    let bs = op.block_size().max(1);
    let nb = n.div_ceil(bs);
    let mut rep = SolveReport { level, n, diag_ratio: lu.diag_ratio, ..Default::default() };
    if level == Diagnostics::Off || n == 0 {
        return rep;
    }
    rep.cond1_estimate = Some(lu.cond1_estimate());
    let blocks: Vec<usize> = match level {
        Diagnostics::Full => (0..nb).collect(),
        _ => {
            let mut v: Vec<usize> = must_rows.iter().map(|r| r / bs).collect();
            let rows = samples.min(n);
            v.extend((0..rows).map(|k| (k * n / rows.max(1)) / bs));
            v.sort_unstable();
            v.dedup();
            v
        }
    };
    let mut res = row_residuals(op, &blocks, x, b);
    rep.refined = vec![false; x.cols];
    if level == Diagnostics::Full {
        // one step of iterative refinement, in f64, against the re-assembled operator
        let mut d = DenseMatrix::zeros(n, x.cols);
        for c in 0..x.cols {
            for (k, &i) in res.rows.iter().enumerate() {
                d.set(i, c, res.r[c][k]);
            }
        }
        lu.solve_many(&mut d);
        let mut x2 = x.clone();
        for (v, dv) in x2.data.iter_mut().zip(&d.data) {
            *v += dv;
        }
        let res2 = row_residuals(op, &blocks, &x2, b);
        for c in 0..x.cols {
            if inf_norm(&res2.r[c]) < inf_norm(&res.r[c]) {
                x.col_mut(c).copy_from_slice(x2.col(c));
                res.r[c].clone_from(&res2.r[c]);
                rep.refined[c] = true;
            }
        }
    }
    rep.residuals = res;
    rep
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
