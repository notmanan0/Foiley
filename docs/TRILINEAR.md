# Linear-doublet triangles

Use `--method trilinear`, `Method::TriLinear`, or `"method": "TriLinear"` in
stored polar definitions. The GUI lists it as linear doublets with Galerkin
testing. Existing `Panel` and `TriPanel` definitions retain their methods.

## Formulation

The solver uses the internal Dirichlet source/doublet representation on thick,
flat triangular panels. There are three **independent** doublet values on each
triangle, in its vertex order, and one constant source strength
`sigma = -normal · Vrel`. Coincident nodes are not welded. This permits the
potential jump between the upper and lower trailing-edge traces.

With barycentric trial and test functions, each target/source pair contributes
a 3 by 3 block. In operator notation:

`(K - M/2 + W T) mu = -S sigma`.

`K` is the principal-value doublet operator, `M` is the P1 mass matrix, `T`
extracts the upper-minus-lower TE trace, and `W` extends it into the wake.
Every row is divided by its target triangle's area. Consequently the self
block is exactly -1/12 on its diagonal and -1/24 off its diagonal. Using a
scalar -1/2 diagonal here would be a different discretization.

Inner source and linear-doublet integrals use solid angles and edge logarithms.
In the smooth far field, degree-five source quadrature avoids cancellation;
it retains the three basis functions rather than replacing them by a monopole.
Outer testing uses a symmetric seven-point degree-five rule. Nearby pairs,
including wake triangles near the target or its ground image, use two uniform
subdivisions (112 points). This is fixed-order
numerical Galerkin integration, not an adaptive error-controlled integrator.

The wake varies linearly across each strip and is constant downstream along
each edge. Its endpoint strengths come from the triangle containing the
**whole TE segment**, with upper and lower vertices matched geometrically.
Neither the opposite vertex nor an average of the two triangles of a quad
enters this condition. Splitting the wake downstream preserves this field.

Six unit solutions cover three translations and three rotations. Ground images
are included in the weak form, including the self-image block. Cancellation is
checked during assembly and before/after factorization; the dense LU call
itself is not interruptible. Triangles are accepted on the mesher's criterion
(area above 1e-14), so slivers from STL/STEP bodies are solved; only non-finite
or zero-area triangles are rejected. If no triangle owns a whole TE segment
(e.g. the mesher dropped a collapsed one), each endpoint takes the nearest
vertex of the remaining TE triangles.

## Loads and output

Surface velocity uses the exact gradient of the linear doublet field within
each triangle. Pressure forces and moments are integrated over that triangle,
including the variation of rotational velocity. Cp output remains one area
average per triangle for the existing display/export interface. Raw doublet
output contains three entries per triangle, and reported unknown counts
include all three.

Lift uses the mean of the two wake endpoint strengths. Trefftz downwash includes
both the endpoint vortices and the distributed streamwise vorticity from the
linear circulation gradient. Strip force/drag reduction uses strip-centre
values. Circulation carried across a fuselage varies linearly between the
inner TE traces of the adjacent strips, so no vortex is shed at the junctions,
and its induced drag uses the same linear Trefftz model (with ground images)
as the strips.
The wake is rigid; the free-wake option remains specific to VLM.

## Validation and limits

Run `cargo test --release -p foiley-plane`; the finest wing mesh and the quad
comparison are an ignored test, run with
`cargo test --release -p foiley-plane --test trilinear -- --ignored`. The
regressions check:

* Analytic kernels against independent surface integration, both one-sided
  jumps, the exact mass matrix, and linear gradients.
* Sphere potential and pressure convergence: from 512 to 2,048 triangles,
  potential RMS error falls from about 0.0040 to 0.0010 and Cp RMS error from
  about 0.087 to 0.042 for a unit sphere and unit freestream.
* An AR=8 rectangular wing with a 1% section gives a lift slope of about
  4.60/rad, consistent with lifting-surface theory.
* The same wing with NACA 0012 at 5 degrees gives CL of 0.4383, 0.4351,
  0.4345 and 0.4342 as the chord/half-span meshes increase from 6/8 through
  10/12 and 14/16 to 24/24. The RMS upper/lower TE-panel Cp difference falls
  from 0.0445 to 0.0083. Integrated pressure lift approaches wake lift;
  their difference is still about 2.9% at the finest tested mesh.
* Refining target integration from 112 to 448 points changes lift by less than
  0.1% on the regression wing. Doubling wake subdivisions leaves lift unchanged
  within 1e-7 relative tolerance.
* A sliver triangle the mesher keeps (1e-13 high) changes lift by less than
  1e-5; a missing TE owner falls back to the nearest vertices.
* Zero-angle and positive/negative-angle symmetry, positive induced drag,
  camber, sweep, ground images, rotational damping, fuselage carry-through,
  and symmetric flaps.

The refined thick-wing CL is about 3.6% above the existing uniform-quad solver's
0.4192. The cross-method regression allows 4%; it does not treat that solver as
an exact reference. These are theory and internal consistency checks, not a
claim of numerical identity with flow5 or validation against measured flight
data. Local pressures and moments converge more slowly than wake lift.

As with other flat-panel methods, refine the surface geometry as well as the
solution. Extremely thin gaps, highly stretched panels, and nearly touching
bodies can require additional integration/refinement beyond these benchmarks.
Three times as many unknowns increases dense matrix storage roughly ninefold
and factorization work roughly 27-fold at equal triangle count. Assembly
writes the matrix in place, without a second dense copy; start with a coarse
mesh and check convergence.
