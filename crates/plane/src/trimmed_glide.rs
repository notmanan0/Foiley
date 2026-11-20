//! Longitudinal glide equilibrium at prescribed airspeeds.

use crate::analysis::{
    envelope_failure, evaluate, make_system, mesh_options, nonfinite_reason, wind_dir, Ctx,
    ElevatorTrim, PlaneOpp, PlanePolar, G,
};
use crate::evaluator::Evaluator;
use crate::mesh::{build, Mesh};
use crate::PlaneError;
use foiley_core::linalg::Diagnostics;
use foiley_core::{Progress, Vec3};
use std::collections::HashMap;

const FORCE_TOL: f64 = 1e-6;
const CM_TOL: f64 = 1e-7;
const ALPHA_LIMIT: f64 = 30.0;

struct Sample {
    x: f64,
    residual: f64,
    point: PlaneOpp,
}

fn same_sign(a: &Sample, b: &Sample) -> bool {
    a.residual.is_sign_positive() == b.residual.is_sign_positive()
}

fn root(
    mut lo: Sample,
    mut hi: Sample,
    tolerance: f64,
    mut sample: impl FnMut(f64) -> Result<Option<Sample>, PlaneError>,
) -> Result<Option<Sample>, PlaneError> {
    if lo.residual.abs() <= tolerance {
        return Ok(Some(lo));
    }
    if hi.residual.abs() <= tolerance {
        return Ok(Some(hi));
    }
    if same_sign(&lo, &hi) {
        return Ok(None);
    }

    for i in 0..60 {
        let x = lo.x - lo.residual * (hi.x - lo.x) / (hi.residual - lo.residual);
        // Periodic bisection preserves progress when the secant repeatedly
        // lands at one end, for example near a 2D polar interpolation knot.
        let x = if i % 3 != 2 && x > lo.x && x < hi.x {
            x
        } else {
            (lo.x + hi.x) * 0.5
        };
        let Some(mid) = sample(x)? else {
            return Ok(None);
        };
        if mid.residual.abs() <= tolerance {
            return Ok(Some(mid));
        }
        if same_sign(&lo, &mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(None)
}

fn mesh_at(ctx: &Ctx<'_>, trim: ElevatorTrim, deg: f64) -> Result<Mesh, PlaneError> {
    let mut options = mesh_options(ctx.def, &HashMap::new(), &HashMap::new());
    options.elevator = Some((trim.wing, deg));
    build(ctx.plane, ctx.db, &options)
}

/// At a fixed elevator angle and airspeed, find the alpha where the
/// aerodynamic resultant balances weight, then return Cm as the outer residual.
fn pitch_sample(
    ctx: &Ctx<'_>,
    ev: &Evaluator<'_>,
    speed: f64,
    deg: f64,
    prog: &Progress,
) -> Result<Option<Sample>, PlaneError> {
    let mass = ctx.inertia.mass;
    let weight = mass * G;
    let cog = ctx.inertia.cog;
    let force_sample = |alpha| {
        prog.check()?;
        let (raw, _) = ev.loads(wind_dir(alpha, 0.0), speed, Vec3::ZERO, prog)?;
        let point = ctx.reduce(ev.mesh(), &raw, cog, alpha, 0.0, speed, mass, None);
        // Signed resultant keeps negative-lift states on the low side of
        // the bracket. Accepted states must have positive lift and drag.
        let residual = point.lift.hypot(point.drag).copysign(point.lift) / weight - 1.0;
        Ok(residual.is_finite().then_some(Sample {
            x: alpha,
            residual,
            point,
        }))
    };
    let (Some(lo), Some(hi)) = (force_sample(-ALPHA_LIMIT)?, force_sample(ALPHA_LIMIT)?) else {
        return Ok(None);
    };
    let Some(balance) = root(lo, hi, FORCE_TOL, force_sample)? else {
        return Ok(None);
    };
    let mut point = balance.point;
    if point.lift <= 0.0 || point.drag <= 0.0 || !point.cm.is_finite() {
        return Ok(None);
    }
    point.elevator_deg = Some(deg);
    point.balance_residual = Some(balance.residual.abs());
    Ok(Some(Sample {
        x: deg,
        residual: point.cm,
        point,
    }))
}

pub(crate) fn run(
    ctx: &Ctx<'_>,
    polar: &mut PlanePolar,
    prog: &Progress,
    on_point: &(dyn Fn(&PlaneOpp) + Sync),
) -> Result<(), PlaneError> {
    let trim = ctx
        .def
        .elevator_trim
        .expect("trimmed glide definition was validated");
    let stations = trim.stations();
    let first_mesh = mesh_at(ctx, trim, stations[0])?;
    trim.validate_mesh(ctx.plane, &first_mesh)?;
    polar.unknowns = first_mesh.panels.len();

    // The rigid VLM2 systems at the bracket stations are independent of
    // speed. Factor them once and reuse them for the full speed sweep.
    let mut meshes = vec![first_mesh];
    prog.stage("Preparing elevator stations");
    for &deg in &stations[1..] {
        prog.check()?;
        meshes.push(mesh_at(ctx, trim, deg)?);
    }
    let prepare = prog.sub(0.0, 0.15);
    let systems = meshes
        .iter()
        .enumerate()
        .map(|(i, mesh)| {
            make_system(
                ctx.def,
                mesh,
                ctx.inertia.cog,
                ctx.c_ref,
                &prepare.part(i, meshes.len()),
            )
            .map(|sys| Evaluator::new(sys, ctx.def, ctx.c_ref))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let speeds = ctx.def.range.values();
    let sweep = prog.sub(0.15, 1.0);
    for (i, &speed) in speeds.iter().enumerate() {
        let sub = sweep.part(i, speeds.len());
        sub.check()?;
        sub.stage(format!(
            "Trimming at {speed:.2} m/s ({}/{})",
            i + 1,
            speeds.len()
        ));
        let solved = solve_speed(ctx, trim, speed, &stations, &systems, &sub)?;
        let result = match solved {
            Some(point) => match nonfinite_reason(&point).or_else(|| envelope_failure(ctx.def, &point)) {
                Some(reason) => Err(reason),
                None => Ok(point),
            },
            None => Err(format!("no glide trim found within elevator travel [{:.2}, {:.2}] deg and alpha ±{ALPHA_LIMIT} deg", trim.min_deg, trim.max_deg)),
        };
        match result {
            Ok(point) => {
                prog.count_ok();
                on_point(&point);
                polar.points.push(point);
            }
            Err(reason) => {
                prog.count_failed();
                polar.failed.push(format!("V={speed:.3} m/s: {reason}"));
            }
        }
        sub.set(1.0);
    }
    Ok(())
}

fn solve_speed(
    ctx: &Ctx<'_>,
    trim: ElevatorTrim,
    speed: f64,
    stations: &[f64],
    systems: &[Evaluator<'_>],
    prog: &Progress,
) -> Result<Option<PlaneOpp>, PlaneError> {
    let at = |deg| {
        prog.check()?;
        let mesh = mesh_at(ctx, trim, deg)?;
        let ev = Evaluator::new(
            make_system(
                ctx.def,
                &mesh,
                ctx.inertia.cog,
                ctx.c_ref,
                &prog.sub(0.0, 0.0),
            )?,
            ctx.def,
            ctx.c_ref,
        );
        pitch_sample(ctx, &ev, speed, deg, prog)
    };
    let mut previous = None;
    let mut solution = None;
    for (&deg, ev) in stations.iter().zip(systems) {
        let Some(sample) = pitch_sample(ctx, ev, speed, deg, prog)? else {
            previous = None;
            continue;
        };
        if sample.residual.abs() <= CM_TOL {
            solution = Some(sample);
            break;
        }
        if let Some(lo) = previous.take().filter(|lo| !same_sign(lo, &sample)) {
            solution = root(lo, sample, CM_TOL, at)?;
            if solution.is_some() {
                break;
            }
        } else {
            previous = Some(sample);
        }
    }
    let Some(Sample {
        x: deg, mut point, ..
    }) = solution
    else {
        return Ok(None);
    };
    if ctx.def.diagnostics != Diagnostics::Off {
        let mesh = mesh_at(ctx, trim, deg)?;
        let ev = Evaluator::new(
            make_system(
                ctx.def,
                &mesh,
                ctx.inertia.cog,
                ctx.c_ref,
                &prog.sub(0.0, 0.0),
            )?,
            ctx.def,
            ctx.c_ref,
        );
        point.diagnostics = evaluate(
            ctx,
            &ev,
            ctx.inertia.cog,
            point.alpha,
            0.0,
            speed,
            ctx.inertia.mass,
            prog,
        )?
        .diagnostics;
    }
    Ok(Some(point))
}
