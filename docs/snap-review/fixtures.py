"""Generate small Foiley projects for the viewport regression scripts."""

import argparse
import json
import subprocess
from pathlib import Path


def write(path, data):
    path.write_text(json.dumps(data), encoding="utf-8")


def solve(exe, source, destination, definition, auto_polars=0):
    job = destination.with_suffix(".job.json")
    write(job, {
        "project": str(source),
        "save": str(destination),
        "steps": [{"plane_polar": {"def": definition, "auto_polars": auto_polars}}],
    })
    subprocess.run([str(exe), "--threads", "2", "-q", "run", str(job)], check=True)
    project = json.loads(destination.read_text(encoding="utf-8"))
    if not project["plane_polars"][-1]["points"]:
        raise RuntimeError(f"{definition['name']} has no solved points")
    return project


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("foiley", type=Path, help="path to a development build")
    parser.add_argument("--out", type=Path, default=Path("target/snap/fixtures"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    exe = args.foiley.resolve()
    rectangle = root / "tests/data/projects/review_rect.foiley"
    solve(exe, rectangle, out / "control.foiley", {
        "name": "Tilt sweep", "kind": "Control", "method": "Vlm2", "viscous": False,
        "alpha": 3.0, "range": {"from": 0.0, "to": 1.0, "step": 0.5},
        "ranges": [{"target": {"Tilt": {"wing": 0}}, "min": -12.0, "max": 16.0}],
    })
    rerun = out / "rerun.foiley"
    data = solve(exe, rectangle, rerun, {
        "name": "Same name", "method": "Vlm2", "viscous": False,
        "range": {"from": 6.0, "to": 6.0, "step": 1.0},
    })
    # Model a replaced result with identical pressure and integrated coefficients.
    # The next actual solve restores the CoG at the plane origin.
    data["plane_polars"][-1]["points"][0]["cog"]["y"] = 2.5
    write(rerun, data)
    solve(exe, rectangle, out / "llt.foiley", {
        "name": "LLT stations", "method": "Llt", "viscous": False,
        "range": {"from": 4.0, "to": 4.0, "step": 1.0},
    }, auto_polars=3)


if __name__ == "__main__":
    main()
