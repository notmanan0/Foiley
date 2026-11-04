"""F001: the 2D bench harness checks every run and every expected point."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_bench as rb  # noqa: E402

FAKE = os.path.join(HERE, "fake_xfoil.py")
CASE = "n0012_trip"  # ASEQ 0 10 1: 11 points


def fake(mode):
    return [sys.executable, FAKE, mode]


def go(tmp_path, mode, reps=1, timeout=30.0):
    rb.EXES["fake"] = fake(mode)
    return rb.run("fake", CASE, reps, out=str(tmp_path), timeout=timeout)


def f001_harness_rejects_partial_polar(tmp_path):
    ref = go(tmp_path, "full")
    part = go(tmp_path, "partial")
    assert ref.ok and part.ok
    oper = rb.CASES[CASE][1]
    s = rb.compare(ref.rows(), part.rows(), oper)
    assert s.expected == 11 and s.ref_converged == 11
    assert len(s.matched) == 5 and len(s.missing_in_foiley) == 6
    assert not s.ok and "coverage" in s.failures[0]
    # the full polar against itself passes
    assert rb.compare(ref.rows(), ref.rows(), oper).ok


def f001_empty_common_set_fails():
    s = rb.compare([], [], rb.CASES[CASE][1])
    assert not s.ok and "no operating point in common" in s.failures


def f001_nonzero_exit_and_missing_polar_fail(tmp_path):
    crash = go(tmp_path, "crash")
    assert not crash.ok and crash.reasons() == ["exit code 3"]
    nopol = go(tmp_path, "nopolar")
    assert not nopol.ok and nopol.reasons() == ["no polar written"]


def f001_timeout_fails(tmp_path):
    r = go(tmp_path, "hang", timeout=1.0)
    assert not r.ok and r.reasons()[0].startswith("timeout")


def f001_nonfinite_points_scored(tmp_path):
    ref = go(tmp_path, "full")
    bad = go(tmp_path, "nonfinite")
    s = rb.compare(ref.rows(), bad.rows(), rb.CASES[CASE][1])
    assert s.nonfinite == [1.0] and not s.ok


def f001_stale_polar_rejected(tmp_path):
    rb.EXES["fake"] = fake("nopolar")
    pol = tmp_path / f"{CASE}.fake.pol"
    pol.write_text("stale")
    os.utime(pol, (0, 0))
    # run_once deletes the old file first, so a run that writes nothing fails
    s = rb.run_once(fake("nopolar"), rb.script(*rb.CASES[CASE][:1], rb.CASES[CASE][1], pol.name, str(tmp_path)), str(pol), str(tmp_path))
    assert not s.ok and s.reason == "no polar written"


def f001_old_polar_rejected(tmp_path):
    r = go(tmp_path, "stale")
    assert not r.ok and r.reasons() == ["polar older than the run"]


def f001_median_iqr_and_same_repetition(tmp_path):
    r = go(tmp_path, "full", reps=4)
    assert len(r.samples) == 4 and r.ok
    times = sorted(s.seconds for s in r.samples)
    assert r.median() == times[1]  # lower median
    assert r.iqr() >= 0.0
    assert r.median_sample().rows is r.rows()


def f001_expected_points_follow_xfoil_count():
    assert rb.expected_points("ASEQ 0 16 0.5")[1][-1] == 16.0
    assert len(rb.expected_points("ASEQ -5 14 0.25")[1]) == 77
    col, v = rb.expected_points("CSEQ 0 1.2 0.1")
    assert col == 1 and len(v) == 13


if __name__ == "__main__":
    import tempfile, pathlib, inspect
    n = 0
    for name, fn in list(globals().items()):
        if name.startswith("f001_") and callable(fn):
            args = [pathlib.Path(tempfile.mkdtemp())] if inspect.signature(fn).parameters else []
            fn(*args)
            n += 1
            print("ok", name)
    print(f"{n} passed")
