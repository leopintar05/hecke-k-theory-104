"""Rebuild the degree-104 matrices and compare their exact determinants."""

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from time import monotonic
from zipfile import ZipFile

import flint
import numpy


CERTIFICATE = "Degree-104-vanishing-certificates.zip"
CERTIFICATE_SHA256 = "8f3829d15aaf8b605c5984c6bd90e160e614b8437ecadb47a7e6c72bdcdd5bc9"
ROOT = Path(__file__).resolve().parent


def run(directory, script, *arguments, capture=False):
    result = subprocess.run(
        [sys.executable, "-u", script, *map(str, arguments)],
        cwd=directory, check=True, text=True, stdout=subprocess.PIPE if capture else None,
    )
    return result.stdout


def require_equal(actual, expected, label):
    if actual != expected:
        raise ValueError(f"Mismatch: {label}")
    print(f"PASS {label}", flush=True)


def first_field_and_classical(root, scratch):
    classical = scratch / "classical"
    with ZipFile(root / "prior/Degree-104-lightweight-checks.zip") as bundle:
        bundle.extractall(classical)
        for script, reference in (
            ("first-field/first_classical.py", "first-field/first.log"),
            ("classical/check.py", "classical/check.log"),
        ):
            output = run(classical, script, capture=True)
            require_equal(output.splitlines(), bundle.read(reference).decode().splitlines(), script)

    first = scratch / "rank324"
    matrices = ("inert/inert-52.matrix.json.gz", "split/split-matrix-52.json.gz")
    results = ("inert/inert-52.json", "split/split-52.json")
    with ZipFile(root / "prior/Degree-104-rank-324-certificates.zip") as bundle:
        # Build into a directory containing neither reference matrix nor either output record.
        bundle.extractall(first, [n for n in bundle.namelist() if n not in matrices + results])
        run(first, "inert/inert104.py", 52, results[0], "--seed", 640813)
        run(first, "split/signed_split.py", 52)
        for name in matrices:
            with bundle.open(name) as stored, gzip.open(stored, "rt") as original:
                with gzip.open(first / name, "rt") as rebuilt:
                    require_equal(json.load(rebuilt), json.load(original), f"all entries of {name}")
    run(first, "verify_support.py")
    run(first, "verify_matrices.py", "--exact")


def replay(scratch, workers, threads):
    start = monotonic()
    print("Start UTC", datetime.now(timezone.utc).isoformat(timespec="seconds"), flush=True)
    print(f"Python {platform.python_version()}, python-flint {flint.__version__}, "
          f"NumPy {numpy.__version__}, {platform.system()} {platform.machine()}", flush=True)
    certificate = ROOT / "certificates" / CERTIFICATE
    require_equal(hashlib.sha256(certificate.read_bytes()).hexdigest(), CERTIFICATE_SHA256,
                  "degree-104 source archive")
    root = scratch / "degree104"
    with ZipFile(certificate) as bundle:
        bundle.extractall(root)
    run(ROOT, "check_adjoint104.py")
    first_field_and_classical(root, scratch)

    field = root / "second-field"
    matrix = field / "generated/rows104.bin"
    require_equal(matrix.exists(), False, "second field starts without a saved matrix")
    run(root, "second-field/run104.py", "--workers", workers)
    run(root, "second-field/prime_det.py", matrix, field / "config104.json",
        "recomputed.json", "--threads", threads)
    require_equal((root / "recomputed.json").read_bytes(), (root / "degree104-result.json").read_bytes(),
                  "four complete determinants and CRT record (byte for byte)")
    run(root, "check_chain.py", "--matrix", matrix)
    rows = run(root, "audit/audit_rows.py", "audit/archived-row-source.zip",
               "second-field/config104.json", "--packed", matrix, capture=True)
    require_equal(rows.splitlines(), (root / "audit/direct-row-digests.txt").read_text().splitlines(),
                  "three direct row checks (8,775 entries)")
    for line in rows.splitlines():
        print("Direct row", line)
    result = json.loads((root / "recomputed.json").read_text())
    print(f"PASS degree 104: rank {result['rank']}, four prime fields, gcd {result['gcd']}.")
    print("End UTC", datetime.now(timezone.utc).isoformat(timespec="seconds"),
          f"({monotonic() - start:.1f} seconds)", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, help="use an empty directory and retain reconstructed files")
    parser.add_argument("--workers", type=int, default=5, help="matrix generation workers")
    parser.add_argument("--threads", type=int, default=3, help="FLINT determinant threads")
    args = parser.parse_args()
    if not __debug__ or sys.flags.optimize or min(args.workers, args.threads) < 1:
        parser.error("Use an unoptimized Python interpreter and positive worker counts.")
    if args.work_dir:
        scratch = args.work_dir.resolve()
        scratch.mkdir(parents=True, exist_ok=True)
        if any(scratch.iterdir()):
            parser.error("--work-dir must be empty")
        replay(scratch, args.workers, args.threads)
    else:
        with tempfile.TemporaryDirectory(prefix="hecke-104-") as directory:
            replay(Path(directory), args.workers, args.threads)


if __name__ == "__main__":
    main()
