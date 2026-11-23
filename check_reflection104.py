"""Check two degree-104 prime exclusions using ordinary Bernoulli numbers.

Run with python-flint 0.9.0 on 64-bit Linux or macOS.
"""

import ctypes
from ctypes.util import find_library
from importlib.metadata import distribution
import sys

import flint


def modular_bernoulli():
    if sys.platform not in ("linux", "darwin") or ctypes.sizeof(ctypes.c_ulong) != 8:
        raise RuntimeError("Use 64-bit Linux or macOS.")
    wheel = distribution("python-flint")
    for entry in wheel.files or []:
        if entry.name.startswith("libflint") and (".so" in entry.name or entry.suffix == ".dylib"):
            library = str(wheel.locate_file(entry))
            break
    else:
        library = find_library("flint")
    if not library:
        raise RuntimeError("Could not locate the FLINT shared library.")
    function = ctypes.CDLL(library).bernoulli_mod_p_harvey
    function.argtypes = (ctypes.c_ulong, ctypes.c_ulong)
    function.restype = ctypes.c_ulong
    return function


def main():
    bernoulli_mod = modular_bernoulli()
    for prime in (101, 1009):
        index = prime - 53
        exact = flint.fmpq.bernoulli(index)
        residue = int(exact.p) * pow(int(exact.q), -1, prime) % prime
        if bernoulli_mod(index, prime) != residue:
            raise ValueError("Exact rational Bernoulli control failed.")

    for prime, expected in (
        (43_279_339_111, 14_311_957_058),
        (490_440_141_803, 449_277_795_842),
    ):
        if not flint.fmpz(prime).is_prime():
            raise ValueError(f"Candidate {prime} is composite.")
        index = prime - 53
        residue = bernoulli_mod(index, prime)
        if residue != expected or residue == 0:
            raise ValueError(f"Bernoulli exclusion failed at {prime}.")
        print(f"p = {prime}, B_{index} mod p = {residue}", flush=True)


if __name__ == "__main__":
    main()
