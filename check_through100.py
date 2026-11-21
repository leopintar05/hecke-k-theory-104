"""Check the archived integers for the degree-100 argument with standard Python."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, gcd, isqrt, prod
from pathlib import Path
import re
from zipfile import ZipFile

ARCHIVE_HASH = "c1a718a1d6d5b0e268c91a28fcfebfc2cc462cb7cd058cde6db470a24243ccc2"


def require(test, message):
    if not test:
        raise ValueError(message)


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def factor_product(factors, upper):
    require(len({p for p, exponent in factors}) == len(factors), "duplicate factor")
    for p, exponent in factors:
        require(2 <= p <= upper and exponent > 0 and trial_prime(p), f"factor {p}^{exponent}")
    return prod(p**exponent for p, exponent in factors)


def divided_bernoulli(maximum):
    # For eta of conductor three, multiplying its generating series by
    # exp(3t)-1 gives t(exp(t)-exp(2t)). Equate coefficients at degree n+1.
    bernoulli = []
    for n in range(maximum + 1):
        lower = sum(comb(n + 1, i) * 3**(n + 1 - i) * v for i, v in enumerate(bernoulli))
        bernoulli.append(Fraction((n + 1) * (1 - 2**n) - lower, 3 * (n + 1)))
    return {j: abs((bernoulli[j] / j).numerator) for j in range(1, maximum + 1)}


def verify_manifest(bundle):
    manifest = "SHA256SUMS"
    entries = bundle.read(manifest).decode().splitlines()
    files = []
    for line in entries:
        expected, name = line.split(None, 1)
        files.append(name)
        require(sha256(bundle.read(name)).hexdigest() == expected, name)
    require(sorted(files + [manifest]) == sorted(bundle.namelist()), "incomplete manifest")
    return len(files)


def check_first_field(support, b):
    require([r["e"] for r in support] == list(range(2, 51, 2)), "first-field degrees")
    for r in support:
        e = r["e"]
        ratio = Fraction(int(r["G"]), b[e + 1] * b[2 * e + 1]**3)
        require(ratio.numerator == factor_product(r["ratio_num"], 337), f"ratio numerator, e={e}")
        require(ratio.denominator == factor_product(r["ratio_den"], 337), f"ratio denominator, e={e}")


def check_block_products(bundle, e, cert):
    log = bundle.read(f"logs/panel-{e}.log").decode()
    blocks = re.findall(r"^pivot verified (\d+) of (\d+) block determinant (\d+) seconds ", log, re.M)
    rank, width, modulus = cert["rank"], cert["block"], cert["modulus"]
    require([int(n) for n, _, _ in blocks] == list(range(width, rank, width)) + [rank], "incomplete blocks")
    require(all(int(r) == rank for _, r, _ in blocks), "inconsistent block dimensions")
    determinant = prod(int(v) for _, _, v in blocks) % modulus
    require(determinant == cert["determinant_mod_M"], f"block product, e={e}")
    require(gcd(determinant, modulus) == cert["gcd"] == 1, f"determinant gcd, e={e}")
    require(f"packed SHA-256 {cert['packed_sha256']}" in log.splitlines(), f"logged matrix hash, e={e}")
    return len(blocks)


def check_upper_rows(bundle, b):
    weights = json.loads(bundle.read("weights.json"))
    certs = json.loads(bundle.read("certificates.json"))
    require(set(weights) == set(certs) == {str(e) for e in range(38, 51, 2)}, "upper degrees")
    for e in range(38, 51, 2):
        w, c = weights[str(e)], certs[str(e)]
        require(2 * e == w["degree"] == c["degree"], f"degree for e={e}")
        removed = [(int(p), exponent) for p, exponent in w["removed"].items()]
        small = factor_product(removed, 2**31 - 1)
        numerators = [b[e + 1], b[2 * e + 1]]
        require(numerators == w["b"] and prod(numerators) == small * c["modulus"], f"modulus for e={e}")
        require(c["modulus"] == w["modulus"] and c["modulus"].bit_length() == w["bits"], "modulus record")
        require(c["rank"] == w["rank"] == comb(e // 2 + 1, 3), "signed rank")
        blocks = check_block_products(bundle, e, c)
        print(f"degree {2*e:3}   rank {c['rank']}   modulus {w['bits']} bits   block product ({blocks:2} blocks)   gcd 1")
    return certs["50"]


def degree_100_crosscheck(bundle, cert):
    log = bundle.read("logs/prime-50.log").decode()
    record = re.search(r"^direct full determinant e rank prime residue (\d+) (\d+) (\d+) (\d+) ", log, re.M)
    require(record is not None, "missing single-prime cross-check")
    e, rank, prime, residue = map(int, record.groups())
    require((e, rank) == (50, cert["rank"]), "single-prime dimensions")
    require(prime == 109534107737962267421 and cert["modulus"] % prime == 0, "degree-100 prime divisor")
    # Lucas primality certificate. These checks give 3 order p-1 modulo p.
    factors = [(2, 2), (5, 1), (23, 1), (47, 1), (6990751, 1), (724719341, 1)]
    require(factor_product(factors, 724719341) == prime - 1, "factorization of p-1")
    require(pow(3, prime - 1, prime) == 1, "Fermat condition")
    require(all(pow(3, (prime - 1) // q, prime) != 1 for q, _ in factors), "order of 3")
    require(cert["determinant_mod_M"] % prime == residue == 14999523935004074291, "single-prime match")
    print(f"Degree 100: M = {cert['modulus']}")
    print(f"Degree 100: determinant modulo M = {cert['determinant_mod_M']}")
    print(f"Degree 100: verified prime {prime}, residue {residue} matches the saved direct determinant.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="K-theory-through-100-certificates.zip")
    args = parser.parse_args()
    require(sha256(args.archive.read_bytes()).hexdigest() == ARCHIVE_HASH, "original archive SHA-256")
    with ZipFile(args.archive) as bundle:
        print(f"Original archive SHA-256 {ARCHIVE_HASH}")
        print(f"Verified {verify_manifest(bundle)} manifest entries.")
        b = divided_bernoulli(101)
        check_first_field(json.loads(bundle.read("base/support-witnesses.json")), b)
        print("Verified all 25 first-field support ratios against the stored G values.")
        degree_100_crosscheck(bundle, check_upper_rows(bundle, b))


if __name__ == "__main__":
    main()
