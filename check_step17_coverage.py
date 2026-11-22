"""Verify the Bernoulli field choices and the deductions from the stored Hecke results."""

from argparse import ArgumentParser
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
import json
from math import comb, gcd, prod
from pathlib import Path
import re
from zipfile import ZipFile


CUTOFF = 2**31
RELEASE_SHA256 = "1b531cc567b90c2ddbacbfdaca7dd2ca1a6ae3df7cc679bb1e64e1f88b6d9be1"
RELEASE_FILE = "hecke-certificates-v1.1.zip"


def check(condition, description):
    if not condition:
        raise ValueError(description)


def bernoulli_quotients(d, last):
    """Entry n is B_(n+1,chi)/(n+1), for chi the Legendre character mod d."""
    residues = {a * a % d for a in range(1, d)}
    terms = [(a, 1 if a in residues else -1) for a in range(1, d)]
    values = []
    for n in range(last + 1):
        # Coefficients of (exp(d*t)-1) times the exponential generating series.
        moment = sum(sign * a ** (n + 1) for a, sign in terms)
        known = sum(comb(n + 1, j) * d**j * values[n + 1 - j]
                    for j in range(2, n + 2))
        values.append(Fraction(moment - known, (n + 1) * d))
    return values


def num_product(values, e):
    return abs(values[e].numerator * values[2 * e].numerator)


def odd_part(value):
    while value % 2 == 0:
        value //= 2
    return value


def field_units():
    b3 = bernoulli_quotients(3, 104)
    b7 = bernoulli_quotients(7, 104)
    check((b3[0], b7[0]) == (Fraction(-1, 3), -1), "first Bernoulli units")
    odd_common = {}
    for e in range(2, 53, 2):
        for values in (b3, b7):
            for i in (0, e, 2 * e):
                check(values[i].denominator < CUTOFF, f"large denominator at e={e}")
        common = gcd(num_product(b3, e), num_product(b7, e))
        if odd_part(common) > 1:
            odd_common[e] = odd_part(common)
    check(odd_common == {12: 47, 18: 37, 36: 37**2, 42: 43**2}, "odd common factors")
    check(gcd(num_product(b3, 52), num_product(b7, 52)) == 4, "degree-104 gcd")
    check(tuple(v[i].denominator for v in (b3, b7) for i in (52, 104))
          == (3, 9, 1, 49), "degree-104 denominators")
    print("Field choice verified for all 26 exponents. At e=52, gcd(A3,A7) = 4.")
    return b3


def read_json(archive, name):
    return json.loads(archive.read(name))


def doubled_target(e):
    return 2 * (2 ** (2 * e) + 2**e + 1)


def log_text(archive, name):
    return archive.read(name).decode()


def check_through72(archive, b3):
    excluded = defaultdict(set)

    def compact(e, modulus, residue):
        check(gcd(doubled_target(e), modulus) == gcd(residue, modulus) == 1,
              f"compact exclusion at e={e} modulo {modulus}")
        excluded[e].add(modulus)

    with ZipFile(BytesIO(archive.read("K-theory-through-40-proof-bundle.zip"))) as first:
        for line in log_text(first, "verification/two-prime-confirmed.txt").splitlines():
            if match := re.match(r"e=(\d+),", line):
                e = int(match[1])
            if match := re.search(r"(\d+): B_(\d+) = (\d+) mod (\d+)", line):
                p, n, residue, modulus = map(int, match.groups())
                if p >= CUTOFF:
                    check(p == modulus and n == p - 1 - e and gcd(residue, p) == 1,
                          f"cyclotomic exclusion at e={e} and p={p}")
                    excluded[e].add(p)
        for e, p in ((16, 25078104767246023), (20, 552972895573825165951583191)):
            record = log_text(first, f"verification/minus7-{e}-confirmed.txt")
            modulus = int(re.search(r"p=(\d+)", record)[1])
            residue = int(re.search(r"norm of Schur determinant modulo p = (\d+)", record)[1])
            check(modulus % p == 0, f"recorded modulus at e={e}")
            compact(e, p, residue)

    with ZipFile(BytesIO(archive.read("K-theory-through-48-certificates.zip"))) as next_cases:
        for name in ("logs/compact-minus-seven.txt", "logs/degree44-compact-native-check.txt"):
            record = log_text(next_cases, name)
            p, e = map(int, re.search(r"p=(\d+) e=(\d+)", record).groups())
            residue = int(re.search(r"norm of Schur determinant modulo p = (\d+)", record)[1])
            compact(e, p, residue)
        lines = log_text(next_cases, "logs/inert-desc48.log").splitlines()
        line = next(line for line in lines if "detSchur, norm = " in line)
        e, p, *_, residue = map(int, line.split(" = ")[1].split())
        compact(e, p, residue)

    # A rational norm or a pair of split-coefficient determinants treats every embedding.
    for line in log_text(archive, "logs/degree-52/packed.log").splitlines():
        if line.startswith("prime,"):
            p, _, _, _, norm = map(int, line.split(": ")[1].split())
            compact(26, p, norm)
    split_residues = defaultdict(dict)
    for name in archive.namelist():
        if not re.match(r"logs/degree-(52|56|60|64|68)/", name):
            continue
        for line in log_text(archive, name).splitlines():
            label, _, data = line.partition(" = ")
            if not re.fullmatch(r"[\d ]+", data) or not ("schur" in label.lower()):
                continue
            if "sign" in label.lower():
                continue
            e, p, *_, residue = map(int, data.split())
            if e not in (26, 28, 30, 32, 34):
                continue
            if "norm" in label.lower():
                compact(e, p, residue)
            elif match := re.search(r"-(0|1)\.log$", name):
                split_residues[e, p][int(match[1])] = residue
    for (e, p), values in split_residues.items():
        if len(values) == 2:
            compact(e, p, values[0] * values[1])
    record = log_text(archive, "logs/degree-72/composite.log")
    for matrix in ("evaluation", "numerator"):
        check(f"{matrix} right inverse verified 1944" in record, f"degree-72 {matrix} inverse")
    match = re.search(r"verified both units: e, modulus, rank, seed = 36 (\d+) 1944 ", record)
    check(match is not None, "degree-72 composite modulus")
    compact(36, int(match[1]), 1)

    # Dividing by these small integers is safe whether or not they are prime.
    small = (*range(2, 501), 3343, 7691, 8609, 13687, 15083, 16519, 18541,
             28447, 36293, 38177, 7547263, 42800473)
    check(max(small) < CUTOFF, "cutoff for lower-degree factors")
    for e in range(2, 37, 2):
        remaining = num_product(b3, e)
        moduli = prod(excluded[e])
        while (part := gcd(remaining, moduli)) > 1:
            remaining //= part
        for q in small:
            while remaining % q == 0:
                remaining //= q
        check(0 < remaining < CUTOFF, f"uncovered lower-degree divisor at e={e}")
    print("Through degree 72, every large divisor is covered by a stored exclusion.")


def check_degree104(archive, b3):
    result = read_json(archive, "degree104-result.json")
    with ZipFile(BytesIO(archive.read("prior/Degree-104-rank-324-certificates.zip"))) as first:
        support = read_json(first, "support.json")
        check(support["rank"] == 324 and result["rank"] == 2925, "matrix dimensions")
        # Decimal parses long integer literals exactly, without a global digit-limit change.
        determinants = [int(Decimal(first.read(name).decode().strip()))
                        for name in ("inert-2-det.txt", "split-7-sum-det.txt")]
    primes = list(map(int, support["large_candidate_primes"]))
    p1, p2, p3, p4 = primes
    check(all(p >= CUTOFF for p in primes), "candidate cutoff")
    check(abs(b3[52].numerator) == 2 * p1 * p2 * p3, "first numerator factorization")
    check(abs(b3[104].numerator) == 2 * 23 * 47 * p4, "second numerator factorization")
    expected_gcd = 3**1578 * 5**33 * 11**31 * 19**18 * 23**6 * 29 * 79 * p1 * p2 * p3
    check(gcd(*determinants) == expected_gcd, "rank-324 determinant gcd")

    modulus = prod(primes)
    determinant = int(result["determinant_mod_M"])
    check(int(result["modulus"]) == modulus, "degree-104 modulus")
    check([int(row["prime"]) for row in result["components"]] == primes, "residue moduli")
    for row in result["components"]:
        check(determinant % int(row["prime"]) == int(row["determinant"]), "residue")
    check(gcd(doubled_target(52), modulus) == 1, "degree-104 elimination scalar")
    check(gcd(determinant, modulus) == 1, "rank-2925 determinant residue")
    print("Degree 104: both stored rank-324 determinants give the stated gcd.")
    print("All four rank-2925 residues agree; the determinant and doubled target are units.")


def check_lower(archive, b3):
    with ZipFile(BytesIO(archive.read("K-theory-through-72-certificates.zip"))) as lower:
        check_through72(lower, b3)
    rows = read_json(archive, "base/support-witnesses.json")
    check([r["e"] for r in rows] == list(range(2, 51, 2)), "first-field degree range")
    for row in rows:
        e = row["e"]
        a = abs(b3[e].numerator)
        b = abs(b3[2 * e].numerator)
        ratio = Fraction(int(row["G"]), a * b**3)
        for factors, value in ((row["ratio_num"], ratio.numerator),
                               (row["ratio_den"], ratio.denominator)):
            check(all(2 <= q <= 337 and v > 0 for q, v in factors), "small factors")
            check(prod(q**v for q, v in factors) == value, f"support ratio at e={e}")
    cases = read_json(archive, "certificates.json")
    weights = read_json(archive, "weights.json")
    check(sorted(map(int, cases)) == list(range(38, 51, 2)), "signed degree range")
    for key, case in cases.items():
        e = int(key)
        modulus = int(case["modulus"])
        check(gcd(doubled_target(e), modulus) == 1, f"elimination scalar at e={e}")
        check(gcd(int(case["determinant_mod_M"]), modulus) == 1, f"determinant at e={e}")
        removed = [(int(q), v) for q, v in weights[key]["removed"].items()]
        check(all(2 <= q < CUTOFF and v > 0 for q, v in removed), f"small factors at e={e}")
        check(num_product(b3, e) == modulus * prod(q**v for q, v in removed),
              f"uncovered large prime at e={e}")
    print("Verified 25 first-field support identities and seven signed upper-degree exclusions.")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("release", nargs="?", type=Path,
                        help="path to the v1.1 ZIP or the repository certificates directory")
    location = parser.parse_args().release
    if location is None:
        here = Path(__file__).resolve().parent
        location = here / RELEASE_FILE
        if not location.is_file():
            location = here / "certificates"
    inputs = (
        ("Degree-104-vanishing-certificates.zip",
         "8f3829d15aaf8b605c5984c6bd90e160e614b8437ecadb47a7e6c72bdcdd5bc9", check_degree104),
        ("K-theory-through-100-certificates.zip",
         "c1a718a1d6d5b0e268c91a28fcfebfc2cc462cb7cd058cde6db470a24243ccc2", check_lower),
    )
    if location.is_dir():
        data = [(location / name).read_bytes() for name, _, _ in inputs]
    elif location.is_file():
        packed = location.read_bytes()
        check(sha256(packed).hexdigest() == RELEASE_SHA256, "unexpected v1.1 release digest")
        with ZipFile(BytesIO(packed)) as release:
            data = [release.read("certificates/" + name) for name, _, _ in inputs]
    else:
        parser.error(f"place {RELEASE_FILE} beside this script or supply the path to it "
                     "or to the repository certificates directory")
    b3 = field_units()
    for (name, digest, verify), payload in zip(inputs, data):
        check(sha256(payload).hexdigest() == digest, f"unexpected digest for {name}")
        with ZipFile(BytesIO(payload)) as archive:
            verify(archive, b3)


if __name__ == "__main__":
    main()
