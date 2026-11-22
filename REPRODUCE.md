# Reproducing the Hecke calculations

For version 3, run `python3 check_step17_coverage.py` from the checkout. For the standalone release, download the verifier and `hecke-certificates-v1.1.zip` into one directory and run the same command. Python 3.9 or later suffices for this check. It verifies the certificate digests, recomputes the field choice for the three Bernoulli units at each exponent, and checks the stored arithmetic through degree 104.

For the earlier, independent degree-104 self-dual case in the first field, run `python3 check_adjoint104.py`. The verifier uses only the Python standard library. It builds the full 16-dimensional weight-53 cusp space from Fourier expansions, computes the trace tests at 2 and 7, and checks the complete factorization of the gcd of their determinants.

Run these commands from the directory containing `reproduce.py`. The degree-104 run uses Python 3.12, a C99 compiler with unsigned 128-bit integer support (Clang or GCC), and the two pinned Python packages. Allow about 2 GB of free scratch disk space.

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -u reproduce.py
```

The last result line starts with **PASS degree 104** and reports rank **2925**, four prime fields and gcd **1**. Each run starts in a fresh temporary directory. To retain the new matrices, pass `--work-dir ./replay-output` with a new or empty directory. Use `--workers 2 --threads 1` on a smaller machine.

The full run includes `check_adjoint104.py`. It rebuilds the other classical tests over both fields, reconstructs both rank-324 matrices from neighbor sums, compares every entry to the distributed copy, and recomputes their integer determinants. It then generates the rank-2,925 matrix from the signed basis, computes its determinant in each of the four prime fields and checks the Chinese remainder calculation. A second implementation evaluates all 2,925 columns in each of three rows directly. The prime-field results and the packed-matrix digest are recorded in `degree104-result.json` inside the degree-104 certificate ZIP.

To check the **stored** earlier arithmetic records, run this additional command. It checks 25 first-field support ratios, the seven block-determinant products for degrees 76 through 100, and a prime-field cross-check in degree 100.

```sh
python check_through100.py certificates/K-theory-through-100-certificates.zip
```

The through-100 archive includes the detailed mathematical reports, matrix generators and logs. Its nested `K-theory-through-72-certificates.zip` contains the earlier stages and their own reconstruction instructions. The degree-104 archive includes the rank-324 and classical certificates under `prior/`, the second-field generator under `second-field/`, and the direct row evaluator under `audit/`.

The two original certificate ZIPs are preserved byte for byte. Their SHA-256 values are

```text
8f3829d15aaf8b605c5984c6bd90e160e614b8437ecadb47a7e6c72bdcdd5bc9  certificates/Degree-104-vanishing-certificates.zip
c1a718a1d6d5b0e268c91a28fcfebfc2cc462cb7cd058cde6db470a24243ccc2  certificates/K-theory-through-100-certificates.zip
```
