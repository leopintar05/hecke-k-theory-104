# Hecke calculations for the conditional vanishing of integer K-groups through degree 104

Leo Pintar, École Polytechnique · Version 1.0 · 21 November 2026

In degree 104, the two rank-324 first-field Hecke tests and the classical calculations leave four possible primes. The rank-2,925 second-field evaluation matrix has nonzero determinant modulo each of them. Together with the earlier certificates this gives $K_{4k}(\mathbb Z)=0$ for $1\le k\le26$ under three arithmetic assumptions. [Two-page preprint](paper/hecke-through-104.pdf) · [LaTeX source](paper/main.tex)

The assumptions concern the strict Selmer group detected by Alpöge's bound, whether the computed Hecke spaces cover every admitted representation, and whether one residual representation at a fixed coefficient prime supplies all the traces compared together within each field. Page 2 states the assumptions and identifies the point where individual trace congruences must hold simultaneously.

## Reproducing the degree-104 result

With Python 3.12 and a C compiler installed, run the following commands from the root of the repository or from the unpacked computational supplement.

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -u reproduce.py
```

The final result line is `PASS degree 104: rank 2925, four prime fields, gcd 1.` The verification run started from the distributed ZIP and a fresh Python environment on Apple Silicon. Both rank-324 matrices matched entry for entry. The four complete prime-field determinants reproduced the reference JSON byte for byte. Direct evaluation of three rows agreed in all 8,775 entries. [Instructions](REPRODUCE.md) · [Run transcript](verification/degree104-replay.log) · [Recomputed residues](verification/degree104-result.json)

The earlier archive contains the first-field restrictions, second-field tests, mathematical reports and reconstruction instructions through degree 100. Its nested through-72 archive includes the earlier stages. For this release the stored first-field support ratios and the degree-76 through degree-100 block products were rechecked. The full matrix regeneration in this release run was for degree 104. [Through-100 certificate archive](certificates/K-theory-through-100-certificates.zip) · [Earlier record checks](verification/through100-records.log)

## Sources and contributions

Alpöge suggested the application to integer K-groups in Step 17, footnote 6. A separate public computation from the Grothendieck–Montpellier archive reports the range through degree 48. This package includes certificates for the fourteen additional multiples of four from 52 through 104. [Alpöge's announcement](https://alpo.ge/remark.pdf#page=4) · [Grothendieck–Montpellier report of 21 November](https://github.com/KokunoYumeto/grothendieck-montpellier-archive/blob/daddc785c90b0447f448ab372b1014efa2afb208/research/grh_20261120/intakes/K17/INTAKE_K17_20261121_part4.md#L69-L101)

Leo Pintar initiated and led this investigation into the $K$-theoretic consequences of Alpöge’s Step 17. OpenAI Codex contributed to the mathematical arguments, wrote and ran the verification code, and drafted the manuscript.
