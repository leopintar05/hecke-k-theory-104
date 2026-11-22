# Hecke calculations for the conditional vanishing of integer K-groups through degree 104

Leo Pintar, École Polytechnique · Paper version 3 · 22 November 2026

We obtain $K_{4k}(\mathbb Z)=0$ for $1\le k\le26$, conditional on the strict Selmer comparison in Alpöge’s full Step 17 construction. The construction supplies one representation at one coefficient prime for all the traces tested in a field. Three Bernoulli unit tests then restrict the possible source to the compact space covered by the matrices. [Two-page paper](paper/hecke-through-104.pdf) · [Detailed proof](paper/coverage.md) · [LaTeX source](paper/main.tex)

In degree 104 the first field excludes every large prime outside a four-prime list using two rank-324 determinants. In the second field the rank-2,925 determinant is coprime to the product of those four primes. The verifier recomputes the Bernoulli field choices for all 26 degrees and checks the stored determinant deductions. It uses the Python standard library.

```sh
python3 check_step17_coverage.py
```

For a fresh reconstruction of the matrices, see the [reproduction instructions](REPRODUCE.md). The full degree-104 replay recorded for v1.1 reproduced both rank-324 matrices and all four prime-field residues. Three rows of the rank-2,925 matrix were also evaluated directly, checking 8,775 entries. [Replay transcript](verification/degree104-replay.log) · [Degree-104 residues](verification/degree104-result.json)

## Sources and contributions

Levent Alpöge suggested the application to $K_{4k}(\mathbb Z)$ in Step 17, footnote 6. The Grothendieck–Montpellier Archive’s *KTOWER*, checkpoint 11, reports finite certificates through degree 80 and supplies the joint-representation and Bernoulli-unit results used here. [Alpöge’s remark](https://alpo.ge/remark.pdf#page=4) · [Montpellier’s KTOWER report](https://github.com/KokunoYumeto/grothendieck-montpellier-archive/blob/edb183a9af46cb9568e0d3caf2248e98b491e6fb/raw/research_returns/KTOWER_checkpoint_11_20261122T0001Z_unzipped/KTOWER.tex)

Leo Pintar initiated and led this investigation into the $K$-theoretic consequences of Alpöge’s Step 17. OpenAI Codex contributed to the mathematical arguments, wrote and ran the verification code, and drafted the manuscript.
