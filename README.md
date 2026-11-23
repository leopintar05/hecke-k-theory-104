# Vanishing of integer K-groups through degree 104 from Alpöge’s Step 17

Leo Pintar, École Polytechnique · Paper version 5 · 23 November 2026

We obtain $K_{4k}(\mathbb Z)=0$ for $1\le k\le26$, assuming the arithmetic construction in Alpöge’s pinned Step 17. At level 2 with trivial character, the source’s zero-length criterion is exactly the vanishing of the rational strict Selmer group used here. The paper proves this comparison for both auxiliary fields. The construction supplies one representation at one coefficient prime for all the traces tested in a field. Three Bernoulli unit tests then restrict the possible source to the compact space covered by the matrices. [Full paper (seven pages)](paper/hecke-through-104.pdf) · [Detailed proof](paper/coverage.md) · [LaTeX source](paper/main.tex)

In degree 104 the Bernoulli criterion assigns four large primes to the second field. Two rank-324 determinants exclude all other large primes in the first field. Three of the four divide their gcd. The fourth requires a second-field test because its first-field Bernoulli condition fails. A rank-2,925 determinant excludes all four. Ordinary Bernoulli residues give independent cyclotomic exclusions at the first two primes. [Degree-104 residues](verification/degree104-result.json) · [Cyclotomic check](check_reflection104.py)

The standard-library verifier recomputes the field choices for all 26 degrees and checks the recorded determinants and prime exclusions.

```sh
python3 check_step17_coverage.py
```

For a fresh reconstruction of the matrices, see the [reproduction instructions](REPRODUCE.md). The full degree-104 replay recorded for v1.1 reproduced both rank-324 matrices and all four prime-field residues. Three rows of the rank-2,925 matrix were also evaluated directly, checking 8,775 entries. [Replay transcript](verification/degree104-replay.log) · [Degree-104 residues](verification/degree104-result.json)

## Sources and contributions

Levent Alpöge suggested the application to $K_{4k}(\mathbb Z)$ in Step 17, footnote 6. The Grothendieck–Montpellier Archive’s *KTOWER*, checkpoint 11, reports finite certificates through degree 80 and supplies the joint-representation and Bernoulli-unit results used here. [Alpöge’s remark](https://alpo.ge/remark.pdf#page=4) · [Montpellier’s KTOWER report](https://github.com/KokunoYumeto/grothendieck-montpellier-archive/blob/edb183a9af46cb9568e0d3caf2248e98b491e6fb/raw/research_returns/KTOWER_checkpoint_11_20261122T0001Z_unzipped/KTOWER.tex)

Leo Pintar initiated and led this investigation into the $K$-theoretic consequences of Alpöge’s Step 17. OpenAI Codex contributed to the mathematical arguments, wrote and ran the verification code, and drafted the manuscript.
