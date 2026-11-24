# When mixing lowers similarity-sensitive entropy

*Exact concavity criteria and sharp bounds on entropy gain*

Four outcomes are safe, five can fail: a finite test for when pooling lowers similarity-sensitive diversity, and how
much an observation can raise it.

Version 0.1.0-candidate. DOI [10.5281/zenodo.24514480](https://doi.org/10.5281/zenodo.24514480) (all versions:
[10.5281/zenodo.24514479](https://doi.org/10.5281/zenodo.24514479)). Release page:
https://evidencepress.org/releases/similarity-entropy-concavity/

**Status: unrefereed research candidate.** Author: Anonymous. Text CC0-1.0; code MIT (see [LICENSE.md](LICENSE.md)).

## What this is about

**The problem.** Ecologists, forecasters and machine-learning researchers often measure the diversity or uncertainty
of a distribution in a way that knows some outcomes are alike: two closely related species count as less diverse than
two unrelated ones. A similarity matrix K records how alike each pair of outcomes is, from 0 (unrelated) to 1
(identical). The standard measure of this kind, the order-one similarity-sensitive entropy of Leinster and Cobbold,

    H_K(p) = - sum_i p_i log (K p)_i ,

averages minus the log of each outcome's *typicality*, its expected similarity to a random draw. With K the identity it
is ordinary Shannon entropy. Users expect it to behave like Shannon entropy: pooling two communities should not lower
diversity, and learning something should not, on average, raise uncertainty. Mathematically both statements are
*concavity* of H_K.

**What was known.** Concavity can fail. Gallego Posada, Vani, Schwarzer and Lacoste-Julien (GAIT, 2020) proved it for
two outcomes and conjectured it whenever K is positive definite and satisfies the multiplicative triangle inequality
K_ik >= K_ij K_jk. Bavaud showed failures without the triangle inequality and asked for exact conditions. Miller (2026)
refuted the conjecture with an 18-outcome example satisfying both conditions and asked for the least number of outcomes
at which concavity can fail. In ecology, a claim that the higher orders of the Ricotta–Szeidl family are concave was
later acknowledged by Ricotta, Szeidl and Pavoine to rest on a faulty proof.

**What is new.**
1. **A finite test** (Theorem 1). For any finite similarity matrix: merge identical rows, check that the reduced matrix
   L is positive definite, and check that one explicit matrix per outcome, T_j = L^-1 D(L e_j) + D(L e_j) L^-1 - e_j e_j^T,
   is positive semidefinite. This decides exactly whether pooling can ever lower the entropy. A simple sign condition
   on L^-1 covers every ultrametric similarity, such as one read from a dendrogram or phylogeny (Corollary 2).
2. **Four is safe, five can fail** (Theorem 4). For positive definite similarities satisfying the multiplicative triangle
   inequality, every matrix on four outcomes is safe (an exact computer-checked polynomial certificate), and a
   five-outcome example built on the complete bipartite graph K_{2,3} fails. This answers Miller's question.
3. **Trees** (Proposition 8, Corollary 9). On weighted trees the test becomes one matrix and runs in linear time. A star
   with equal edges needs at least 13 outcomes to fail, and minimal failing examples can be arbitrarily large.
4. **How big a failure can be** (Theorem 10). If K is positive semidefinite or satisfies the triangle inequality,
   an observation can raise expected entropy by at most half the expected posterior entropy, hence by at most
   (1/2) log n for n outcomes. Over positive semidefinite similarities the exact worst case is
   (1/2) log n - log log n + log 2 - 1 + o(1). A by-product: Bavaud's effective entropy, which is always concave, lies
   between one half and all of H_K (Proposition 12).
5. **Subadditivity** (Theorem 13). Joint entropy of paired variables under a product similarity is subadditive for
   every joint law only in degenerate cases.
6. **Other indices** (Theorem 14, Corollary 15). Among average-rarity indices, only Rao's quadratic entropy is concave
   for every positive definite metric similarity; a published ecological recommendation that orders two and above are
   safe must be narrowed to order two.
7. **Infinite spaces** (Theorem 16). On standard Borel spaces, concavity is decided by finite restrictions.

**The idea.** Measure a change in the distribution by the relative changes it causes in the typicalities instead of by
the change in probabilities. In those coordinates the curvature matrix of H_K depends linearly on the distribution,
so it is enough to check the corners of the probability simplex: one matrix per outcome. For four outcomes, every
four-point distance can be written as a sum of "cuts", and the corner conditions become seven polynomials with only
positive coefficients.

**Limits.** Failures are rare and usually tiny: the five-outcome example fails throughout a window of edge similarities
of width about 0.033 just below the positive definiteness threshold, by about 0.0003 nats at a = 12/17, and the exact
test finds no failure on a grid below that window; none of 2000 seeded random five-outcome examples failed (a numerical
observation, not a certificate). The worst case nevertheless grows with the number of outcomes. The four-outcome step
rests on an exact computer certificate, not a conceptual proof. Several equivalences (concavity, nonnegative expected
information gain, proper scoring) were already noted by Miller and by GAIT and are credited in the paper. The
asymptotic statements converge slowly (Remark 11 gives their effective ranges).

## Contents

| Path | What |
|---|---|
| `manuscript.pdf`, `manuscript.tex` | The paper (19 pages); Section 1 is a plain introduction with a plan of the paper; Section 11 gives practical guidance |
| `figures/mixing-example.pdf` | Figure 1: the K_{2,3} example and its mixing gap |
| `code/classify.py` | The exact rational concavity test of Theorem 1 for any similarity matrix given as JSON (`python classify.py kernel.json`) |
| `code/five-state.json` | The five-outcome example at a = 7/10, as input for `classify.py` |
| `code/verify.py` | The four-outcome polynomial certificate (Table 1) and the five-outcome interior curvature |
| `code/verify_independent.py` | A second, differently normalised four-outcome certificate (needs SymPy) |
| `code/verify_mixing.py` | The five-outcome Jensen gap 999949/999364 and a 2050-outcome tree example with gain above log 2 |
| `code/check_classifier.py` | Tests `classify.py`: star thresholds, the five-outcome example, the order-three check |
| `code/check_extensions.py` | Cut decomposition, the K_{2,3} failure interval, the star minimum, ultrametric sign checks, the seeded random experiment and the asymptotic table (needs SymPy and mpmath) |
| `code/figure.py` | Draws Figure 1 (needs NumPy and Matplotlib) |
| `results/` | Recorded outputs of the checks |
| `requirements.txt` | Pinned versions of the optional packages |
| `refs.bib` | The manuscript's eleven references |
| `AI_INDEX.md` | Machine-oriented map of claims, evidence and limits |
| `RESEARCH_GATES.json` | Prior-art searches, sources read, limitation triage |
| `REVIEW_SUMMARY.md` | Internal review: objections and how they were handled |
| `REPLAY_RECEIPT.md`, `ENVIRONMENT.txt` | Release replay and the cross-platform comparison policy |
| `MANIFEST.sha256` | SHA-256 of every file |

## Replay

From `code/`, with CPython 3.12 or later. `classify.py`, `verify.py`, `verify_mixing.py` and `check_classifier.py` use
the standard library only; the other scripts need the packages in `requirements.txt`.

```sh
cd code
python -B -O verify.py                 > ../results/verify.txt
python -B -O verify_independent.py     > ../results/verify_independent.txt      # SymPy
python -B -O verify_mixing.py          > ../results/verify_mixing.txt
python -B -O check_classifier.py       > ../results/check_classifier.txt
python -B -O classify.py five-state.json > ../results/classify-five-state.json
python -B -O check_extensions.py --output ../results/extensions.json           # SymPy + mpmath, about 75 s
python -B figure.py                                                            # NumPy + Matplotlib
```

`check_classifier.py` and `check_extensions.py` import `classify` from their own folder, so run them from `code/` and
without `-I`. Every script exits with status 1 on any failed check, also under `python -O` (no `assert` statements).

**Cross-platform policy** (declared before tagging; the same in CI, ENVIRONMENT.txt and REPLAY_RECEIPT.md): the four
text outputs and `results/classify-five-state.json` are exact and are compared byte for byte. In
`results/extensions.json` every exact value (integers, rationals, strings, booleans, structure) is compared exactly,
and floating-point fields (the mpmath asymptotic table at 50 digits and printed summaries such as the root of the
K_{2,3} determinant polynomial and the interval width) are compared to a relative tolerance of 1e-9; an interpreter
field, if present, is excluded. `figures/mixing-example.pdf` is regenerated by `figure.py`; only its embedded creation
date changes. On macOS arm64 under CPython 3.12.11 every output matched, and `extensions.json` was byte-identical.
GitHub Actions on Linux runs the exact checks and three negative controls as a blocking job, and `check_extensions.py`
with the figure as a separate, non-blocking job.

## Assurance

Written proofs with exact computational certificates. Internal AI-assisted review only: five Claude reviewer roles
(editor, methodology, domain, perspective, devil's advocate) run as two subagents; the lead's proof audit and clean
replay with negative controls; one revision round with extensions; then a re-review and a proof audit of the changed
material. No second-model (Astra/Codex) review took place. This is not peer review, independent reproduction or formal
verification. See [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).

## Cite

Anonymous (2026). *When mixing lowers similarity-sensitive entropy: exact concavity criteria and sharp bounds on
entropy gain* (Version 0.1.0-candidate). Evidence Press. https://doi.org/10.5281/zenodo.24514480
