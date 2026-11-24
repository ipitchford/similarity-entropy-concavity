# Review summary (public, 0.1.0-candidate)

This records the internal review of the candidate: who looked, what they raised, and what was done. It is
AI-assisted internal review, not peer review or independent reproduction.

## What ran

1. **Panel (round 1).** Five Claude reviewer roles, editor, methodology (proof audit of every section and the
   scripts), domain (literature, attribution and novelty), perspective (readability and significance) and devil's
   advocate, run as two subagents. No reviewer found a mathematical error that breaks a theorem.
2. **Clean replay and proof audit by the lead.** Every bundle script rerun from a fresh copy; outputs equal to the
   recorded ones. The replay found that three scripts used `assert` statements, which pass silently under
   `python -O`; the methodology reviewer confirmed this with deliberately falsified copies.
3. **One revision round** with extensions: the sign condition covering every ultrametric similarity (Corollary 2);
   the factor-two comparison with Bavaud's effective entropy (Proposition 12); an explicit constant for the fair
   double star (Theorem 10(b)); practical guidance with a table of safe similarity families and quantified rarity of
   failures (Section 11); and the four/five parallel with Rao's quadratic entropy (Remark 7).
4. **Verification.** A re-review of everything that changed, checked against the manuscript rather than the fix list,
   and a proof audit of the new and rewritten proofs (Corollary 2, Lemmas 5-6, Theorem 4, Theorem 10(b),
   Proposition 12, Corollary 15 and the changed attributions), with a fresh rerun of every script and a negative
   control. No correctness issue remained; the minor wording items found were applied.

No second-model (Astra/Codex) adversarial review took place for this release.

## Sources

All eleven cited works were obtained and read: Leinster–Cobbold, GAIT (Gallego Posada et al.), Bavaud, Miller
(arXiv v3), Ricotta–Szeidl 2006, Ricotta–Szeidl–Pavoine 2021 (its preprint text, which contains the acknowledgement of
the 2006 error), Meckes, Leinster (magnitude), Botta-Dukát, Reeve et al. and Deza–Laurent (cited only for the classical
cut decomposition, whose coefficients the paper states and the code checks symbolically). Attributions were checked
against those texts: the Hessian is GAIT's (Supplement C.2), the equivalences of Section 4 are Miller's and GAIT's, the
sign identity of Corollary 2 is Miller's (Appendix B.1), the Stieltjes property of ultrametric inverses is Bavaud's
(Theorem 6c), positive definiteness of four-point metric kernels is Meckes's, and the K_{2,3} singularity is noted by
Leinster and Meckes.

## Main objections and dispositions

| Objection | Disposition |
|---|---|
| Correctness | No error that breaks a theorem was found in round 1, in the re-review or in the proof audit. |
| A sign sentence in the five-state interior check was reversed | Corrected: the positive number is the second derivative of H_K; the Hessian formula gives its negative (Section 5.1). |
| The four-state proof used positive definiteness before proving it | New Lemma 5 (positive definiteness of metric similarities on at most four points, credited to Meckes, with a self-contained cut proof) and Lemma 6 (blocks); the case analysis is now explicit. |
| The minimum theorem conflated "PD and MTI" with "every MTI matrix of order at most four" | Theorem 4 split into (a), for MTI matrices of order at most four, and (b), the five-state PD MTI failure. |
| "Factor one half is optimal" was only asymptotic | Theorem 10(b) says the ratio tends to one half; Remark 11 gives effective ranges computed by `check_extensions.py`. |
| Known equivalences presented as results | Section 4 rewritten as a credited recap of Miller and GAIT; only the decidability by Theorem 1 is claimed. |
| Credit for Ricotta–Szeidl and its later correction | Section 9 states what Ricotta and Szeidl claimed and that Ricotta, Szeidl and Pavoine acknowledged the error; Corollary 15 gives explicit counterexamples to their Theorem 2 for every order above two. |
| Missing related work | Meckes (four-point negative type and the K_{3,2} example), Leinster (magnitude), Reeve et al. (partitioning diversity) and Botta-Dukát (whose order-two-and-above recommendation is narrowed to order two) are cited. |
| Readability | New context-first introduction (problem, why concavity matters, what was known, numbered contributions with credit, roadmap); a notation section; a purpose sentence opens every section; figure and table captions explain themselves. |
| Checks could pass vacuously under optimisation | `assert` statements replaced by explicit checks that exit with status 1, also under `python -O`; new `check_extensions.py`; negative controls confirm detection. |
| Floating-point evidence for Gaussian kernels from reviewers | Not adopted: the reviewers' float-only observations were not certified, so the paper does not claim them. Random-kernel counts were replaced by a seeded experiment run by the package. |

## Limits that remain

The open problems in Section 12 of the paper: a structural proof of four-point safety, a pentagonal characterisation
at five points, the constant term for positive definite MTI similarities with fair observations, the smallest failing
weighted tree, lowering order seven to six in Theorem 14, and whether the factor two in Proposition 12 is attained. The
literature search is recorded in RESEARCH_GATES.json.
