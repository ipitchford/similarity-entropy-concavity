# AI index: When mixing lowers similarity-sensitive entropy (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction with a plan of the paper; Section 11 gives practical guidance). DOI 10.5281/zenodo.24514480
(concept 10.5281/zenodo.24514479). Status: unrefereed candidate; internal AI-assisted review only.

## Normalisation

n states (outcomes). Similarity matrix K: symmetric, entries in [0,1], unit diagonal. Law p in the simplex Delta_n.
Typicality q = Kp. Order-one similarity-sensitive entropy H_K(p) = -sum_{i: p_i > 0} p_i log q_i (natural logs, nats);
K = I gives Shannon entropy. Concave = H_K(lambda p + (1-lambda) r) >= lambda H_K(p) + (1-lambda) H_K(r).
PD/PSD = positive (semi)definite. MTI = multiplicative triangle inequality K_ik >= K_ij K_jk; with positive entries
K = exp(-d) for a metric d ("metric similarity"). Reduction: merge identical rows; L = reduced matrix of order m.
Vertex matrices T_j = L^{-1} D(L e_j) + D(L e_j) L^{-1} - e_j e_j^T. Observation: random posterior P with E P = p;
G = E H_K(P) - H_K(p) (excess of expected posterior entropy over prior entropy). K_{2,3}(a): five states in parts of
two and three, a across parts, a^2 within a part. Theorem numbering follows the manuscript (one shared counter).
GAIT = Gallego Posada, Vani, Schwarzer, Lacoste-Julien, PMLR 108 (2020). Miller = arXiv:2601.03064v3.
Bavaud = doi:10.1007/s00357-021-09395-4. RS = Ricotta-Szeidl, doi:10.1016/j.tpb.2006.06.003; RSP =
Ricotta-Szeidl-Pavoine, doi:10.1016/j.ecolind.2021.107971. LC = Leinster-Cobbold, doi:10.1890/10-2402.1.

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A (main) | H_K concave on Delta_n iff L is PD and every vertex matrix T_j is PSD; then H_L strictly concave; equality iff class masses agree | Theorem 1; Section 3 | written proof; Hessian from GAIT Supplement C.2 (credited); affine reparametrisation new |
| B | Sign condition: L^{-1} with nonpositive off-diagonal and nonnegative row sums implies concavity; covers every ultrametric similarity and Laplace kernels on the line | Corollary 2 | written proof; Miller's identity (App. B.1) and Bavaud Theorem 6c credited |
| C | Concavity is equivalent to nonnegative expected information gain, nonnegative GAIT divergence and proper tangent scoring; Theorem 1 decides each | Section 4 | credited recap (Miller Sec. 4.1, Lemma A.1, Sec. 6.2; GAIT Def. 2); only the decidability is new |
| D | Every similarity with 0 < K_ij < 1 off the diagonal is induced by a deterministic map from Laplace-kernel points | Proposition 3 | written proof; bears on Miller's representation-change question (Section 6.2) |
| E (main) | (a) Every MTI similarity of order at most four has concave entropy; (b) a PD MTI similarity of order five (K_{2,3}(a), a = 12/17 or 7/10) is not concave; hence five is the minimum order for which PD + MTI fail to guarantee concavity | Theorem 4; Section 5; Lemmas 5-6; Table 1 | (a) exact polynomial certificate (two independent normalisations) plus written reduction; (b) exact rational values; answers Miller's minimum-cardinality question; PD of four-point metric kernels credited to Meckes |
| F | Weighted trees: one matrix decides concavity, linear time; equal-edge star with m leaves concave iff m <= C(a) = 2(1+a)/(a(1-a)); min C = 6 + 4 sqrt 2 at a = sqrt 2 - 1, so the least failing equal-edge star has 13 states; minimal failing examples of unbounded order | Proposition 8; Corollary 9 | written proof; checked symbolically and by the classifier |
| G (main) | If K is PSD or MTI: H_K(p) >= (1/2) E H_K(P), so G <= (1/2) log n. (a) Over PSD similarities sup G = (1/2) log n - log log n + log 2 - 1 + o(1), attained by binary observations; (b) over strictly positive PD tree similarities with fair binary observations sup G >= (1/2) log n - log log n - (5/2) log 2 - 1/2 + o(1), and the factor 1/2 is asymptotically optimal; (c) over MTI similarities without PSD, sup G = (1/2) log n + O(1) for fair binary observations | Theorem 10; Remark 11 | written proof; Remark 11 effective ranges from check_extensions.py |
| H | Bavaud's effective entropy E satisfies (1/2) H_K <= E <= H_K under PSD or MTI | Proposition 12 | written proof; upper bound is Bavaud's |
| I (main) | Joint entropy under K tensor L subadditive for every joint law iff one factor is all-ones, or a partition with at least two classes and the other entropy concave | Theorem 13 | written proof; answers the second part of Bavaud's Section 6 question |
| J (main) | Among indices sum_i p_i phi((Kp)_i), only affine nonincreasing phi are concave for every PD MTI similarity (order seven suffices); in the RS family only order two (Rao) is safe; explicit six-state counterexamples to RS Theorem 2 for every order above two | Theorem 14; Corollary 15 | written proof; RS 2006 error acknowledged in RSP 2021 (credited); order-three exact check 1557/25000 and gap 1557/50000000000 |
| K | On standard Borel spaces concavity holds iff it holds on every finite restriction; explicit plug-in bias bound | Theorem 16 | written proof; sampling device credited to Miller's proof of Lemma A.1 |

## Evidence

| Check | Script | Output | Scope |
|---|---|---|---|
| Four-state certificate: seven Schur-complement polynomials in cut weights with all coefficients positive; supports 255/63/63/1557/1557/1245/15172, minima 2..8, maxima up to 464400; five-state interior curvature +2748940385160/93148730546881 at a = 7/10 | [code/verify.py](code/verify.py) | [results/verify.txt](results/verify.txt) | exact integers and rationals; standard library |
| Second four-state certificate with the diagonally congruent normalisation (e.g. 30,562 terms for {1,2,3}), all positive | [code/verify_independent.py](code/verify_independent.py) | [results/verify_independent.txt](results/verify_independent.txt) | exact (SymPy polynomial rings over ZZ) |
| K_{2,3}(12/17): 1 - 2a^2 = 1/289, Jensen ratio 999949/999364 > 1 (J about 0.00029260 nats); 2050-state double star with exact exp(4G) > 16 | [code/verify_mixing.py](code/verify_mixing.py) | [results/verify_mixing.txt](results/verify_mixing.txt) | exact rationals |
| Classifier (Theorem 1 in rational arithmetic) against 15 star-threshold comparisons, the five-state example and the RS order-three gap 1557/50000000000 | [code/check_classifier.py](code/check_classifier.py), [code/classify.py](code/classify.py) | [results/check_classifier.txt](results/check_classifier.txt), [results/classify-five-state.json](results/classify-five-state.json) | exact rationals |
| Cut weights recover all six distances; K_{2,3}(a) part-block determinant (1-a^2)(1-2a^2), Schur block, root 0.6736480786 of the determinant polynomial, nonconcave and PD throughout (0.673648, 1/sqrt 2), width about 0.033; classifier agrees on the grid 3/5..7/10 (concave up to 67/100); star minimum 6 + 4 sqrt 2; 60 random rational ultrametrics satisfy the sign condition; 0 of 2000 seeded random five-state PD MTI kernels fail; asymptotic table (Remark 11) | [code/check_extensions.py](code/check_extensions.py) | [results/extensions.json](results/extensions.json) | symbolic/rational, except the asymptotic table (mpmath, 50 digits) and float summaries, compared to 1e-9 relative; the random experiment is a numerical observation, not a certificate |
| Figure 1 | [code/figure.py](code/figure.py) | [figures/mixing-example.pdf](figures/mixing-example.pdf) | illustration |

Numbers quoted on the release page and in the audio: four and five outcomes; 18 (Miller's example); window of width
about 0.033; about 0.0003 nats; 0 of 2000; half the log of the number of outcomes; 13 states for the smallest failing
equal-edge star; factor two for the effective entropy. All in results/ or proved.

## Limits

The four-state step rests on an exact computer certificate, not a conceptual proof. Failures found are small and
confined to narrow parameter ranges, but the worst case grows like (1/2) log n. The random five-state experiment is
seeded and uncertified. The asymptotic statements converge slowly (Remark 11). The equivalences of Section 4 are
credited to Miller and GAIT. Open problems in Section 12: a structural four-point proof, a pentagonal characterisation
at five points, the constant term for PD MTI fair observations, the smallest failing tree, order six in Theorem 14, and
whether the factor two in Proposition 12 is attained.

## Replay and records

Replay commands and cross-platform policy (exact outputs byte for byte; float fields of extensions.json to 1e-9
relative): [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art
and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Dependencies:
[requirements.txt](requirements.txt). Licences: [LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff).
Hashes: [MANIFEST.sha256](MANIFEST.sha256).
