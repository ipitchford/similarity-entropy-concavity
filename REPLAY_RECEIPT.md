# Replay receipt (0.1.0-candidate, 24 November 2026)

**Environment.** macOS arm64; CPython 3.12.11 (uv-managed, `uv run --no-project --python 3.12 [--with ...] python`;
the version was printed inside the run with `sys.version`); SymPy 1.14.0 and mpmath 1.3.0 where needed (versions
printed inside the run); NumPy 2.3.5 and Matplotlib 3.10.8 for the figure. Each run used a fresh copy of the approved
candidate, cwd = its `code/` folder, `PYTHONDONTWRITEBYTECODE=1`.

| Run | Command | Exit | Output vs `results/` |
|---|---|---|---|
| four-state certificate, plain and -O | `python [-O] -I -B verify.py` | 0, 0 | byte-identical to `verify.txt` |
| second four-state certificate, -O | `python -O -I -B verify_independent.py` (SymPy) | 0 | byte-identical to `verify_independent.txt` |
| five-state gap and tree example, plain and -O | `python [-O] -I -B verify_mixing.py` | 0, 0 | byte-identical to `verify_mixing.txt` |
| classifier checks, plain and -O | `python [-O] -B check_classifier.py` | 0, 0 | byte-identical to `check_classifier.txt` |
| classifier on the five-state example, -O | `python -O -I -B classify.py five-state.json` | 0 | byte-identical to `classify-five-state.json` |
| extensions, -O (about 71 s on this host) | `python -O -B check_extensions.py --output new.json` (SymPy, mpmath) | 0 | ALL CHECKS PASSED; JSON byte-identical to `extensions.json` |
| figure | `python -B figure.py` (NumPy, Matplotlib) | 0 | `figures/mixing-example.pdf` regenerated; only /CreationDate differs |

**Negative controls (under -O).**
1. `check_classifier.py` with the expected order-three gap 1557/50000000000 changed to 1558/50000000000:
   "CHECK FAILED: gap == F(1558, 50000000000)", exit 1.
2. `verify.py` with the expected Table 1 count for {1,2,3} changed from 15172 to 15173:
   "CHECK FAILED: ((1, 2, 3), (15172, 8, 464400))", exit 1.
3. `verify_mixing.py` with the expected Jensen ratio 999949/999364 changed to 999949/999365:
   "CHECK FAILED: 1-2*a*a == F(1, 289) and ratio == F(999949, 999365)", exit 1.

**After the replay** only the title-page DOI line of `manuscript.tex`, the rebuilt PDF (three pdflatex passes, 19
pages; a word-level text comparison with the approved PDF shows only the DOI line) and the release records changed; the
code, figure and recorded outputs are byte-identical to the replayed candidate.

**Cross-platform policy** (declared before tagging; see ENVIRONMENT.txt): exact outputs are compared byte for byte.
In `results/extensions.json` exact values are compared exactly and floating-point fields to a relative tolerance of
1e-9. The GitHub Actions workflow `.github/workflows/replay.yml` reruns the exact checks under `python -O` on Linux,
compares byte for byte and runs the three negative controls (blocking); `check_extensions.py` and the figure run in a
separate non-blocking job under the same policy.

**Scope.** The computations certify the four-state polynomial positivity, the five-state failure and the stated exact
values. The seeded random experiment (0 of 2000) is a numerical observation, not a certificate. The proofs carry the
general theorems.
