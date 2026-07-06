# Planned contrasts: rationale and implementation

## What the reviewer suggested

The reviewer did not propose changing the statistical model. Instead, they proposed
**changing the contrast coding** of the factors so that each model coefficient maps
directly onto one of the study's hypotheses.

Instead of fitting with default treatment coding and then using `emmeans` to answer
H1 and H2, the coefficients themselves become the tests of interest.

---

## Does this make sense?

Yes, for three reasons:

1. **The contrast scheme is statistically valid.** Each contrast vector sums to zero
   (required for any contrast), and the two `word_type` contrasts are orthogonal
   (see proof below).

2. **Nothing about the model changes.** Log-likelihood, fitted values, residuals,
   AIC/BIC, and predictions are identical. It is a pure reparametrization.

3. **It is standard practice in confirmatory research.** When hypotheses are
   pre-specified, planned contrasts make the analysis more transparent: the
   p-value for H1 and H2 can be read directly from the coefficient table.

---

## Should we re-run model selection?

No. Model reduction is driven by log-likelihood comparisons, AIC/BIC, and
convergence. Since contrast coding leaves all of these unchanged, the same
reduction decisions would be reached regardless of which coding was used.
The already-reduced models are valid — apply the contrasts directly to them.

---

## The contrast scheme

### `word_type` (3 levels: L1, L2-Remote, L2-Recent)

| Contrast | L1 | L2-Remote | L2-Recent | Role |
|---|---|---|---|---|
| Remote_vs_Recent | 0 | +1/2 | −1/2 | H1 |
| L1_vs_L2 | +2/3 | −1/3 | −1/3 | control |

### `relatedness` (2 levels: Unrelated, Related)

| Contrast | Unrelated | Related | Role |
|---|---|---|---|
| Unrelated_vs_Related | +1/2 | −1/2 | control |

With this coding, β > 0 means Unrelated > Related (e.g., longer RT or lower
accuracy for unrelated words).

---

## Proof of orthogonality

Two contrasts are orthogonal if their dot product equals zero: **Σ(c1ᵢ × c2ᵢ) = 0**.

| | L1 | L2-Remote | L2-Recent |
|---|---|---|---|
| Remote_vs_Recent (c1) | 0 | +1/2 | −1/2 |
| L1_vs_L2 (c2) | +2/3 | −1/3 | −1/3 |
| c1 × c2 | 0 | −1/6 | +1/6 |

**Σ = 0 + (−1/6) + (1/6) = 0** ✓

Orthogonality means the two contrasts carry non-overlapping information: the test
of H1 (Remote vs Recent) is statistically independent of the L1 vs L2 control.

For `relatedness` there is only one contrast (2-level factor), so orthogonality
does not apply.

---

## Coefficient table: before vs after

Both coding schemes produce **6 coefficients** — reparametrization changes what
each coefficient means, not how many there are. A 3-level factor always requires
2 columns in the model matrix regardless of coding.

### Before (default treatment coding)

| Coefficient | What it tests |
|---|---|
| Intercept | mean of L1/Related (reference level) |
| word_typeL2-Remote | L2-Remote minus L1 |
| word_typeL2-Recent | L2-Recent minus L1 |
| relatednessUnrelated | Unrelated minus Related |
| word_typeL2-Remote:relatednessUnrelated | interaction for L2-Remote, relative to L1 |
| word_typeL2-Recent:relatednessUnrelated | interaction for L2-Recent, relative to L1 |

H1 and H2 required `emmeans` post-hoc comparisons to be extracted.

### After (planned contrasts)

| Coefficient | What it tests |
|---|---|
| Intercept | grand mean |
| Remote_vs_Recent | L2-Remote minus L2-Recent → **H1** |
| L1_vs_L2 | L1 minus average of both L2s → control |
| Unrelated_vs_Related | Unrelated minus Related → control |
| Remote_vs_Recent:Unrelated_vs_Related | difference in priming between L2-Remote and L2-Recent → **H2** |
| L1_vs_L2:Unrelated_vs_Related | difference in priming between L1 and L2s → control |

H1 and H2 are readable directly from the coefficient table.
