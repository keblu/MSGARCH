Twenty-five fixes to the layer around the likelihood — inference, persistence, forecasting,
time-index handling, input validation and numerical robustness — with regression tests for
each. 18 commits, one per defect or coherent group, so any of them can be reviewed or dropped
on its own.

**The likelihood itself is untouched.** On the SMI example the log-likelihood at a fixed
parameter vector is bit-identical to `master` (`-3389.2962430913985`), the estimates agree to
eight significant figures, and AIC/BIC for an unconstrained model are unchanged. What changes
is what `summary()` prints, what `predict()` returns for a Bayesian fit, whether a saved fit
can be reloaded, and what happens at the edges.

`R CMD check --as-cran` on the branch: **1 WARNING, 1 NOTE**, both environmental (an Apple
clang warning raised inside R's own `R_ext/Boolean.h`, and "unable to verify current time"
offline). `master` reports 2 WARNINGs and 1 NOTE; the version bump clears one WARNING and
removing dead code clears the S3 NOTE.

---

## 1. Wrong numbers reported to the user

### Standard errors used the delta method transposed — `R/Inference.R:38`

```r
mSandwitch <- t(mJacob) %*% mInvHessian %*% mJacob   # -> mJacob %*% mInvHessian %*% t(mJacob)
```

`numDeriv::jacobian` returns `∂f_i/∂x_j`, so `Var(g(θ̂)) = J V J'`. The two orientations agree
only if `J` is symmetric, and it is not: the working→natural map is triangular inside each
regime (the sGARCH bound on `beta` is `0.9999 − alpha1`; the gjr/tGARCH bounds also involve
`alpha2` and the shape/skew parameters), and the transition-probability block is *anti*-diagonal.
Checked against the observed information computed directly in the natural parameterisation:

| | natural-scale `H` | this PR | `master` |
|---|---|---|---|
| `alpha1_1` | 0.01606 | 0.015111 | **0.034126** |
| `beta_1`   | 0.02197 | 0.020912 | **0.009584** |
| `alpha1_2` | 0.00431 | 0.004373 | **0.006103** |
| `beta_2`   | 0.00416 | 0.004256 | **0.000492** |
| `P_1_1`    | 0.00982 | 0.009727 | **0.000595** |
| `P_2_1`    | (at bound) | 0.030683 | **0.502501** |

Six of eight were wrong, by factors from 0.06× to 16×, and because the P block is anti-diagonal
the two transition probabilities had their standard errors **exchanged**. Only `alpha0_k`, whose
map is a plain `exp`, was unaffected.

### `Pr(>|t|)` was one-sided, then cancelled in the tail — `R/Inference.R:43`

`1 - pnorm(abs(t))` under a two-sided label. Now `2 * pnorm(-abs(t))`, evaluated on the lower
tail: `2 * (1 - pnorm(abs(t)))` cancels to exactly zero once `|t|` exceeds about 8.3, while the
lower tail stays representable to about 38. The single-regime GARCH fit used in the tests has a
`t` of 26.9, whose p-value was reported as 0 and is now 1.6e-159.

### `predict()` on an MCMC fit returned draw #1 — `R/CondVol.R:32`

`vol` is `(T+1) × ndraw` and the one-step value was `vol[dim(PredProb)[1]]` — a single index
into a matrix is linear indexing, i.e. the last row of the *first column*. On a 100-draw chain
fitted to SMI the reported value was 1.019753 (the first draw) against a posterior mean of
1.039776, with a spread of 0.999837–1.097337 across draws. `Volatility()` already averaged
correctly, so the two methods disagreed on the same fit.

### `AIC`/`BIC` mis-counted `regime.const.pars` — `R/Utils.R:420`

A regime-constant parameter removes `K − 1` degrees of freedom, not one (cf. `f_rename_par`,
which strips `name_2 … name_K`). Right only at `K = 2`: with a regime-constant shape parameter
the df was 17 instead of 16 at `K = 3` and 27 instead of 25 at `K = 4`, always over-penalising
the constrained model.

---

## 2. Failures and silent data corruption

### Saved specs and fits could not be reloaded — `R/Utils.R:349-353`

R serializes external pointers as `NULL`, so a `saveRDS`-ed spec or fit comes back with dead
pointers. `f_check_spec` exists to rebuild them, but the rebuild branch called
`spec$rcpp.func$get_mean()` / `get_sd()` — the very pointer whose failure had just triggered it:

```
Volatility(fit) : Error in .External(...): NULL value passed as symbol address
```

Those two values are already held on the R side in `spec$prior.mean` / `spec$prior.sd`, which
the next two lines were using anyway, so the C++ round trip was dead code. This is the ordinary
workflow of fitting a model, saving it, and analysing it later — or sending a spec to a
`parLapply` worker.

### `ts`/`zoo` methods recycled observations when `newdata` was supplied

`Volatility`, `predict`, `PIT`, `PredPdf` and `Risk` concatenated `object$data` with `newdata`,
then derived the index from the *already-concatenated* series and appended another
`length(newdata)` points. The index was longer than the values, so `zooreg()` recycled: a
200-point monthly series plus two new returns produced **204** observations ending in two values
copied from the start of the sample, and the model conditioned on them. The same code discarded
the original `start` and `frequency` even when `newdata` was `NULL` (a monthly series became
annual), and forecast indexes advanced by one index unit rather than by the series' own step.
Three helpers in `Utils.R` — `f_combine_data`, `f_future_index`, `f_index_result` — now handle
this in one place, replacing ten duplicated blocks. `ts` input now gives numerically identical
results to the equivalent numeric input.

### A singular stationary solve destroyed the whole fit — `src/MSgarch.h:283`

`loadparam` obtained the stationary distribution with a plain matrix inverse, on every
likelihood evaluation. The plain parameter mapping bounds each free transition probability
separately, so with `K ≥ 3` a row can leave the simplex and the matrix is singular; the uncaught
Armadillo exception then killed the run. `FitMCMC` always uses that mapping, so every chain with
three or more regimes was exposed, as was any ML fit using `fixed` or `regime.const` parameters
— `FitML` on a `K = 3` spec with `regime.const = "nu"` died deterministically at
`P = [[1, 0.999973, −0.999973], [0,1,0], [0,0,1]]`. The solve is now guarded and falls back to
the uniform distribution, which `calc_prior` rejects a moment later anyway.

### `do.sort` silently broke `constraint.spec$fixed`

The identification sort relabels regimes by unconditional variance, which moves a parameter
fixed in one regime into another. With the default `do.sort = TRUE`, a parameter fixed at 0.8
came back taking **sixteen different values** across a 100-draw chain. The sort is now skipped
when parameters are fixed, with a message.

### `FitMCMC` mis-mapped a user `ctr$par0` under `fixed.pars`

Unlike `FitML`, `FitMCMC` never dropped the fixed entries, so the vector handed to the sampler
was too long, `f_rename_par` left a trailing `NA` name, and `f_mapPar` looked bounds up by name
— mapping `alpha0_2` to 80.008 and the last parameter to `NA`.

### `constraint.spec$fixed` accepted transition probabilities that never worked

`CreateSpec` validates against `out$label`, which includes `P_1_1`, but the starting-value
routine hands `"P_1_1"` to a single-regime `CreateSpec` that has no such parameter, and the
prior correction in `Kernel()` indexes `prior.mean`, which covers only the within-regime
coefficients — so the log-posterior became `NA` and was floored to `-1e10`. Now refused with an
explanation.

---

## 3. Numerical robustness

### The Hamilton filter overflowed to `NaN` — `src/MSgarch.h`

Both filters shifted each column of regime log densities by its *smallest* entry, and only when
that fell below `log(DBL_MIN)`. That guards underflow but creates an overflow: the largest
exponent becomes `max − min − 707`, so once the regimes differ by more than about 1400 in log
density the exponential returns `Inf` and normalising gives `NaN`. The gap is reachable with
parameters the package itself accepts — `alpha0 = 3.3e-4` against `alpha0 = 1` is enough — after
which the likelihood floors to `-1e10` and `State()`, `PredPdf()` and `Risk()` return `NaN` or
fail. Both now use the standard log-sum-exp anchor, the largest entry.

Verified against the analytic limit rather than merely "no longer `NaN`": as one regime's density
vanishes the likelihood converges to `log P(other regime) + log f(y)`.

### The multi-regime in-sample CDF wrote the first slice transposed — `src/MSgarch.h:497`

`f_cdf_its` filled the `t = 0` slice of `arma::cube tmp(ny, nx, K)` with `tmp(ix, 0, s)` where
everything else in the file uses `(0, ix, s)`. With a grid shorter than the sample the first
observation's CDF was zero for every point but the first; with a longer grid the write ran off
the cube. The new test checks the first row against a closed form — at `t = 1` each regime's
conditional variance is its unconditional variance and the predictive state distribution is the
ergodic one — and it now matches exactly.

### The native log branches returned the last regime — `src/MSgarch.h:399,474`

`f_pdf` and `f_cdf` accumulate the state-weighted mixture in `out` but under `is_log` overwrote
it with `log(tmp[i])`, where `tmp` holds only the regime evaluated last: `exp()` of the returned
value differed from the mixture by factors of 12 and 26 in the tails. Not reachable from
`PredPdf`/`PIT`, which pass `is_log = FALSE`, but both are exposed as module methods on
`spec$rcpp.func`.

---

## 4. Validation and diagnostics

- **`FitML` reported every failure identically.** The guard tested `llk == 1e+10`, but `f_nll`
  returns `+1e10` so a failure arrives as `-1e10`; `f_OptimFUNDefault` also wraps `optim` in
  `try()`, so `optimizer$value` errored first. Everything — a singular transition matrix, a
  malformed spec, one `NA` in the data — surfaced as
  `$ operator is invalid for atomic vectors`.
- **Data and parameters.** `f_check_y`/`f_check_par` rejected input only when it was *entirely*
  `NaN`, so a single `NA`, `NaN` or `Inf` reached the compiled code and `Volatility()` returned
  a complete, plausible-looking series computed from corrupt data. Both now require finite
  values, and the data must hold at least two observations. `f_nll`/`f_posterior` treat a
  non-finite mapped parameter vector as an infeasible point so strictness cannot abort a fit.
- **`Risk()`.** `alpha`, `nahead` and `ctr$nmesh` are validated (`alpha = 1.5` used to return
  6.6745 silently; `alpha = 0` gave `ES = -Inf`). The evaluation grid spans the observed data
  range, so it can miss part of the predictive distribution; when the omitted mass exceeds the
  requested tail probability, the "quantile" is the grid boundary, and that now warns. On a
  normal fit to SMI the grid's left endpoint has CDF 7e-8 and nothing changes; on a four-point
  sample it covers 15% of the distribution and warns.
- **`prior$sd`** was validated by the *mean* checker; `sd = 0` and `sd = -1` were accepted and
  reached C++. Now checked for finiteness and positivity, and both prior validators produce
  usable messages instead of `stop(cat(...))`.
- **`CreateSpec`** validates `switch.spec$K`, and when expanding one regime through `K` rejects
  an explicitly heterogeneous distribution vector — the guard tested
  `distribution.spec$model`, which does not exist, so `model = "sGARCH"`,
  `distribution = c("norm","std")`, `K = 3` silently produced `norm, std, norm`.
- **`simulate()`** accepts `nburn = 0`, which used to drop the first draw (`1:0` is `c(1, 0)`)
  and then fail on the dimnames.
- **`UncVol()`** averages the horizons *after* the burn-in; `nburn:nahead` is a descending range
  under the shipped defaults and in general keeps part of the transient.

---

## 5. Housekeeping

- The `Sim.MSGARCH_ML_FIT` / `Sim.MSGARCH_MCMC_FIT` methods behind the standing `R CMD check`
  S3 NOTE are unreachable — `Sim` is not exported, every internal caller passes a spec, and
  `simulate.*_FIT` pass `object$spec` explicitly. Removed; the NOTE is gone.
- `Rcpp:::LdFlags()` dropped from both `Makevars`; Rcpp has not needed it since 2013.
- The shipped BIC test asserted `abs(exp.BIC - exp.BIC) < tol`, which is zero by construction.
  (The value it meant to check was correct, so this is a test fix, not a bug fix.)
- The unused log-likelihood accumulator in `f_get_Pstate` is removed.
- The mixture transition matrix built for Viterbi decoding is now row-stochastic. **The decoded
  path is unchanged** — the misplaced factor is constant in the index being maximised over, so
  it cancels — but the matrix should not have rows summing to 1.6 and 0.4.

---

## Tests

Five new files and additions to two existing ones: 28 blocks and 103 assertions. The suite
goes from 23 blocks / 23 assertions to 52 blocks / 128 assertions, adding about 25s of
check time.

- `test_Inference.R` — pins the delta method against an independently computed observed
  information. Anchored on a single-regime GARCH(1,1)-Normal because every parameter there is
  interior, so a central-difference Hessian of the natural-scale negative log-likelihood is well
  conditioned: it agrees with `J V J'` to 1e-4 relative while the transposed sandwich is off by
  44% and 85%. Also the two-sided and tail-precision checks, and the degrees-of-freedom
  arithmetic for `K = 2, 3, 4`.
- `test_Serialization.R` — round-trips a spec, an ML fit and an MCMC fit through
  `saveRDS`/`readRDS` and requires identical results with priors preserved.
- `test_NativeDensity.R` — the in-sample CDF against a closed form, a grid longer than the
  sample, and the log branches against `log()` of the linear ones.
- `test_Robustness.R` — the filter under a log-density gap over 1000 including the analytic
  limit; a singular transition matrix; non-finite and too-short input; `nburn = 0`; the `UncVol`
  window; the mixture Viterbi path against per-observation MAP decoding; `Risk`'s guards and
  warning; prior validation; and that a fixed parameter survives `FitML`, a user `par0` and the
  identification sort.
- `test_TimeSeriesIndex.R` — index/value alignment, `ts` and `zoo` against numeric, and forecast
  index spacing.

Every block opens with a guard asserting its own precondition — the two sandwich orientations
really differ for this model, the round trip really did invalidate the pointers, the posterior
mean really differs from the first draw — so none of them can pass vacuously if the surrounding
code changes. All fail on the commit preceding their fix.

---

## Compatibility

Deliberate changes in output, all of them corrections:

- Standard errors and p-values from `summary()` change for every model with a `beta` or a
  transition probability.
- `predict()` on an MCMC fit returns the posterior mean rather than the first draw.
- `AIC`/`BIC` change for `K ≥ 3` models using `regime.const`.
- `ts` and `zoo` users get results of the correct length, with the correct index, from
  `Volatility`, `predict`, `PIT`, `PredPdf` and `Risk`.
- Models using `constraint.spec$fixed` are no longer regime-sorted, so regime labelling may
  differ from previous fits.

Newly rejected input that used to be accepted: data containing `NA`, `NaN` or `Inf`; data
shorter than two observations; non-finite parameters; `alpha` outside `(0,1)`; `ctr$nmesh < 2`;
non-positive `prior$sd`; `constraint.spec$fixed` on a transition probability; and a
heterogeneous distribution vector combined with `switch.spec$K`.

Because inference and prediction results change, the reverse imports (`MSGARCHelm`, `SBAGM`)
are worth a check before release.

The version bump to 2.52 and the `Date` refresh are in their own commit (`c93e218`) — drop it if
the release number should be set separately; nothing else depends on it.

`REVIEW.md` and `REVIEW_codex.md` in the diff are the working notes behind these fixes: a
defect-by-defect write-up with reproductions, and the verbatim output of an independent
read-only audit used to cross-check it. Happy to drop them from the branch if you would rather
they not live in the repository.
