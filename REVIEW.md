# Review of the `MSGARCH` package (v2.51)

Reviewed 2026-08-10 on R 4.5.2 (aarch64-apple-darwin20), Apple clang 21.
Source: `package MSGARCH/Package/`, git `17017fa` (2022-12-05, maintainer K. Bluteau,
David = aut). Fifth in the R-package sweep after `AdMit`, `bayesGARCH`, `DEoptim`, `GAS`.

**Method.** (i) `R CMD check --as-cran` on a freshly built tarball; (ii) line-by-line read of
all 31 `.R` files and all 30 `src/` files; (iii) analytic re-derivation of every conditional
density, CDF, quantile function and truncated moment in `Normal.h` / `Student.h` / `Ged.h` /
`Symmetric.h` / `Skewed.h`, and of the stationarity conditions in `sGARCH.h` / `gjrGARCH.h` /
`tGARCH.h` / `eGARCH.h` / `sARCH.h`; (iv) numerical reproduction of every defect below on
`data("SMI")` with the shipped binary. Everything reported here is reproduced, not inferred.

> **Status (2026-08-10).** Items 1–4 of the work list at the bottom — **A1, A2, A3, A4 and
> B1** — are **fixed** in `Package/` (7 changed lines across `R/Inference.R`, `R/Utils.R`,
> `R/CondVol.R`). Each fix is verified below in the relevant section, `R CMD check --as-cran`
> is byte-for-byte unchanged (2 WARNINGs, 1 NOTE; `testthat` OK), and all eight specifications
> in a regression sweep (SR / MS-2 / MS-3 / mixture / heterogeneous / skew-t / GED / both
> constraint paths) fit with finite, correct standard errors. Regression cover was added for
> all five — `tests/testthat/test_Inference.R` (new), `tests/testthat/test_Serialization.R`
> (new) and one block appended to `tests/testthat/test_Volatility.R`: 11 blocks / 34
> assertions, all passing on the patched build and all 11 failing on stock 2.51.
> `DESCRIPTION`'s version and `NEWS` are untouched — that is a release decision for the
> maintainer. Everything from **B2** onwards is still open.

> **Independent cross-check (Codex).** The whole of this document — the 14 claims, the five
> fixes and the tests — was re-audited read-only by `codex-cli 0.139.0` working from source
> alone (`codex exec --sandbox read-only`); its verbatim output is in `REVIEW_codex.md`.
> It returned **CONFIRM on 13 of 14 claims, PARTIALLY CONFIRM on B3**, found no error in the
> parts this review calls correct, and marked all five fixes OK. Its three pushbacks:
>
> 1. **B3's "trailing NA name" mechanism is overstated** — Codex expected `f_rename_par` to
>    error on a too-long vector, or the fixed dimension to simply stay in the sampler.
>    **Not upheld.** `names(x) <- <shorter>` pads with `NA` in R rather than erroring, and the
>    trace is reproducible: `f_rename_par` returns
>    `alpha0_1, alpha1_1, alpha0_2, alpha1_2, beta_2, P_1_1, P_2_1, NA` and `f_mapPar` then
>    yields `0.1, 0.1, 80.008, 0.001, 0.1, 0.8001, 0.5, NA`. Codex's alternative also fails:
>    in the working (no-`par0`) path the fixed parameter *is* respected — `beta_1` is exactly
>    0.8 across all 300 draws. B3 stands as written.
> 2. **The AIC/BIC test's only real fit is K = 2, where both dof formulas coincide** — fair,
>    and already stated in that section; the K = 3/4 legs are deliberately synthetic because
>    B2 makes a constrained K ≥ 3 fit impossible.
> 3. **"All 11 blocks fail on stock 2.51" is not literally demonstrated** — it is, by
>    measurement rather than by reading: the AIC/BIC block reports `pass=3 fail=6` on stock.
>    Wording sharpened below.

**`R CMD check --as-cran`: 2 WARNINGs, 1 NOTE — none of them is any of the defects below.**
The WARNINGs are the usual pair (CRAN incoming feasibility: version not bumped + `Date`
over a month old; and one spurious clang warning from R's own `R_ext/Boolean.h`). The NOTE is
a genuine but harmless S3 signature mismatch on `Sim`/`Sim.MSGARCH_ML_FIT`. `tests/testthat/`
(9 files as shipped) passes and covers none of this.

**The conditional distributions are clean.** Unlike `GAS`, every density/CDF/quantile/moment
checked out analytically: the standardised Student-*t* and GED constants, `E|z|` for all three
families, the Fernández–Steel construction in `Skewed.h` (kernel, CDF, inverse-CDF and all four
truncated moments `Eabsz`, `EzIpos`, `EzIneg`, `Ez2Ineg` — I re-derived each and they match,
including both `xi >= 1` and `xi < 1` branches), the Hamilton filter and its over/underflow
bookkeeping, Kim's smoother, the ergodic-distribution formula, the transition-matrix memory
layout, and the tGARCH/gjrGARCH second-moment stationarity conditions. **The likelihood itself
is correct.** The defects are in the layer around it: inference, persistence, forecasting for
Bayesian fits, and the constrained-estimation paths.

---

## Tier 1 — wrong numbers reported to the user

### A1. Every standard error is computed with the delta method transposed — `R/Inference.R:38` — **FIXED**

```r
mJacob     <- numDeriv::jacobian(f_mapPar, vPw_mod, spec = spec, do.plm = do.plm)
mSandwitch <- t(mJacob) %*% mInvHessian %*% mJacob      # <- should be mJacob %*% ... %*% t(mJacob)
```

`numDeriv::jacobian(f, x)[i, j] = ∂f_i/∂x_j`, so with `θ_natural = g(θ_working)` the delta
method is `Var(g) = J V J'`. The package computes `J' V J`. That is identical only when `J` is
symmetric, and here it is not: the working→natural map is triangular within each regime
(the sGARCH upper bound on `beta` is `0.9999 − alpha1`, gjr/tGARCH bounds on `beta` depend on
`alpha1`, `alpha2` **and** on the shape/skew parameters), and the transition-probability block
is *anti*-diagonal, because `f_mapGamma`/`f_unmapGamma`
(`R/ParameterTransformation.R:180-213`) enumerate the off-diagonal entries in column-major
order but attach the row-major parameter names to them. The Jacobian at the default
MS(2)-GARCH(1,1)-Normal fit on `SMI`:

```
         alpha0_1 alpha1_1  beta_1 alpha0_2 alpha1_2 beta_2   P_1_1   P_2_1
alpha1_1        0   0.0795  0.0000        0        0      0       0       0
beta_1          0  -0.0767  0.0303        0        0      0       0       0
beta_2          0        0  0.0000        0  -0.0054  5e-04       0       0
P_1_1           0        0  0.0000        0        0      0  0.0000 -0.0212
P_2_1           0        0  0.0000        0        0      0  0.0013  0.0000
```

**Evidence.** I computed the observed information directly in the natural parameterisation
(central differences on `Kernel(spec, par, y, log = TRUE, do.prior = FALSE)`, which reproduces
`fit$loglik` to machine precision) and compared:

| | natural-scale `H` | delta method `J V J'` | **package `J' V J`** | package / correct |
|---|---|---|---|---|
| `alpha0_1` | 0.00744 | 0.00725 | 0.00725 | 1.00 |
| `alpha1_1` | 0.01606 | 0.01511 | **0.03413** | **2.26×** |
| `beta_1`   | 0.02197 | 0.02091 | **0.00958** | **0.46×** |
| `alpha0_2` | 0.01688 | 0.01769 | 0.01769 | 1.00 |
| `alpha1_2` | 0.00431 | 0.00437 | **0.00610** | **1.40×** |
| `beta_2`   | 0.00416 | 0.00426 | **0.00049** | **0.12×** |
| `P_1_1`    | 0.00982 | 0.00973 | **0.00059** | **0.06×** |
| `P_2_1`    | (boundary, 0.9987) | 0.03072 | **0.50250** | **16×** |

The natural-scale information agrees with `J V J'` to 2–6% on the seven interior parameters
(`P_2_1` sits at 0.9987 where the natural-scale Hessian is not usable). The shipped numbers
are wrong on **six of eight** parameters, by factors from 0.06× to 16×. Because the P-block of
`J` is anti-diagonal, transposing literally **swaps the two transition probabilities' standard
errors**: `summary(fit)` reports `P_1_1 = 0.978 (s.e. 0.0006)` and
`P_2_1 = 0.9987 (s.e. 0.503)` when the correct values are ≈ 0.0097 and ≈ 0.031.

Only `alpha0_k`, whose map is a plain `exp`, is unaffected. Every published table produced by
`summary()` on a model with a `beta` or a transition probability is affected. **Fix**: one
transpose. Verified.

**After the fix** (`mSandwitch <- mJacob %*% mInvHessian %*% t(mJacob)`), the same fit:

```
         Estimate  Std. Error   t value  Pr(>|t|)   natural-scale H
alpha0_1 0.021632    0.007245     2.986  0.002830          0.007444
alpha1_1 0.087024    0.015111     5.759  0.000000          0.016063
beta_1   0.881494    0.020912    42.152  0.000000          0.021975
alpha0_2 0.020660    0.017687     1.168  0.242771          0.016884
alpha1_2 0.005396    0.004373     1.234  0.217243          0.004309
beta_2   0.994041    0.004256   233.588  0.000000          0.004161
P_1_1    0.978348    0.009727   100.577  0.000000          0.009823
P_2_1    0.998703    0.030720    32.509  0.000000    (boundary)
```

and across a sweep of eight specifications (single-regime; MS-2 sGARCH-Normal; MS-2
gjrGARCH/tGARCH with sstd/std; MS-2 eGARCH-GED; MIX-2 sGARCH-std; MS-3; `fixed.pars`;
`regime.const.pars`) every fit returns finite positive standard errors that reproduce
`J V J'` to 1e-8 and p-values inside [0, 1].

### A2. `Pr(>|t|)` is a one-sided p-value — `R/Inference.R:43` — **FIXED**

```r
vPvalues <- 1 - pnorm(abs(vTest))          # should be 2 * (1 - pnorm(abs(vTest)))
```

The column is labelled `Pr(>|t|)` and printed by `summary()`. Reported values are exactly half
the two-sided p-value: on the `SMI` fit, `alpha1_1` shows `0.00539` where the two-sided value
is `0.01077`, and `P_2_1` shows `0.0234` vs `0.0469`. Combined with A1 the two errors do not
cancel — they compound.

(Separately: t-tests against zero are meaningless for `beta`, for `nu`, and for the transition
probabilities, all of which have non-zero-centred supports. A note in `?FitML` would help.)

**After the fix**: the `Pr(>|t|)` column equals `2 * (1 - pnorm(abs(t)))` exactly
(`all.equal` TRUE) on every specification in the sweep.

### A3. `predict()` on an MCMC fit returns draw #1, not the posterior mean — `R/CondVol.R:32` — **FIXED**

```r
vol <- matrix(NA, nrow = dim(PredProb)[1], ncol = nrow(par.check))
...
tmp <- mean(vol[dim(PredProb)[1]])          # linear index -> row T+1 of COLUMN 1 only
```

`vol` is `(T+1) × ndraw`. A single index into a matrix is linear indexing, so `vol[T+1]` is the
last row of the *first* column. `mean()` of a scalar is a no-op. The intended expression is
`mean(vol[dim(PredProb)[1], ])`.

Reproduced on `FitMCMC(CreateSpec(), SMI, ctr = list(nburn=500, nmcmc=1000, nthin=10))`,
100 retained draws:

```
predict(mc, nahead = 1)$vol : 1.019753
draw #1 only                : 1.019753   <- exact match
posterior mean over draws   : 1.039776
range across draws          : 0.999837 – 1.097337
```

The reported one-step-ahead volatility is whichever value the first retained draw happens to
give. `Volatility()` (the in-sample path) averages across draws correctly, so the two functions
are mutually inconsistent on the same fit. This also feeds the `h = 1` row of
`predict(..., nahead = h)` and of `UncVol()` for MCMC fits.

**After the fix**, same chain: `predict(mc, nahead = 1)$vol = 1.039776`, exactly the posterior
mean, and no longer draw #1's `1.019753`. The ML path is unchanged (one column, so
`mean(vol[N, ])` and the old `vol[N]` coincide): `predict(fit, nahead = 1)$vol = 1.030426`
before and after.

### A4. `AIC`/`BIC` degrees of freedom undercount `regime.const.pars` — `R/Utils.R:420` — **FIXED**

```r
dofMSGARCH = function(object){
  length(object$spec$par0) - length(object$spec[["regime.const.pars"]]) - length(object$spec[["fixed.pars"]])
}
```

Each regime-constant parameter removes `K − 1` free parameters, not one (`f_rename_par` strips
`name_2 … name_K`). Correct only at `K = 2`:

```
K=2  par0=10  truly free= 9  dofMSGARCH= 9   ok
K=3  par0=18  truly free=16  dofMSGARCH=17   MISMATCH
K=4  par0=28  truly free=25  dofMSGARCH=27   MISMATCH
```

`stats::AIC`/`BIC` use this via `logLik.MSGARCH_ML_FIT`, so every K ≥ 3 model-selection table
built with `regime.const` is penalised wrongly — and always in the direction that *disfavours*
the constrained model. Fix: `- length(regime.const.pars) * (K - 1)`.

**After the fix**: `K=2 -> 9`, `K=3 -> 16`, `K=4 -> 25`, all matching the true free-parameter
count. The unconstrained (`8`) and `fixed.pars` (`7`) cases are unchanged.

---

## Tier 2 — hard failures

### B1. A saved spec or fit cannot be reloaded, and the recovery path that exists for exactly this is dead code — `R/Utils.R:340-359` — **FIXED**

`MSGARCH_SPEC` holds Rcpp module objects, so `saveRDS`/`readRDS` across sessions leaves a stale
external pointer. `f_check_spec` exists to detect that and rebuild — but the rebuild branch
re-dereferences the *same dead pointer* on its second line:

```r
is.OK = tryCatch({ spec$rcpp.func$get_sd(); TRUE }, error = function(e) FALSE)
if (!isTRUE(is.OK)) {
  spec.new   = f_spec(models = spec$name, do.mix = spec$is.mix)
  prior.mean = spec$rcpp.func$get_mean()     # <- line 349: dead pointer again, uncaught
  prior.sd   = spec$rcpp.func$get_sd()       # <- line 350: same
  ...
  spec$rcpp.func$set_mean(spec$prior.mean)   # <- and then uses spec$prior.mean, not the
  spec$rcpp.func$set_sd(spec$prior.sd)       #    prior.mean/prior.sd just computed
}
```

In a fresh session, on a `fit` written with `saveRDS`:

```
Volatility(fit)   : Error in .External(...): NULL value passed as symbol address
State(fit)        : Error in .External(...): NULL value passed as symbol address
predict(fit)      : Error in .External(...): NULL value passed as symbol address
FitML(saved spec) : Error in .External(...): NULL value passed as symbol address
```

**Fix applied**: drop lines 349–353 (they recompute, from the dead pointer, exactly the named
vectors already stored in `spec$prior.mean` / `spec$prior.sd`, and then throw the result away),
and keep lines 355–356, which already read the R-side copies. In a fresh session on a
`readRDS`-ed fit, `Volatility`, `State`, `predict`, `summary`, `Risk`, `PIT`, `simulate` and
`AIC` all now succeed, and return values identical to the in-session ones (`Volatility` head
`1.2045, 1.2530, 1.2372`; `predict` h=1 `1.030426` both ways). Custom priors survive the
rebuild — a spec created with `prior = list(mean = list(beta_1 = 0.7), sd = list(beta_1 = 0.1))`
comes back with mean `0.7` / sd `0.1` on `beta_1` and defaults elsewhere. Saved
`MSGARCH_MCMC_FIT` objects (`predict`, `DIC`) work too. This unblocks the single most common
workflow in applied use: fit overnight, save, analyse later; `parLapply` over a rolling window;
caching a fit in a knitr chunk.

### B2. Singular ergodic-distribution inverse aborts the whole fit — `src/MSgarch.h:277-286`

`loadparam` recomputes the stationary distribution on *every* likelihood evaluation with a raw
Armadillo inverse:

```cpp
arma::mat foo   = (I - as<arma::mat>(P_mat) + Umat).t();
arma::vec delta = (foo).i() * Uvec;          // throws if singular
```

With `do.plm = TRUE` — forced whenever `fixed.pars` or `regime.const.pars` is set, and
hard-coded in `FitMCMC` (`R/FitMCMC.R:134`) — the free transition entries are mapped into
`(0,1)` *independently*, so for `K ≥ 3` a row can sum to more than 1 and `extract_P_it` then
produces a negative last entry. `I − P + U` can be exactly singular, and the Armadillo
exception propagates out of `Kernel`, out of `f_nll`, and kills the run.

Reproduced deterministically. `CreateSpec(model="sGARCH", distribution="std", K=3,
constraint.spec=list(regime.const="nu"))` + `FitML(..., SMI)` with `set.seed(1)`; instrumenting
`f_nll` to record the last parameter point reached by BFGS gives

```
     [,1]     [,2]      [,3]
[1,]    1 0.999973 -0.999973
[2,]    0 1.000000  0.000000
[3,]    0 0.000000  1.000000
det(I - P + U) = 0    rcond = 0
Kernel  -> ERROR: matrix multiplication: problem with matrix inverse
```

`K = 3` *without* `regime.const` fits fine (the `do.plm = FALSE` map always yields a proper
stochastic matrix); `K = 2` *with* `regime.const` fits fine (one free entry per row cannot
overflow the simplex). So the failure is specific to `K ≥ 3` on the `do.plm` path — which
includes **all** `FitMCMC` runs with three or more regimes, where a single unlucky proposal
destroys the entire chain with no partial output.

Fix: guard the inverse (`arma::solve` with `arma::solve_opts::no_error`, or check
`arma::inv(...)`'s bool return) and return the model's own "infeasible" signal (`-1e10`) instead
of throwing; the constraint `all(0 < P_it < 1)` already exists in `calc_prior` but is evaluated
*after* `loadparam`, so it never gets the chance to reject the point.

### B3. `FitMCMC` mis-aligns the parameter vector when `fixed.pars` is combined with `ctr$par0` — `R/FitMCMC.R:140-151`

`FitML` removes the fixed entries from the starting vector (`f_remove_fixedpar`,
`R/FitML.R:126`); `FitMCMC` substitutes their values but **never removes them**, so `par0` is
passed to the sampler with `d` elements when the sampler expects `d − n_fixed`. `f_rename_par`
then labels the over-long vector with the short name list, leaving a trailing `NA` name, and
`f_mapPar` looks the bounds up by name — so every parameter after the fixed one is mapped with
the *wrong* bounds and the last one becomes `NA`:

```
spec: 8 parameters, beta_1 fixed -> sampler expects 7
length(par0) inside FitMCMC = 8
names from f_rename_par     : alpha0_1,alpha1_1,alpha0_2,alpha1_2,beta_2,P_1_1,P_2_1,NA
f_mapPar gives              : 0.1, 0.1, 80.008, 0.001, 0.1, 0.8001, 0.5, NA
FitMCMC(..., ctr=list(par0=spec$par0)) -> Error : matrix multiplication: problem with matrix inverse
```

`alpha0_2` becomes 80, `alpha1_2` becomes 0.001, and the `NA` then detonates the inverse of
B2. Without a user `par0` the same spec samples fine (`accept = 0.287`), so this is purely
the `ctr$par0` branch. Fix: add `par0 <- f_remove_fixedpar(par0, spec$fixed.pars)` after
line 148, mirroring `FitML`.

### B4. `constraint.spec = list(fixed = list(P_i_j = …))` is accepted but unusable

`CreateSpec` validates fixed-parameter names against `out$label`, which includes `P_1_1`, so
fixing a transition probability is accepted. Two independent things then break:

1. `f_recover_fixedpar_SR` (`R/ParameterConstraints.R:34-56`) splits fixed parameters by
   regime with `gsub("_k", "", name)`; `"P_1_1"` matches the regime-1 test, is handed to a
   *single-regime* `CreateSpec` whose labels are only `alpha0_1, alpha1_1, beta_1`, and dies:
   `Wrong name in fixed.pars: P_1_1`. So `FitML(spec, data)` and `FitMCMC(spec, data)` both
   fail before the optimiser starts. (With an explicit `ctr$par0`, `FitML` bypasses the
   starting-value routine and works — `loglik = -3390.59`.)
2. The prior correction in `Kernel` (`R/Kernel.R`) subtracts
   `dnorm(par[, names(fixed.pars)], prior.mean[names(fixed.pars)], …)`, but `prior.mean` only
   covers the `sum(NbParams)` within-regime coefficients — the transition probabilities have a
   uniform prior and no entry. `prior.mean["P_1_1"]` is `NA`, the whole log-posterior becomes
   `NA` and is floored to `-1e10`:

```
length(prior.mean) = 6, length(label) = 8   (P_1_1, P_2_1 not covered)
log-posterior  (do.prior=TRUE)  : -1e+10
log-likelihood (do.prior=FALSE) : -3487.971
```

so even reaching the sampler would give a frozen chain. Either reject `P_*` in
`f_check_parameterConstraints` with a clear message, or handle it in both places.

### B5. The `FitML` failure guard has the wrong sign and is unreachable anyway — `R/FitML.R:129-137`

```r
optimizer <- ctr$OptimFUN(vPw, f_nll, spec, data_, ctr$do.plm)
llk <- -optimizer$value
if (llk == 1e+10) { f_error("FitML -> Error during optimization"); stop() }
```

`f_nll` returns `+1e10` on failure, so `llk` is `−1e10`; the test can never be true. And
`f_OptimFUNDefault` wraps `optim` in `try()`, so on any error `optimizer` is a character
`try-error` and `optimizer$value` throws first. Every failure mode therefore surfaces as

```
Error in optimizer$value : $ operator is invalid for atomic vectors
```

— which is what a user sees for B2, for B4, and for a single `NA` in the data
(`f_check_y`, `R/Utils.R:177`, only rejects data that is *entirely* `NaN`, so
`y[100] <- NA` sails through and produces this message). Three unrelated problems, one
uninterpretable error. Fix: test `inherits(optimizer, "try-error")` first, then `llk == -1e10`;
and make `f_check_y` reject `any(!is.finite(y))`.

---

## Tier 3 — silently wrong results in specific calls

### C1. `pdf_Rcpp` / `cdf_Rcpp` with `is_log = TRUE` return the last regime, not the mixture — `src/MSgarch.h:397-401, 472-476`

```cpp
for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
  for (int i = 0; i < nx; i++) {
    tmp[i] = (*it)->spec_calc_pdf(x[i] / sig) / sig;
    out[i] = out[i] + tmp[i] * PLast[s];          // out = correct mixture
  }
  s++;
}
if (is_log) { for (int i = 0; i < nx; i++) out[i] = log(tmp[i]); }   // tmp = LAST regime only
```

`tmp` is overwritten each regime, so the mixture in `out` is thrown away and replaced by the
log density of regime *K* alone:

```
x                        -3       -1        0        1        3
mixture pdf         0.00561  0.24042  0.39670  0.24042  0.00561
returned (log=TRUE) -2.68223 -1.77355 -1.65997 -1.77355 -2.68223
exp(returned)/pdf     12.199    0.706    0.479    0.706   12.199   <- should be all 1
cdf, same ratio       26.480    1.972    1.000    0.814    0.926
```

Not reachable from `PredPdf`/`PIT` (both always pass `FALSE` and take the log in R), but these
are live methods on `spec$rcpp.func` and the flag is part of the C++ signature. Fix:
accumulate into a scratch vector and take `log(out[i])`.

### C2. `MSgarch::f_cdf_its` writes the first observation transposed — `src/MSgarch.h:497`

```cpp
tmp(ix, 0, s) = (*it)->spec_calc_cdf(x(ix, 0) / sig);   // t=0 block
...
tmp(i, ix, s) = (*it)->spec_calc_cdf(x(ix, i) / sig);   // t>=1 loop — correct orientation
```

The cube is `(ny, nx, K)`. `SingleRegime::f_cdf_its` and both `f_pdf_its` write `(0, ix, s)`;
only the multi-regime CDF got it backwards. Consequences: the `t = 1` row is left at zero for
every grid point except the first, and if `nx > T` the write runs off the cube.

```
cdf_its, t=1, regime 1: 0.007916 0 0 0 0        <- zeros
cdf_its, t=2, regime 1: 0.013265 0.133676 0.5 0.866324 0.986735
pdf_its, t=1, regime 1: 0.0262 0.232471 0.481275 0.232471 0.0262   (pdf path is fine)

PIT(spec, x = c(-2,-1,0,1,2), par, data = SMI, do.its = TRUE)[1:2, ]
1990-11-12  0.015613  0.000000  0.0  0.000000  0.000000     <- t=1 wrong
1990-11-13  0.020841  0.140055  0.5  0.859945  0.979159

PIT(..., x = <length T+10>, do.its = TRUE) -> Error: Cube::operator(): index out of bounds
```

Harmless in the common `PIT(fit, do.its = TRUE)` call (there `nx = 1`), wrong whenever a user
supplies an evaluation grid.

### C3. `UncVol` averages the burn-in window and ignores the converged tail — `R/UncVol.R:88` + `R/Utils.R:160`

The documentation says `nahead = 5000L`, `nburn = 1000L`, "simulating nsim paths up to
`nburn + nahead` … discarding the first `nburn` … and computing the mean of the remaining".
The code has the two defaults **swapped** (`f_process_ctr(type = 2)`:
`nburn = 5000L, nahead = 1000L`) and averages

```r
out <- mean(tmp[ctr$nburn:ctr$nahead])        # should be mean(tmp[(nburn+1):(nburn+nahead)])
```

With the shipped defaults that is `tmp[5000:1000]` — a *descending* sequence covering horizons
1000–5000 out of the 6000 simulated, i.e. the best-converged 1000 horizons are the ones thrown
away. It happens to land far enough out that the default answer is roughly right; it stops
being right the moment a user passes their own control list:

```
UncVol(fit, ctr = list(nsim=100, nburn=400, nahead=100))
  what the code averages   (horizons 100..400) : 1.36344
  what the docs promise    (horizons 401..500) : 1.06793      <- 28% apart
two calls, identical settings                  : 1.36344 / 1.33357
```

The last line is a second issue: `UncVol` is a Monte Carlo estimate with no `seed` argument, so
it is not reproducible, and it is driven from a two-point fake sample `data = c(1, 1)`
(`R/UncVol.R:84`). Worth documenting at minimum.

### C4. `State()$Viterbi` uses a transposed transition matrix for mixtures — `R/State.R:40`

```r
P <- TransMat(object, par = par[i, ], nahead = 1)          # 1 x K row of mixture weights
if (isTRUE(object$is.mix)) P <- matrix(rep(P, object$K), nrow = object$K, ncol = object$K)
```

`matrix()` fills column-major, so row *i* becomes `(p_i, …, p_i)` instead of `(p_1, …, p_K)`.
With mixture weights `(0.8, 0.2)`:

```
passed to Viterbi():   [,1] [,2]        correct:   [,1] [,2]
                 [1,]   0.8  0.8                    0.8  0.2
                 [2,]   0.2  0.2                    0.8  0.2
row sums: 1.6, 0.4   <- not a stochastic matrix
```

The decoded path for every `do.mix = TRUE` model is computed against this. Fix: `byrow = TRUE`.
(Filtered/smoothed/predicted probabilities are computed in C++ from the correct matrix and are
unaffected — only `$Viterbi`.)

### C5. `CreateSpec`'s `K`-expansion guard reads a field that does not exist — `R/CreateSpec.R:213`

```r
if (length(variance.spec$model) > 1 | length(distribution.spec$model) > 1)
#                                                             ^^^^^ should be $distribution
```

`distribution.spec$model` is always `NULL`, so `length(NULL) > 1` is always `FALSE` and the
check never fires for the distribution vector. The subsequent `rep(...)` then recycles and
`distribution[1:length(model)]` truncates:

```r
CreateSpec(variance.spec = list(model = "sGARCH"),
           distribution.spec = list(distribution = c("std","norm")),
           switch.spec = list(K = 3))
$name  "sGARCH_std"  "sGARCH_norm"  "sGARCH_std"
```

The user asked for something contradictory and silently got a scrambled three-regime spec
instead of the intended error.

---

## Tier 4 — minor, hygiene, and dead code

- **`R/Utils.R:294` `f_getGamma`** scrambles the transition matrix for `K ≥ 3` (it reads the
  row-major parameter vector with a column-major `matrix()` fill). Currently **unused** — but
  it is a landmine for anyone who reaches for it.
- **`R/Utils.R:2` `f_GammaParNames`** emits `P_1_1, P_2_1, P_3_1, P_1_2, …` (column-major),
  whereas the labels actually built by `f_spec:398-401` and consumed by
  `MSgarch::extract_P_it` are `P_1_1, P_1_2, P_2_1, …` (row-major). Also unused; delete or fix,
  don't leave two contradictory definitions of the parameter order in the same file.
- **`src/Utils.h:42`** `adaptiveSimpsonsAux` guards recursion with `bottom <= 100` instead of
  `bottom <= 0`, so any call with `maxRecursion ≤ 100` returns after a single Simpson step with
  no adaptation. Unreachable — `adaptiveSimpsons` is never called — but it should either be
  fixed or removed. (The integrator actually in use, `Skewed::compositeSimpsons` with
  `Nsi = 5`, is fine: 10 panels over a short interval on a smooth integrand.)
- **`src/MSgarch.h:71`** `P_mean = 1 / K` is integer division on an `int` — always 0 for
  `K ≥ 2`. Both `P_mean` and `P_sd` are write-only; delete them.
- **`src/MSgarch.h:768,774`** `eval_model` calls `calc_prior(theta_j)` twice per parameter
  vector and discards the first result. Cost is small next to the filter (measured ~0.8 ms per
  kernel evaluation at `T = 2500`), but it is free to remove.
- **`src/adaptMCMC.cpp:82`** a full `eig_sym` is computed and immediately discarded every MCMC
  iteration (`makePositiveDefinite` recomputes it on line 83). Also `sigma.eye(l_param,
  l_param)` on line 81 zeroes `sigma` *in place* mid-expression; it happens to be harmless
  because `S` was cached on line 77, but it is a trap for the next editor.
- **`R/CreateSpec.R:312`** validates `prior.sd` with `f_check_parameterPriorMean`, so a bad
  name in `prior$sd` reports "Wrong name in prior.mean".
- **`R/Utils.R:177`** `f_check_y` only rejects data that is entirely `NaN`; a single `NA`
  passes and surfaces later as the B5 error.
- **`logLik.MSGARCH_ML_FIT`** sets `nobs = length(data)`, but `calc_lndMat` evaluates `T − 1`
  terms (observation 1 initialises the recursion). `BIC` therefore uses `log(T)` rather than
  `log(T − 1)` — negligible, but worth a line in `?FitML` since it also explains the
  `State()$Viterbi[1] = Viterbi[2]` fudge at `R/State.R:99`.
- **`eGARCH::set_vol`** initialises `h = exp(α₀/(1−β))`, i.e. `exp(E[ln h])`, and
  **`tGARCH::set_vol`** initialises `h = (E[σ])²`, not `E[h]`. Defensible as a recursion
  starting value, but `UncVol`'s per-regime `unc_vol_Rcpp` inherits the same convention and
  therefore reports something other than the unconditional variance for those two models.
  Document it.
- **Dead/unused**: `f_map_deriv`, `f_rbindrep`, `f_getGamma`, `f_GammaParNames`,
  `Decoding_HMM` (exported to R, called nowhere), `TransMat`'s `object$is.shape.ind` branch
  (`is.shape.ind` is never set by `CreateSpec`), and the `UnmapParameters_univ(x, "norm",
  FALSE)` path in `Mapping.cpp`, which returns an uninitialised vector — currently guarded by
  the `vDist[k] != "norm"` test at `R/StartingValues.R:189`, so unreachable.
- **`R CMD check` NOTE**: register `Sim` methods or align their signatures with the generic.
  The two WARNINGs are the version/date pair and a clang artefact, both cosmetic.

---

## Suggested order of work

1. ~~`R/Inference.R:38` (transpose) and `:43` (×2) — one line each, fixes every published
   standard error and p-value. **A1/A2.**~~ **DONE**
2. ~~`R/Utils.R:349-350` — two lines, makes saved fits reloadable. **B1.**~~ **DONE**
3. ~~`R/CondVol.R:32` — add the missing comma. **A3.**~~ **DONE**
4. ~~`R/Utils.R:420` — `* (K - 1)`. **A4.**~~ **DONE**
5. Guard the inverse in `src/MSgarch.h:283` and return `-1e10` instead of throwing; then fix
   the `FitML` guard at `R/FitML.R:133`. **B2/B5.**
6. `R/FitMCMC.R:148` + `f_check_parameterConstraints` — the two `fixed.pars` paths. **B3/B4.**
7. `src/MSgarch.h:497`, `:399`, `:474`; `R/State.R:40`; `R/UncVol.R:88`;
   `R/CreateSpec.R:213`. **C1–C5.**

Items 1–4 were seven lines of R and are done, with regression cover added for the two
most damaging defects:

- **`tests/testthat/test_Inference.R`** — pins the delta method against an independently
  computed natural-scale observed information (single-regime GARCH-Normal, where every
  parameter is interior and the numerical Hessian is well conditioned: the correct sandwich
  agrees to 1e-4 relative, the transposed one is off by 44% and 85%), asserts the shipped
  standard errors equal `sqrt(diag(J V J'))` to 1e-8 on both a single-regime and an MS(2)
  fit, and checks that `Pr(>|t|)` is two-sided. Each orientation test first asserts that the
  two orientations actually *differ* for that model, so it cannot silently go vacuous.
- **`tests/testthat/test_Serialization.R`** — round-trips a spec, an `MSGARCH_ML_FIT` and an
  `MSGARCH_MCMC_FIT` through `saveRDS`/`readRDS` (R restores external pointers as `NULL`, so
  this reproduces the cross-session failure within one session) and requires `Volatility`,
  `State`, `predict`, `AIC`, `summary`, `DIC` and the log-kernel to return values identical
  to the originals, with user-supplied priors preserved. It too opens with a guard asserting
  that the round trip really did invalidate the pointers.

- **`tests/testthat/test_Inference.R`**, two further blocks for **A4** — `dofMSGARCH` for
  K = 2, 3, 4 with a regime-constant `nu` (K = 2 cannot separate "one" from "K − 1", which is
  why the shipped code looked right), plus the unconstrained and `fixed.pars` cases as
  invariants, and an end-to-end check that `logLik()`'s `df` attribute and the `AIC`/`BIC`
  arithmetic follow it. The K = 3 and K = 4 legs run against a fit-shaped list rather than a
  real fit, because a constrained K ≥ 3 model cannot currently be fitted at all — that is
  **B2**, still open.
- **`tests/testthat/test_Volatility.R`**, one appended block for **A3** — recomputes the
  one-step-ahead volatility draw by draw through the public interface (each call carries a
  single parameter vector, so it cannot depend on how the draws are pooled) and requires
  `predict()` on the MCMC fit to equal their mean to 1e-10. Guarded by an assertion that the
  posterior mean is distinguishable from the first draw.

11 test blocks, 34 assertions, +13s of check time. All pass on the patched package and all 11
fail on stock 2.51 — block-level, measured, not inferred. Within the failing blocks some
individual assertions still pass on stock, by construction: the AIC/BIC block reports
`pass=3 fail=6` there because its K = 2 leg is exactly the case the two dof formulas agree on,
and the dof block reports `pass=3 fail=2` because K = 2 plus the two invariants are meant to
hold either way. The only *additions* that pass wholly on both are the deliberate
preconditions (the round trip really invalidates the pointers; the two sandwich orientations
really differ; the posterior mean really differs from draw #1), which exist so the tests
cannot quietly go vacuous.
