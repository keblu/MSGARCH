# Independent Codex cross-check of REVIEW.md

Run 2026-08-10, `codex-cli 0.139.0`, read-only over `package MSGARCH/`:

```bash
codex exec --sandbox read-only --skip-git-repo-check "<audit prompt>"
```

Codex audited (1) the 14 load-bearing claims, (2) the five applied fixes, and (3) the new
tests, by reading the source only — it could not run R. Verbatim output below; see the
"Independent cross-check" section of `REVIEW.md` for how the three pushbacks were resolved.

---

## Part 1 — claim verdicts

| # | claim | verdict | evidence (file:line) |
|---|---|---|---|
| A1 | delta-method sandwich transposed; P-block anti-diagonal | CONFIRM | Original code uses `t(mJacob) %*% mInvHessian %*% mJacob` while `jacobian(f_mapPar, ...)` returns output-by-input derivatives ([`Inference.R:36-38`](Package/R/Inference.R:36)). The P block is anti-diagonal for `K=2`: `f_mapGamma()` fills off-diagonals by `mGamma[!mGamma]` then returns `c(t(mGamma[, -K]))`, keeping the original names ([`ParameterTransformation.R:204-209`](Package/R/ParameterTransformation.R:204)). |
| A2 | p-value is one-sided under `Pr(>|t|)` | CONFIRM | The column is named `Pr(>|t|)` at construction, but the value is `1 - pnorm(abs(vTest))` ([`Inference.R:4-5`, `:43`](Package/R/Inference.R:4)). |
| A3 | MCMC forecast uses linear index draw #1 | CONFIRM | `vol` is allocated as `(T+1) x nrow(par.check)` and line 32 indexes it as `vol[dim(PredProb)[1]]`, a single linear index ([`CondVol.R:19`, `:32`](Package/R/CondVol.R:19)). |
| A4 | `dofMSGARCH` subtracts one instead of `K-1` per regime-constant parameter | CONFIRM | `dofMSGARCH` subtracts `length(regime.const.pars)` only, while `f_rename_par()` removes `name_2 ... name_K` for each constrained base name ([`Utils.R:420-421`, `:364-367`](Package/R/Utils.R:420)). |
| B1 | stale-pointer recovery re-dereferences dead pointer; lines 349-353 unused | CONFIRM | The failed branch calls `spec$rcpp.func$get_mean()` and `get_sd()` immediately after `get_sd()` failed; computed `prior.mean/prior.sd` are not used by lines 355-356 ([`Utils.R:340-356`](Package/R/Utils.R:340)). |
| B2 | raw `.i()` can abort under invalid K>=3 PLM transition rows; K=2 nonsingular | CONFIRM | `loadparam()` forms rows via `extract_P_it()` with last entry `1 - sum(P_it)`, then calls `(foo).i()` without a guard ([`MSgarch.h:193-198`, `:281-283`](Package/src/MSgarch.h:193)). For `K=2`, `det(I-P+U)=2(1+P_2_1-P_1_1)>0` when both free entries are in `(0,1)`. |
| B3 | FitMCMC custom `par0` fixed-parameter path too long/mis-mapped | PARTIALLY CONFIRM | It is true FitMCMC substitutes fixed values but does not remove them before `f_unmapPar()` ([`FitMCMC.R:143-150`](Package/R/FitMCMC.R:143)), unlike FitML’s `f_remove_fixedpar()` ([`FitML.R:125-126`](Package/R/FitML.R:125)). But the exact “trailing NA name” mechanism is overstated: with names present the fixed dimension remains in the sampler, and with names absent `f_rename_par()` is more likely to fail length/name assignment than silently trail `NA` ([`Utils.R:370-374`](Package/R/Utils.R:370)). |
| B4 | fixed transition probability accepted but breaks fitters by two mechanisms | CONFIRM | `CreateSpec()` validates against `out$label`, which includes transition labels ([`CreateSpec.R:278`, `:397-400`](Package/R/CreateSpec.R:397)); starting values pass recovered fixed names into single-regime `CreateSpec()` ([`StartingValues.R:201-206`](Package/R/StartingValues.R:201)). Separately, prior correction indexes `prior.mean[names(fixed.pars)]`, but MS priors come only from within-regime model params ([`Kernel.R:50-57`; `MSgarch.h:146-171`](Package/R/Kernel.R:50)). |
| B5 | `llk == 1e+10` check wrong and try-error not handled | CONFIRM | `f_nll()` returns `+1e10` on likelihood failure, so `llk <- -optimizer$value` is `-1e10`; if `optim` itself fails, `f_OptimFUNDefault()` returns a `"try-error"` and `optimizer$value` is invalid ([`Likelihood.R:20-24`; `FitML.R:129-136`; `Optimizers.R:2-4`](Package/R/FitML.R:129)). |
| C1 | MS `f_pdf`/`f_cdf` log path logs last-regime tmp, not mixture | CONFIRM | Both functions accumulate into `out`, then under `is_log` assign `out[i] = log(tmp[i])`, where `tmp` is overwritten in the regime loop ([`MSgarch.h:387-400`, `:462-475`](Package/src/MSgarch.h:387)). |
| C2 | `f_cdf_its` first slice has swapped cube indices | CONFIRM | The initial CDF slice writes `tmp(ix, 0, s)` while later rows and `f_pdf_its` write `tmp(i, ix, s)` / `tmp(0, ix, s)` ([`MSgarch.h:490-507`; compare `:420-435`](Package/src/MSgarch.h:490)). |
| C3 | `UncVol` window wrong and defaults swapped vs docs | CONFIRM | Docs say `nahead=5000`, `nburn=1000`, but `f_process_ctr(type=2)` sets `nburn=5000`, `nahead=1000`; `UncVol` then averages `tmp[ctr$nburn:ctr$nahead]` ([`UncVol.R:16-21`, `:88`; `Utils.R:159-160`](Package/R/UncVol.R:16)). |
| C4 | mixture Viterbi transition matrix transposed/row sums wrong | CONFIRM | For mixtures, `TransMat()` returns a one-row probability vector; `State()` expands with `matrix(rep(P, K), nrow=K, ncol=K)`, which fills columns, so rows become repeated single probabilities rather than the mixture vector ([`State.R:90-94`; `Transmat.R:53-58`](Package/R/State.R:90)). |
| C5 | K-expansion guard checks nonexistent `distribution.spec$model` | CONFIRM | The guard tests `length(distribution.spec$model) > 1`, but the checked object uses `$distribution` ([`CreateSpec.R:205`, `:212-219`](Package/R/CreateSpec.R:205)). |

On the negative claims: I found no source-level contradiction to REVIEW.md’s claim that the base conditional densities/CDFs/quantiles, Fernandez-Steel moments, Hamilton filter, Kim smoother, valid-matrix ergodic formula, transition-memory layout, and tGARCH/gjrGARCH stationarity formulas are correct. Caveat: the ergodic formula is mathematically fine for valid stochastic matrices, but B2 shows the implementation does not guard invalid PLM matrices before inverting.

## Part 2 — fixes

A1: `mJacob %*% mInvHessian %*% t(mJacob)` is the correct orientation for `numDeriv::jacobian(f_mapPar, ...)`. The boundary adjustment changes the evaluation point to a nearby interior `vPw_mod`; it does not alter the orientation argument, and the Hessian is also computed at `vPw_mod` when recomputed. OK

A2: `2 * (1 - pnorm(abs(vTest)))` matches the printed `Pr(>|t|)` label. The broader statistical caveat remains that zero-null z/t tests are not meaningful for all bounded/support-shifted parameters, but this fix corrects the advertised two-sided calculation. OK

A3: `mean(vol[dim(PredProb)[1], ])` preserves ML behavior because a one-draw matrix row averages to the same scalar. Averaging posterior volatility rather than variance is a modeling/reporting choice already implied by `Volatility()`’s row-mean-on-volatility path ([`CondVol.R:53-55`](Package/R/CondVol.R:53)); it is separate from the indexing bug. OK

A4: `- length(regime.const.pars) * (K - 1L)` is right for allowed specs: `K=1` rejects non-null regime constants, and `CreateSpec()` forbids using `fixed.pars` and `regime.const.pars` together ([`ParameterConstraints.R:115-121`; `CreateSpec.R:293-294`](Package/R/CreateSpec.R:293)). It also works for mixtures because the constrained duplicated within-regime parameters are still K copies collapsed to one. OK

B1: Rebuilding only `rcpp.func` is sufficient from this source: all external-pointer-backed callables are collected there by `f_spec()` ([`CreateSpec.R:350-366`](Package/R/CreateSpec.R:350)), while `func` contains plain R closures. The `!is.null` guards should not affect normal `f_spec()` / `ExtractStateFit()` specs because `prior.mean` and `prior.sd` are initialized non-null ([`CreateSpec.R:412-416`](Package/R/CreateSpec.R:412)). OK

## Part 3 — tests

`test_Inference.R` “Standard errors use J V J'”: fails on original, has a non-vacuous guard, but partly restates implementation via `fit$Inference$Hessian`. OK

`test_Inference.R` “observed information in natural scale”: independent and would fail on original; finite-difference tolerance and full `SMI` ML fit make it somewhat machine/runtime sensitive. WEAK

`test_Inference.R` “Transition-probability standard errors”: fails on original and directly covers the P-block; uses fixed historical parameters, so it is reasonably stable. OK

`test_Inference.R` “Pr(>|t|)”: fails on original and is not vacuous. OK

`test_Inference.R` “Degrees of freedom”: catches original for `K=3/4`; independent enough. OK

`test_Inference.R` “AIC and BIC”: the real fitted case is only `K=2`, where old and new formulas coincide; affected `K=3/4` cases are synthetic objects that test plumbing, not fitting. WEAK

`test_Serialization.R` “round trip invalidates pointers”: useful anti-vacuity guard, but depends on current Rcpp serialization behavior. OK

`test_Serialization.R` “reloaded spec usable”: fails on original and checks outputs against a live spec. OK

`test_Serialization.R` “keeps priors”: fails on original and checks rebuilt C++ prior state. OK

`test_Serialization.R` “reloaded ML fit usable”: fails on original; runtime acceptable but includes a full ML fit. OK

`test_Serialization.R` “reloaded MCMC fit usable”: fails on original but is RNG/runtime-sensitive and uses a very short chain. WEAK

`test_Volatility.R` “MCMC forecast averages draws”: should fail on original and has a non-vacuity guard; fragile because it runs MCMC in unit tests and depends on the first draw differing from the mean. WEAK

## Missed or overstated

B3 is overstated. The core bug is real: FitMCMC leaves the fixed parameter in user-supplied `ctr$par0` ([`FitMCMC.R:143-150`](Package/R/FitMCMC.R:143)). But the claimed “trailing NA name” is not the clean source-level failure mode: if names survive, the sampler has an extra non-fixed dimension; if names are stripped, `f_rename_par()` is asked to assign a shorter name vector ([`Utils.R:370-374`](Package/R/Utils.R:370)).

The AIC/BIC new test overclaims end-to-end coverage for the actual affected case. Its only real fit is `K=2` ([`test_Inference.R:138-143`](Package/tests/testthat/test_Inference.R:138)), where subtracting `1` and subtracting `K-1` are identical.

REVIEW.md says all 11 new blocks fail on stock 2.51. Source-wise that is plausible for the five bug areas, but not literally demonstrated by every block: the AIC/BIC fitted subcase at `K=2` would pass under the old formula, and only the synthetic `K=3/4` part fails.
tokens used
97,458
## Part 1 — claim verdicts

| # | claim | verdict | evidence (file:line) |
|---|---|---|---|
| A1 | delta-method sandwich transposed; P-block anti-diagonal | CONFIRM | Original code uses `t(mJacob) %*% mInvHessian %*% mJacob` while `jacobian(f_mapPar, ...)` returns output-by-input derivatives ([`Inference.R:36-38`](Package/R/Inference.R:36)). The P block is anti-diagonal for `K=2`: `f_mapGamma()` fills off-diagonals by `mGamma[!mGamma]` then returns `c(t(mGamma[, -K]))`, keeping the original names ([`ParameterTransformation.R:204-209`](Package/R/ParameterTransformation.R:204)). |
| A2 | p-value is one-sided under `Pr(>|t|)` | CONFIRM | The column is named `Pr(>|t|)` at construction, but the value is `1 - pnorm(abs(vTest))` ([`Inference.R:4-5`, `:43`](Package/R/Inference.R:4)). |
| A3 | MCMC forecast uses linear index draw #1 | CONFIRM | `vol` is allocated as `(T+1) x nrow(par.check)` and line 32 indexes it as `vol[dim(PredProb)[1]]`, a single linear index ([`CondVol.R:19`, `:32`](Package/R/CondVol.R:19)). |
| A4 | `dofMSGARCH` subtracts one instead of `K-1` per regime-constant parameter | CONFIRM | `dofMSGARCH` subtracts `length(regime.const.pars)` only, while `f_rename_par()` removes `name_2 ... name_K` for each constrained base name ([`Utils.R:420-421`, `:364-367`](Package/R/Utils.R:420)). |
| B1 | stale-pointer recovery re-dereferences dead pointer; lines 349-353 unused | CONFIRM | The failed branch calls `spec$rcpp.func$get_mean()` and `get_sd()` immediately after `get_sd()` failed; computed `prior.mean/prior.sd` are not used by lines 355-356 ([`Utils.R:340-356`](Package/R/Utils.R:340)). |
| B2 | raw `.i()` can abort under invalid K>=3 PLM transition rows; K=2 nonsingular | CONFIRM | `loadparam()` forms rows via `extract_P_it()` with last entry `1 - sum(P_it)`, then calls `(foo).i()` without a guard ([`MSgarch.h:193-198`, `:281-283`](Package/src/MSgarch.h:193)). For `K=2`, `det(I-P+U)=2(1+P_2_1-P_1_1)>0` when both free entries are in `(0,1)`. |
| B3 | FitMCMC custom `par0` fixed-parameter path too long/mis-mapped | PARTIALLY CONFIRM | It is true FitMCMC substitutes fixed values but does not remove them before `f_unmapPar()` ([`FitMCMC.R:143-150`](Package/R/FitMCMC.R:143)), unlike FitML’s `f_remove_fixedpar()` ([`FitML.R:125-126`](Package/R/FitML.R:125)). But the exact “trailing NA name” mechanism is overstated: with names present the fixed dimension remains in the sampler, and with names absent `f_rename_par()` is more likely to fail length/name assignment than silently trail `NA` ([`Utils.R:370-374`](Package/R/Utils.R:370)). |
| B4 | fixed transition probability accepted but breaks fitters by two mechanisms | CONFIRM | `CreateSpec()` validates against `out$label`, which includes transition labels ([`CreateSpec.R:278`, `:397-400`](Package/R/CreateSpec.R:397)); starting values pass recovered fixed names into single-regime `CreateSpec()` ([`StartingValues.R:201-206`](Package/R/StartingValues.R:201)). Separately, prior correction indexes `prior.mean[names(fixed.pars)]`, but MS priors come only from within-regime model params ([`Kernel.R:50-57`; `MSgarch.h:146-171`](Package/R/Kernel.R:50)). |
| B5 | `llk == 1e+10` check wrong and try-error not handled | CONFIRM | `f_nll()` returns `+1e10` on likelihood failure, so `llk <- -optimizer$value` is `-1e10`; if `optim` itself fails, `f_OptimFUNDefault()` returns a `"try-error"` and `optimizer$value` is invalid ([`Likelihood.R:20-24`; `FitML.R:129-136`; `Optimizers.R:2-4`](Package/R/FitML.R:129)). |
| C1 | MS `f_pdf`/`f_cdf` log path logs last-regime tmp, not mixture | CONFIRM | Both functions accumulate into `out`, then under `is_log` assign `out[i] = log(tmp[i])`, where `tmp` is overwritten in the regime loop ([`MSgarch.h:387-400`, `:462-475`](Package/src/MSgarch.h:387)). |
| C2 | `f_cdf_its` first slice has swapped cube indices | CONFIRM | The initial CDF slice writes `tmp(ix, 0, s)` while later rows and `f_pdf_its` write `tmp(i, ix, s)` / `tmp(0, ix, s)` ([`MSgarch.h:490-507`; compare `:420-435`](Package/src/MSgarch.h:490)). |
| C3 | `UncVol` window wrong and defaults swapped vs docs | CONFIRM | Docs say `nahead=5000`, `nburn=1000`, but `f_process_ctr(type=2)` sets `nburn=5000`, `nahead=1000`; `UncVol` then averages `tmp[ctr$nburn:ctr$nahead]` ([`UncVol.R:16-21`, `:88`; `Utils.R:159-160`](Package/R/UncVol.R:16)). |
| C4 | mixture Viterbi transition matrix transposed/row sums wrong | CONFIRM | For mixtures, `TransMat()` returns a one-row probability vector; `State()` expands with `matrix(rep(P, K), nrow=K, ncol=K)`, which fills columns, so rows become repeated single probabilities rather than the mixture vector ([`State.R:90-94`; `Transmat.R:53-58`](Package/R/State.R:90)). |
| C5 | K-expansion guard checks nonexistent `distribution.spec$model` | CONFIRM | The guard tests `length(distribution.spec$model) > 1`, but the checked object uses `$distribution` ([`CreateSpec.R:205`, `:212-219`](Package/R/CreateSpec.R:205)). |

On the negative claims: I found no source-level contradiction to REVIEW.md’s claim that the base conditional densities/CDFs/quantiles, Fernandez-Steel moments, Hamilton filter, Kim smoother, valid-matrix ergodic formula, transition-memory layout, and tGARCH/gjrGARCH stationarity formulas are correct. Caveat: the ergodic formula is mathematically fine for valid stochastic matrices, but B2 shows the implementation does not guard invalid PLM matrices before inverting.

## Part 2 — fixes

A1: `mJacob %*% mInvHessian %*% t(mJacob)` is the correct orientation for `numDeriv::jacobian(f_mapPar, ...)`. The boundary adjustment changes the evaluation point to a nearby interior `vPw_mod`; it does not alter the orientation argument, and the Hessian is also computed at `vPw_mod` when recomputed. OK

A2: `2 * (1 - pnorm(abs(vTest)))` matches the printed `Pr(>|t|)` label. The broader statistical caveat remains that zero-null z/t tests are not meaningful for all bounded/support-shifted parameters, but this fix corrects the advertised two-sided calculation. OK

A3: `mean(vol[dim(PredProb)[1], ])` preserves ML behavior because a one-draw matrix row averages to the same scalar. Averaging posterior volatility rather than variance is a modeling/reporting choice already implied by `Volatility()`’s row-mean-on-volatility path ([`CondVol.R:53-55`](Package/R/CondVol.R:53)); it is separate from the indexing bug. OK

A4: `- length(regime.const.pars) * (K - 1L)` is right for allowed specs: `K=1` rejects non-null regime constants, and `CreateSpec()` forbids using `fixed.pars` and `regime.const.pars` together ([`ParameterConstraints.R:115-121`; `CreateSpec.R:293-294`](Package/R/CreateSpec.R:293)). It also works for mixtures because the constrained duplicated within-regime parameters are still K copies collapsed to one. OK

B1: Rebuilding only `rcpp.func` is sufficient from this source: all external-pointer-backed callables are collected there by `f_spec()` ([`CreateSpec.R:350-366`](Package/R/CreateSpec.R:350)), while `func` contains plain R closures. The `!is.null` guards should not affect normal `f_spec()` / `ExtractStateFit()` specs because `prior.mean` and `prior.sd` are initialized non-null ([`CreateSpec.R:412-416`](Package/R/CreateSpec.R:412)). OK

## Part 3 — tests

`test_Inference.R` “Standard errors use J V J'”: fails on original, has a non-vacuous guard, but partly restates implementation via `fit$Inference$Hessian`. OK

`test_Inference.R` “observed information in natural scale”: independent and would fail on original; finite-difference tolerance and full `SMI` ML fit make it somewhat machine/runtime sensitive. WEAK

`test_Inference.R` “Transition-probability standard errors”: fails on original and directly covers the P-block; uses fixed historical parameters, so it is reasonably stable. OK

`test_Inference.R` “Pr(>|t|)”: fails on original and is not vacuous. OK

`test_Inference.R` “Degrees of freedom”: catches original for `K=3/4`; independent enough. OK

`test_Inference.R` “AIC and BIC”: the real fitted case is only `K=2`, where old and new formulas coincide; affected `K=3/4` cases are synthetic objects that test plumbing, not fitting. WEAK

`test_Serialization.R` “round trip invalidates pointers”: useful anti-vacuity guard, but depends on current Rcpp serialization behavior. OK

`test_Serialization.R` “reloaded spec usable”: fails on original and checks outputs against a live spec. OK

`test_Serialization.R` “keeps priors”: fails on original and checks rebuilt C++ prior state. OK

`test_Serialization.R` “reloaded ML fit usable”: fails on original; runtime acceptable but includes a full ML fit. OK

`test_Serialization.R` “reloaded MCMC fit usable”: fails on original but is RNG/runtime-sensitive and uses a very short chain. WEAK

`test_Volatility.R` “MCMC forecast averages draws”: should fail on original and has a non-vacuity guard; fragile because it runs MCMC in unit tests and depends on the first draw differing from the mean. WEAK

## Missed or overstated

B3 is overstated. The core bug is real: FitMCMC leaves the fixed parameter in user-supplied `ctr$par0` ([`FitMCMC.R:143-150`](Package/R/FitMCMC.R:143)). But the claimed “trailing NA name” is not the clean source-level failure mode: if names survive, the sampler has an extra non-fixed dimension; if names are stripped, `f_rename_par()` is asked to assign a shorter name vector ([`Utils.R:370-374`](Package/R/Utils.R:370)).

The AIC/BIC new test overclaims end-to-end coverage for the actual affected case. Its only real fit is `K=2` ([`test_Inference.R:138-143`](Package/tests/testthat/test_Inference.R:138)), where subtracting `1` and subtracting `K-1` are identical.

REVIEW.md says all 11 new blocks fail on stock 2.51. Source-wise that is plausible for the five bug areas, but not literally demonstrated by every block: the AIC/BIC fitted subcase at `K=2` would pass under the old formula, and only the synthetic `K=3/4` part fails.
