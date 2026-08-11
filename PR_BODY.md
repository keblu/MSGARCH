Five bugs in the reporting and persistence layer around the likelihood, plus regression
tests for each. The likelihood itself is not touched: no estimate, log-likelihood,
conditional variance or state probability changes anywhere in this branch. What changes is
what `summary()` prints, what `AIC`/`BIC` count, what `predict()` returns for a Bayesian
fit, and whether a saved fit can be reloaded at all.

Everything below is reproduced on `data("SMI")` with the shipped code. `R CMD check
--as-cran` is unchanged by the branch (2 WARNINGs, 1 NOTE, `testthat` OK — and the version
bump in the last commit clears one of the two WARNINGs).

## 1. Standard errors used the delta method transposed — `R/Inference.R:38`

```r
mSandwitch <- t(mJacob) %*% mInvHessian %*% mJacob   # -> mJacob %*% mInvHessian %*% t(mJacob)
```

`numDeriv::jacobian` returns `∂f_i/∂x_j`, so `Var(g(θ̂)) = J V J'`. The two orientations
agree only if `J` is symmetric, and it is not: the working→natural map is triangular inside
each regime (the sGARCH bound on `beta` is `0.9999 − alpha1`; the gjrGARCH and tGARCH bounds
on `beta` also involve `alpha2` and the shape/skew parameters), and the transition-probability
block is *anti*-diagonal.

On the default MS(2)-GARCH(1,1)-Normal fit to `SMI`, checked against the observed information
computed directly in the natural parameterisation:

| | natural-scale `H` | fixed (`J V J'`) | before (`J' V J`) |
|---|---|---|---|
| `alpha1_1` | 0.01606 | 0.01511 | **0.03413** |
| `beta_1`   | 0.02197 | 0.02091 | **0.00958** |
| `alpha1_2` | 0.00431 | 0.00437 | **0.00610** |
| `beta_2`   | 0.00416 | 0.00426 | **0.00049** |
| `P_1_1`    | 0.00982 | 0.00973 | **0.00059** |
| `P_2_1`    | (at bound) | 0.03072 | **0.50250** |

Six of eight were wrong, by factors from 0.06× to 16×. Because the P block is anti-diagonal,
transposing **exchanged the two transition probabilities' standard errors**. Only `alpha0_k`,
whose map is a plain `exp`, was unaffected.

## 2. `Pr(>|t|)` was one-sided — `R/Inference.R:43`

`1 - pnorm(abs(t))` under a two-sided label; now `2 * (1 - pnorm(abs(t)))`.

## 3. `AIC`/`BIC` mis-counted `regime.const.pars` — `R/Utils.R:420`

A regime-constant parameter leaves one free value where there were `K`, so it removes `K − 1`
degrees of freedom (cf. `f_rename_par`, which strips `name_2 … name_K`). `dofMSGARCH`
subtracted one, which is right only at `K = 2`: with a regime-constant shape parameter the df
was 17 instead of 16 at `K = 3` and 27 instead of 25 at `K = 4`, always over-penalising the
constrained model.

## 4. Saved specs and fits could not be reloaded — `R/Utils.R:349-353`

R serializes external pointers as `NULL`, so a `saveRDS`-ed spec or fit comes back with dead
pointers. `f_check_spec` exists to rebuild them, but the rebuild branch called
`spec$rcpp.func$get_mean()` / `get_sd()` — the very pointer whose failure had just triggered
the branch:

```
Volatility(fit) : Error in .External(...): NULL value passed as symbol address
State(fit)      : Error in .External(...): NULL value passed as symbol address
predict(fit)    : Error in .External(...): NULL value passed as symbol address
```

The two values read there are already held on the R side in `spec$prior.mean` /
`spec$prior.sd`, which the next two lines were using anyway, so the C++ round trip was dead
code. Removing it makes the rebuild work; a reloaded spec, ML fit or MCMC fit now returns
values identical to before it was saved, with user priors preserved.

This is the ordinary workflow of fitting a model, saving it, and analysing it in a later
session — or shipping a spec to a `parLapply` worker.

## 5. `predict()` on an MCMC fit returned draw #1 — `R/CondVol.R:32`

`vol` is `(T+1) × ndraw` and the one-step-ahead value was `vol[dim(PredProb)[1]]` — a single
index into a matrix is linear indexing, i.e. the last row of the *first column*. On a
100-draw chain fitted to `SMI` the reported value was 1.019753 (the first draw) against a
posterior mean of 1.039776, with a spread of 0.999837–1.097337 across draws. `Volatility()`
already averaged correctly, so the two methods disagreed on the same fit. The
single-parameter (ML) path is unchanged.

## Tests

`test_Inference.R` and `test_Serialization.R` are new; `test_Volatility.R` gains one block.
11 blocks, 34 assertions, +13s of check time. All 11 fail on the current code.

Two things worth pointing out, since they are what makes the tests worth having:

- The standard-error test is anchored on a **single-regime** GARCH(1,1)-Normal, not on the
  MS(2) default. Every parameter there is interior, so a central-difference Hessian of the
  natural-scale negative log-likelihood is well conditioned: it agrees with `J V J'` to 1e-4
  in relative terms while the transposed sandwich is off by 44% and 85%. That is an
  independent check on the *value*, not a restatement of the formula.
- Every block opens with a guard asserting its own precondition — the two sandwich
  orientations really differ for this model; the round trip really did invalidate the
  pointers; the posterior mean really differs from the first draw — so none of them can pass
  vacuously if the surrounding code changes.

The `K = 3` and `K = 4` legs of the degrees-of-freedom test run against a fit-shaped list
rather than a real fit, because a constrained `K ≥ 3` model cannot currently be fitted at
all (see below).

## Not in this branch

The last commit (version bump + `NEWS`) is separable — drop it if the release number should
be decided elsewhere; nothing depends on it.

A review of the package turned up nine further issues that are **not** addressed here,
several of them more serious than some of the above. The two worth flagging now:

- **`src/MSgarch.h:283`** computes the ergodic distribution with a raw Armadillo `.i()` on
  `I − P + U`, on every likelihood evaluation. When `do.plm = TRUE` — forced by `fixed.pars`
  and `regime.const.pars`, and hard-coded in `FitMCMC` — the free transition entries are
  mapped into `(0,1)` independently, so for `K ≥ 3` a row can leave the simplex, the matrix
  can be exactly singular, and the uncaught exception aborts the entire run. Reproducible:
  `FitML` on a `K = 3` spec with `regime.const = "nu"` dies at
  `P = [[1, 0.999973, −0.999973], [0,1,0], [0,0,1]]`. This is why the tests above cannot fit
  a constrained `K ≥ 3` model.
- **`R/FitML.R:133`** guards optimisation failure with `if (llk == 1e+10)`, but `f_nll`
  returns `+1e10` so `llk` is `−1e10`; and `f_OptimFUNDefault` wraps `optim` in `try()`, so
  `optimizer$value` errors first. Every failure mode — including a single `NA` in the data,
  which `f_check_y` lets through — therefore surfaces as
  `Error in optimizer$value : $ operator is invalid for atomic vectors`.

Happy to open these as separate issues or as a follow-up PR, whichever you prefer.
