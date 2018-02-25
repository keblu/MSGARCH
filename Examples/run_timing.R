###############################################################
## TIMIMNG COMPARISON
## - ML estimation vs. MCMC estimation
## - MSGARCH package vs. bayesGARCH package for GARCH(1,1)-N

# > sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=French_Switzerland.1252  LC_CTYPE=French_Switzerland.1252    LC_MONETARY=French_Switzerland.1252
# [4] LC_NUMERIC=C                        LC_TIME=French_Switzerland.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] microbenchmark_1.4-4 zoo_1.8-1            compiler_3.4.3       Matrix_1.2-12        MSGARCH_2.0         
# [6] tools_3.4.3          fanplot_3.4.1        expm_0.999-2         coda_0.19-1          Rcpp_0.12.15        
# [11] mvtnorm_1.0-7        ellipse_0.4.1        codetools_0.2-15     grid_3.4.3           bayesGARCH_2.1.3    
# [16] lattice_0.20-35      mcmcse_1.3-2        

###############################################################
## INSTALL PACKAGES
# install.packages("MSGARCH")
# install.packages("bayesGARCH")
# install.packages("microbenchmark")
# install.packages("mcmcse")

###############################################################
## ML vs MCMC timing comparison for MS2-GARCH(1,1)-N model

rm(list = ls())
library("MSGARCH")
library("microbenchmark")

data("SMI", package = "MSGARCH")
y <- as.vector(SMI)

spec <- MSGARCH::CreateSpec()

f.ML <- function() {
  MSGARCH::FitML(spec = spec, data = y)
}

f.MCMC <- function() {
  MSGARCH::FitMCMC(spec = spec, data = y)
}

set.seed(1234)
microbenchmark::microbenchmark(f.ML(), f.MCMC(), times = 4L)
# Unit: seconds
# expr     min        lq      mean    median        uq       max neval
# f.ML()    2.5490  2.638808  2.985428  2.904798  3.332047  3.583115     4
# f.MCMC() 24.7685 25.580564 26.924724 26.427310 28.268884 30.075782     4


###############################################################
## MSGARCH vs bayesGARCH timing comparison for GARCH(1,1)-N

rm(list = ls())
library("MSGARCH")
library("microbenchmark")
library("bayesGARCH")
library("mcmcse")

data("SMI", package = "MSGARCH")
y <- as.vector(SMI)

# same MCMC setup
n.bi    <- 5000
n.chain <- 10000
n.thin  <- 10

f.MSGARCH <- function() {
  spec <- MSGARCH::CreateSpec(variance.spec = list(model = "sGARCH"),
                              distribution.spec = list(distribution = "norm"),
                              switch.spec = list(K = 1))
  
  fit <- MSGARCH::FitMCMC(spec = spec, data = y, 
                          ctr = list(nburn = n.bi, nmcmc = n.chain, nthin = n.thin))
  
  out <- mcmcse::ess(fit$par)
  return(out)
}

f.bayesGARCH <- function() {
  addPriorConditions <- function(psi) {
    psi[2] + psi[3] < 1
  }
  
  fit <- bayesGARCH::bayesGARCH(y = y, lambda = 100, delta = 200, 
                                control = list(n.chain = 1, l.chain = n.bi + n.chain, 
                                               addPriorConditions = addPriorConditions, 
                                               refresh = 0))
  invisible(capture.output(smpl <- bayesGARCH::formSmpl(fit, l.bi = n.bi,
                                         batch.size = n.thin)))
  out <- mcmcse::ess(smpl)
  return(out)
}


set.seed(1234)
microbenchmark::microbenchmark(f.MSGARCH(), f.bayesGARCH(), times = 4L)
# Unit: seconds
# expr       min        lq      mean    median        uq       max neval
# f.MSGARCH()     5.737974  5.776521  5.860372  5.864642  5.944222  5.974228     4
# f.bayesGARCH() 62.773728 63.024597 64.590266 64.213374 66.155935 67.160588     4

# Estimate effective sample size (ESS) as described in Gong and Flegal (2015)
# ESS is the size of an iid sample with the same variance as the current sample

ess.MSGARCH    <- f.MSGARCH()
ess.bayesGARCH <- f.bayesGARCH()

ess.MSGARCH
# alpha0_1 alpha1_1   beta_1 
# 775.5998 943.7135 703.0423 

ess.bayesGARCH
# alpha0   alpha1     beta       nu
# 250.5667 275.4786 204.2709 979.6102 


