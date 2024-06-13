# This file demonstrates an implementation of the model in ortho_F.stan
# to simulate draws from the matrix von Mises-Fisher distribution
# with parameter Fmat.
#
# It is provided since there appears to be no such code available despite
# some discussion online and mention in some unpublished manuscripts.
#
# The draws are stored as the transformed parameter Q

setwd("~/set/your/directory")
library(rstan)

n = 20
p = 5
Fmat = matrix(rnorm(n*p),nrow=n)

ortho_data <- list(
  n = nrow(Fmat),
  r = ncol(Fmat),
  F = Fmat
)

t0=Sys.time()
fit <- stan(
  file = "ortho_F.stan",
  data = ortho_data,
  chains = 4,
  warmup = 20,
  iter = 1000,
  cores = 1,
  refresh = 200
)
t1=Sys.time()
t1-t0
print(fit, probs=c(.1,.5,.9))
# range of Rhat values
range(summary(fit)$summary[,10])

# extract the MvMF draws
draws_Q = extract(fit)$Q
dim(draws_Q) # 3920 x n x p
# one specific draw
draws_Q[1000,,]
