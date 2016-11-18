library(StudyPrior)

set.seed(3)
n <- sample(10:100, 5, replace = TRUE)
z <- runif(5)
x <- mapply(rbinom, size=n, n=1, prob=z)

Xs <- rbinom(1,100,0.7)
Ns <- 200

F.MFP <- binom.prior("MAP.FB", x = x, n=n)

F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=100, mc.cores=30)

## With feedback
F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=30)
F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:Ns, N=Ns, verbose=TRUE,  mc.cores=30)
F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:Ns, N=Ns)

# Approximations

C.MFP <- list(conj.approx2(F.MFP, "beta", max.degree=4, return.value=0.3, length.fit = 500))
C.MFP <- C.MFP[rep(1,201)]

C.MEP <- mclapply(mc.cores=30,0:200, function(X) {
  conj.approx2(function(p) F.MEP(p,X), "beta", max.degree=4, return.value=0.3, length.fit = 500)})

C.PFP <- list(conj.approx2(F.PFP, "beta", max.degree=4, return.value=0.3, length.fit = 500))
C.PFP <- C.PFP[rep(1,201)]

C.PEP <- mclapply(mc.cores=30,0:200, function(X) {
  conj.approx2(function(p) F.PEP(p,X), "beta", max.degree=4, return.value=0.3, length.fit = 500)})

C.PSP <- mclapply(mc.cores=30,0:200, function(X) {
  conj.approx2(function(p) F.PSP(p,X), "beta", max.degree=4, return.value=0.3, length.fit = 500)})



###############################################################

ESSMAT <-
  sapply(list(C.MFP,
              C.MEP,
              C.PFP,
              C.PEP,
              C.PSP),
         function(fl) sapply(0:Ns, function(Xs) ess.fun.list(posterior.fun.list(Xs,Ns, fl[[Xs+1]]))))

PROBMAT <- sapply(0:Ns, function(Xs) sapply(seq(0,1,len=100), function(p) dbinom(Xs, Ns, p)))

ess <- lapply(apply(PROBMAT %*% ESSMAT,2,list), `[[`, x=1)


mse.f <- lapply(list(
  F.MFP,
  F.MEP,
  F.PFP,
  F.PEP,
  F.PSP),
  function(pr) calc.MSE(prior=pr, prob.range=c(0,1), length = 100, mc.cores=30, n.binom=Ns))

mse <- lapply(list(
  C.MFP,
  C.MEP,
  C.PFP,
  C.PEP,
  C.PSP),
function(pr) calc.MSE(prior=pr, prob.range=c(0,1), length = 100, mc.cores=30, n.binom=Ns))

tf <- system.time(
bias.f <- mclapply(list(
  F.MFP,
  F.MEP,
  F.PFP,
  F.PEP,
  F.PSP),
  function(pr) calc.bias(prior=pr, prob.range=c(0,1), length = 100, n.binom=Ns), mc.cores=10)
)

tc <- system.time(
bias <- mclapply(list(
  C.MFP,
  C.MEP,
  C.PFP,
  C.PEP,
  C.PSP),
  function(pr) calc.bias(prior=pr, prob.range=c(0,1), length = 100, n.binom=Ns), mc.cores=10)
)

do.sigmat <- function(pr) sig.matrix(prior=pr, n.control=200, n.treatment = 200, check.xt=00:200, check.xs=0:200,
                        level=0.975, mc.cores=30, debug=FALSE)


  SIGMAT.f <- SIGMAT <- list()

  system.time(SIGMAT.f$F.MFP <- do.sigmat(F.MFP))
  system.time(SIGMAT.f$F.MEP <- do.sigmat(F.MEP))
  system.time(SIGMAT.f$F.PFP <- do.sigmat(F.PFP))
  system.time(SIGMAT.f$F.PEP <- do.sigmat(F.PEP))
  system.time(SIGMAT.f$F.PSP <- do.sigmat(F.PSP))

  system.time(SIGMAT$C.MFP <- do.sigmat(C.MFP))
  system.time(SIGMAT$C.MEP <- do.sigmat(C.MEP))
  system.time(SIGMAT$C.PFP <- do.sigmat(C.PFP))
  system.time(SIGMAT$C.PEP <- do.sigmat(C.PEP))
  system.time(SIGMAT$C.PSP <- do.sigmat(C.PSP))


pow.f <- mclapply(SIGMAT.f,
                function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.85), length =200, treatment.difference=0.12))

pow <- mclapply(SIGMAT,
                  function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.85), length =200, treatment.difference=0.12))

t1e <- mclapply(SIGMAT,
                function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.9), length = 200, treatment.difference = 0))


t1e.f <- mclapply(SIGMAT.f,
                function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.9), length = 200, treatment.difference = 0))


cover.f <-  mclapply(list(F.MFP,
                        F.MEP,
                        F.PFP,
                        F.PEP,
                        F.PSP),
                     function(pr) calc.coverage(prior=pr, level = 0.95, n.control = 200, smooth = 0.03), mc.cores=30)


cover <-  mclapply(list(C.MFP,
                        C.MEP,
                        C.PFP,
                        C.PEP,
                        C.PSP),
                   function(pr) calc.coverage(prior=pr, level = 0.95, n.control = 200, smooth = 0.03), mc.cores=30)

x.3 <- x
n.3 <- n

save(n.3, x.3, ess, bias,  cover, t1e, pow, mse,SIGMAT,  file="study_3.rda")
save(cover.f, bias.f, t1e.f, pow.f, mse.f,SIGMAT.f,  file="study_3_f.rda")
save(F.MFP, F.MEP,F.PFP,F.PEP,F.PSP,
     C.MFP, C.MEP,C.PFP,C.PEP,C.PSP,
     file="Study3_models.rda")
