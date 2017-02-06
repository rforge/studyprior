
library(StudyPrior)

set.seed(3)
n <- sample(10:100, 10, replace = TRUE)
z <- 0.6 + runif(10)/5
x <- mapply(rbinom, size=n, n=1, prob=z)

Xs <- rbinom(1,100,0.6)
Ns <- 200

F.MFP <- binom.prior("MAP.FB", x = x, n=n)
F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=30)
F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=100, mc.cores=30)
F.PJP <- binom.PP.FB.MC.BE(x,n,d.prior.a = .5,d.prior.b = 5,samples=5000, length=100, mc.cores=30)

## With feedback
F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:Ns, N=Ns, verbose=TRUE, mc.cores=30)
F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:Ns, N=Ns)

## Approximations

C.MFP <- list(conj.approx2(F.MFP, "beta", max.degree=4, return.value=0.1, length.fit = 300))

C.MFP <- C.MFP[rep(1,201)]

C.MEP <- vector("list", length=201)
C.MEP[[101]] <- conj.approx2(function(p) F.MEP(p,100), "beta",min.degree = 3, max.degree=3, return.value=0.05, length.fit = 200)


for(X in 101:200){
  starts <-  c(attr(C.MEP[[X]],"weights")[1:3],attr(C.MEP[[X]],"pars")[1:6])
  print(X)
  # print(starts)
  C.MEP[[X+1]] <- conj.approx2(function(p) F.MEP(p,X), "beta",min.degree = 3, max.degree=3, return.value=0.05, length.fit = 200, starts=starts)
}

for(X in 100:0){
  starts <-  c(attr(C.MEP[[X+2]],"weights")[1:3],attr(C.MEP[[X+2]],"pars")[1:6])
  print(X)
  # print(starts)
  C.MEP[[X+1]] <- conj.approx2(function(p) F.MEP(p,X), "beta",min.degree = 3, max.degree=3, return.value=0.05, length.fit = 200, starts=starts)
}

C.PFP <- list(conj.approx2(F.PFP, "beta", max.degree=4, return.value=0.1, length.fit = 300))

C.PFP <- C.PFP[rep(1,201)]


C.PEP <- mclapply(mc.cores=30,0:200, function(X) {
  conj.approx2(function(p) F.PEP(p,X), "beta", max.degree=4, return.value=0.1, length.fit = 300)})

C.PSP <- mclapply(mc.cores=30,0:200, function(X) {
  conj.approx2(function(p) F.PSP(p,X), "beta", max.degree=5, return.value=0.1, length.fit = 300)})

####################################################

ESSMAT <-
  sapply(list(C.MFP,
              C.MEP,
              C.PFP,
              C.PEP,
              C.PSP),
         function(fl) sapply(0:Ns, function(Xs) ess.fun.list(posterior.fun.list(Xs,Ns, fl[[Xs+1]]))))

PROBMAT <- sapply(0:Ns, function(Xs) sapply(seq(0,1,len=100), function(p) dbinom(Xs, Ns, p)))

ess <- lapply(apply(PROBMAT %*% ESSMAT,2,list), `[[`, x=1)



mse <- lapply(list(
  C.MFP,
  C.MEP,
  C.PFP,
  C.PEP,
  C.PSP),
  function(pr) calc.MSE(prior=pr, prob.range=c(0,1), length = 100, mc.cores=30, n.binom=Ns))


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


SIGMAT <- list()

system.time(SIGMAT$C.MFP <- do.sigmat(C.MFP))
system.time(SIGMAT$C.MEP <- do.sigmat(C.MEP))
system.time(SIGMAT$C.PFP <- do.sigmat(C.PFP))
system.time(SIGMAT$C.PEP <- do.sigmat(C.PEP))
system.time(SIGMAT$C.PSP <- do.sigmat(C.PSP))



pow <- mclapply(SIGMAT,
                function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.85), length =200, treatment.difference=0.12))

t1e <- mclapply(SIGMAT,
                function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.9), length = 200, treatment.difference = 0))


cover <-  mclapply(list(C.MFP,
                        C.MEP,
                        C.PFP,
                        C.PEP,
                        C.PSP),
                   function(pr) calc.coverage(prior=pr, level = 0.95, n.control = 200, smooth = 0.03), mc.cores=30)

n.4 <- n
x.4 <- x
save(n.4, x.4, ess, bias,  cover, t1e, pow, mse,SIGMAT,  file="study_4.rda")
save(F.MFP, F.MEP,F.PFP,F.PEP,F.PSP,
     C.MFP, C.MEP,C.PFP,C.PEP,C.PSP,
     file="Study4_models.rda")
