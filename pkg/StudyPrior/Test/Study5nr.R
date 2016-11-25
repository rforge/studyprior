
library(StudyPrior)
recalc <- 1:100

mclapply(mc.cores=30, recalc, function(i){
# for(i in recalc){
  try({
    set.seed(300*i)
    N <- sample(6:10,1)
    n <- sample(50:100, N, replace = TRUE)
    z <- rnorm(N, 0.65, 0.1)
    x <- mapply(rbinom, size=n, n=1, prob=z)
    x/n
    Xs <- rbinom(1,100,0.6)
    Ns <- 200

    F.MFP <- binom.prior("MAP.FB", x = x, n=n)
    F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1)
    F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=100, mc.cores=1)
    ## With feedback
    F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:Ns, N=Ns, verbose=TRUE, mc.cores=1)
    F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:Ns, N=Ns)

    ## Approximations

    C.MFP <- list(conj.approx2(F.MFP, "beta", max.degree=4, return.value=0.01, length.fit = 100))

    C.MFP <- C.MFP[rep(1,201)]

    C.MEP <- vector("list", length=201)
    C.MEP[[101]] <- conj.approx2(function(p) F.MEP(p,100), "beta", max.degree=4, return.value=0.05, length.fit = 100)

    for(X in 101:200){
      starts <-  c(attr(C.MEP[[X]],"weights"),attr(C.MEP[[X]],"pars"))
      print(X)
      print(starts)
      C.MEP[[X+1]] <- conj.approx2(function(p) F.MEP(p,X), "beta", max.degree=4, return.value=0.05, length.fit = 100, starts=starts)
    }

    for(X in 100:0){
      starts <-  c(attr(C.MEP[[X+2]],"weights"),attr(C.MEP[[X+2]],"pars"))
      print(X)
      # print(starts)
      C.MEP[[X+1]] <- conj.approx2(function(p) F.MEP(p,X), "beta", max.degree=4, return.value=0.05, length.fit = 100, starts=starts)
    }

    C.PFP <- list(conj.approx2(F.PFP, "beta", max.degree=4, return.value=0.01, length.fit = 100))

    C.PFP <- C.PFP[rep(1,201)]


    C.PEP <- mclapply(mc.cores=1,0:200, function(X) {
      conj.approx2(function(p) F.PEP(p,X), "beta", max.degree=4, return.value=0.01, length.fit = 100)})

    C.PSP <- mclapply(mc.cores=1,0:200, function(X) {
      conj.approx2(function(p) F.PSP(p,X), "beta", max.degree=4, return.value=0.01, length.fit = 100)})


    save(F.MFP, F.MEP,F.PFP,F.PEP,F.PSP, n,x,
         C.MFP, C.MEP,C.PFP,C.PEP,C.PSP,
         file=paste0("5/d_study_5_models_",i,".rda"))

  })})

####################################################
recalc <- numeric()
for(i in 1:100){
  if(!file.exists(file=paste0("5/d_study_5_models_",i,".rda"))) recalc <- c(recalc,i)
}

for(i in 1:100){
  print(i)
  load(file=paste0("5/d_study_5_models_",i,".rda"))
  try({



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

    n.5 <- n
    x.5 <- x
    save(n.5, x.5, ess, bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("5/c_study_5_",i,".rda"))

  })
}
