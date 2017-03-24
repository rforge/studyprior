
library(foreach)
library(StudyPrior)
inla.setOption(num.threads=2)
recalc <- 1:10

mclapply(mc.cores=20, recalc, function(i){
  # lapply( recalc, function(i){
  # for(i in recalc){
  print(i)

  set.seed(300*i)
  N <- 5
  n <- rep(50,N)
  z <- 0.65
  x <- mapply(rbinom, size=n, n=1, prob=z)
  x/n
  Xs <- rbinom(1,100,0.6)
  Ns <- 200
  cat(i,".\n")
  F.MFP <- binom.prior("MAP.FB", x = x, n=n)
  cat(i,"..\n")
  F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1, verbose=FALSE)
  cat(i,"...\n")
  F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=512, mc.cores=1, verbose=FALSE)
  ## With feedback
  cat(i,"....\n")
  F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:Ns, N=Ns, verbose=FALSE, mc.cores=1)
  cat(i,".....\n")
  F.FX0 <- binom.prior("PP.Fix", x = x, n=n, d=0)
  cat(i,"......\n")
  F.COR <- binom.prior("PP.Cor", x = x, n=n, d.prior.cor=0.5, samples=5000, length=512)
  cat(i,".......\n")
  F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:Ns, N=Ns)
  cat(i,"........\n")
  F.SGL <- binom.prior("PP.EB", x = sum(x), n=sum(n), X=0:Ns, N=Ns, verbose=FALSE, mc.cores=1)
  cat(i,".........\n")
  save(F.MFP, F.MEP,F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL, n,x,
       file=paste0("Fix-5/models_f_",i,".rda"))
}
)


mclapply(mc.cores=25, recalc, function(i){
  load(file=paste0("Fix-5/models_f_",i,".rda"))
  print(i)

  set.seed(300*i)
  N <- 5
  n <- rep(50,N)
  z <- 0.65
  x <- mapply(rbinom, size=n, n=1, prob=z)
  x/n
  Xs <- rbinom(1,100,0.6)
  Ns <- 200
  cat(i,".\n")
  F.MFP <- binom.prior("MAP.FB", x = x, n=n)
  cat(i,"..\n")
  F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1, verbose=FALSE)
  F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=512, mc.cores=1, verbose=FALSE)
  save(F.MFP, F.MEP,F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL, n,x,
       file=paste0("Fix-5/models_f_",i,".rda"))
}
)


Calc.posterior <- function(prior, X, N){
  pr <- function(p) prior(p,X)
  k <- adaptIntegrate(function(p) pr(p)*dbinom(X,N,p),
                      lowerLimit = 0,
                      upperLimit = 1)$integral
  function(p) pr(p)*dbinom(X,N,p)/k
}

Calc.posterior.all <- function(prior, N, mc.cores){
  mclapply(0:N, function(X) Calc.posterior(prior, X, N), mc.cores = mc.cores)
}



library(foreach)
library(doParallel)

####################################################

for(i in 1:10){
  print(i)
  try({
    load(file=paste0("Fix-5/models_f_",i,".rda"))

    Ns <- 200

    posteriors <-
    lapply(list(
      F.MFP, F.MEP,F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL)[2],
      Calc.posterior.all, N=Ns, mc.cores=34
    )


    mse <- lapply(posteriors,
      function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=30, n.binom=Ns))


    tc <- system.time(
      bias <- lapply(posteriors,
        function(pr) {
         print("Starting!")
           calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=30)
        }
    ))

    do.sigmat <- function(pr) sig.matrix(posterior = pr, n.control=200, n.treatment = 200,
                                         check.xt=00:200, check.xs=0:200,
                                         level=0.975, mc.cores=34, debug=FALSE)


    SIGMAT <- lapply(posteriors, function(p){
      print("SIGMATING!")
        do.sigmat(p)
    } )



    pow <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.85), length =200, treatment.difference=0.12))

    t1e <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = 200, prob.range = c(0,0.9), length = 200, treatment.difference = 0))


    cover <-  mclapply(posteriors,
                       function(pr) calc.coverage(posterior=pr, level = 0.95, n.control = 200, smooth = 0.03), mc.cores=30)

    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("Fix-5/study_fix_",i,".rda"))

  })
}
