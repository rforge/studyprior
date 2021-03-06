
library(foreach)
library(StudyPrior)
inla.setOption(num.threads=2)
recalc <- 1:1000
library(parallel)

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
  F.CR9 <- binom.prior("PP.Cor", x = x, n=n, d.prior.cor=0.9, samples=5000, length=512)
  # F.C95 <- binom.prior("PP.Cor", x = x, n=n, d.prior.cor=0.95, samples=5000, length=512)
  F.SUM <- binom.prior("PP.FB", x = sum(x), n=sum(n), samples=5000, length=512, mc.cores=1, verbose=FALSE)
  F.ROB <- conj.approx2(distr = binom.prior("MAP.FB", x = x, n=n),
                       type = "beta",
                       robust = 0.1)
  save(F.CR9,F.SUM,F.ROB, n, x,
       file=paste0("Fix-5/models_g_",i,".rda"))
})

mclapply(mc.cores=35, recalc, function(i){
  print(i)

  set.seed(300*i)
  N <- 5
  n <- rep(50,N)
  z <- 0.65
  x <- mapply(rbinom, size=n, n=1, prob=z)
  x/n
  Xs <- rbinom(1,100,0.6)
  Ns <- 200
  F.PP0 <- binom.PP.FB.COR(x,n,mixture.size = 8^5, d.prior.cor = 0)
  F.PP5 <- binom.PP.FB.COR(x,n,mixture.size = 8^5, d.prior.cor = 0.5)
  F.PP1 <- binom.PP.FB.COR(x,n,mixture.size = 1000, d.prior.cor = 1)
  save(F.PP0,F.PP5,F.PP1, n, x,
       file=paste0("Fix-5/models_h_",i,".rda"))
})

#
# mclapply(mc.cores=25, recalc, function(i){
#   load(file=paste0("Fix-5/models_f_",i,".rda"))
#   print(i)
#
#   set.seed(300*i)
#   N <- 5
#   n <- rep(50,N)
#   z <- 0.65
#   x <- mapply(rbinom, size=n, n=1, prob=z)
#   x/n
#   Xs <- rbinom(1,100,0.6)
#   Ns <- 200
#   cat(i,".\n")
#   F.MFP <- binom.prior("MAP.FB", x = x, n=n)
#   cat(i,"..\n")
#   F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1, verbose=FALSE)
#   F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=512, mc.cores=1, verbose=FALSE)
#   save(F.MFP, F.MEP,F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL, n,x,
#        file=paste0("Fix-5/models_f_",i,".rda"))
# }
# )


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



# library(foreach)
# library(doParallel)

####################################################
CORES <- 20
recalc <- 1:1000

for(i in recalc){
  print(i)
  try({
    load(file=paste0("Fix-5/models_f_",i,".rda"))

    Ns <- 200

    posteriors <-
    lapply(list(
      F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL),
      Calc.posterior.all, N=Ns, mc.cores=CORES
    )

    posteriors2 <-
      lapply(list(
        F.MFP),
        Calc.posterior.all, N=Ns, mc.cores=CORES
      )

    mse <- lapply(posteriors,
      function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns))


    tc <- system.time(
      bias <- lapply(posteriors,
        function(pr) {
         print("Starting!")
           calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES)
        }
    ))

    do.sigmat <- function(pr) sig.matrix(posterior = pr, n.control=200, n.treatment = 200,
                                         check.xt=00:200, check.xs=0:200,
                                         level=0.975, mc.cores=CORES, debug=FALSE)


    SIGMAT <- lapply(posteriors, function(p){
      print("SIGMATING!")
        do.sigmat(p)
    } )



    pow <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = 200,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12),
                    mc.cores=CORES)

    t1e <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = 200,
                                           prob.range = c(0,0.9), length = 200, treatment.difference = 0),
                    mc.cores=CORES)


    cover <-  mclapply(posteriors,
                       function(pr) calc.coverage(posterior=pr, level = 0.95, n.control = 200, smooth = 0.03),
                       mc.cores=CORES)

    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("Fix-5/study_fix_",i,".rda"))

  })
}



# For new priors  -------------------------------------------------
CORES <- 1
recalc <- 1:1000
# recalc <- 501:750
# recalc <- 751:1000

calc <- function(i){
  print(i)
  try({
    load(file=paste0("Fix-5/models_g_",i,".rda"))

    Ns <- 200

    posteriors <-lapply(list(F.CR9,F.SUM),
                          Calc.posterior.all, N=Ns, mc.cores=CORES)

    mse <- lapply(posteriors,
                  function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns))

    mse <- c(mse, list(
      calc.MSE.mean(prior=F.ROB, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns)
    ))

    tc <- system.time(
      bias <- c(lapply(posteriors,
                     function(pr) {
                       print("Starting!")
                       calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES)
                     }
      ), list(calc.bias(prior = F.ROB, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES))
      )
      )

    do.sigmat <- function(pr) sig.matrix(posterior = pr, n.control=200, n.treatment = 200,
                                         check.xt=00:200, check.xs=0:200,
                                         level=0.975, mc.cores=CORES, debug=FALSE)


    SIGMAT <- c(lapply(posteriors, function(p){
      print("SIGMATING!")
      do.sigmat(p)
    } ),
    list(sig.matrix(prior = F.ROB, n.control=200, n.treatment = 200,
                    check.xt=00:200, check.xs=0:200,
                    level=0.975, mc.cores=CORES, debug=FALSE))
    )



    pow <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = 200,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12),
                    mc.cores=CORES)

    t1e <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = 200,
                                           prob.range = c(0,0.9), length = 200, treatment.difference = 0),
                    mc.cores=CORES)


    cover <-  c(mclapply(posteriors,
                       function(pr) calc.coverage(posterior=pr, level = 0.95, n.control = 200, smooth = 0.03),
                       mc.cores=CORES),
                list(calc.coverage(prior = F.ROB, level = 0.95, n.control = 200, smooth = 0.03)))

    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("Fix-5/study_gix_",i,".rda"))

  })
}



mclapply(recalc, calc, mc.cores=40)



# For new priors  h -------------------------------------------------
CORES <- 1
recalc <- 1:1000


calc <- function(j){
  print(j)
  try({
    load(file=paste0("Fix-5/models_h_",j,".rda"))

    Ns <- 200

    priors <- list(F.PP0 = environment(F.PP0)$mix,
                   F.PP1 = environment(F.PP1)$mix,
                   F.PP5 = environment(F.PP5)$mix)

    cover <- t1e <- pow <- SIGMAT <- bias <- mse <- list()

    for(i in 1:3){
      system.time(
        posteriors <- lapply(priors[i], function(pr){
        mclapply(0:200, posterior.mixture.prior, ns=200, mixture.prior=pr)
      })
  )
      system.time(
      mse[i] <- lapply(posteriors, function(pr)
        calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns)
      )
      )
      system.time(
      bias[i] <- lapply(posteriors, function(pr)
        calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES))

      )
      system.time(

      SIGMAT[i] <- c(lapply(posteriors, function(pr)
        sig.matrix(posterior = pr, n.control=200, n.treatment = 200,
                   check.xt=0:200, check.xs=0:200,
                   level=0.975, mc.cores=CORES))
      ) )
      system.time(

      pow[i] <- mclapply(SIGMAT[i],
                      function(S) calc.power(sig.mat=S, n.binom.control = 200,
                                             prob.range = c(0,0.85), length =200, treatment.difference=0.12),
                      mc.cores=CORES)
      )
      system.time(
      t1e[i] <- mclapply(SIGMAT[i],
                      function(S) calc.power(sig.mat=S, n.binom.control = 200,
                                             prob.range = c(0,0.9), length = 200, treatment.difference = 0),
                      mc.cores=CORES)

      )
      system.time(
      cover[i] <-  lapply(posteriors, function(pr)
        calc.coverage(posterior = pr, level = 0.95, n.control = 200, smooth = 0.03))
      )
      rm(posteriors)
      gc()

    }


    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("Fix-5/study_hix_",j,".rda"))

  })
}



mclapply(recalc, calc, mc.cores=30)
