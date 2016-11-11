
#' Coverage using prior using the Bayari thing
#'
#' @param prior
#' @param length
#'
#' @return
#' @export
#'
#' @examples
calc.coverage <- function(prior, level, n.control, smooth, ...){


  CIs <-  sapply(0:n.control, function(Xs){


      if(inherits(prior,"function")){
      post <- function(p,k=1) prior(p, X=Xs)*dbinom(x=Xs, size=n.control, prob=p)/k
      # print(Xs)
      K <- integrate(post, lower=0, upper=1)$value

    } else if (inherits(prior,"mixture.list")){
      post <-function(p) eval.fun.list(p, posterior.fun.list(Xs, n.control, prior))
    } else if (inherits(prior,"list")){
      post <-function(p) eval.fun.list(p, posterior.fun.list(Xs, n.control, prior[[Xs+1]]))
    }

    # get the confidence intervals
    P <- seq(0,1,length.out = 1024)
    marg <- list(x=P, y=post(P))
    CI <- inla.hpdmarginal(level, marg)

    # apply the smoothing
    return(CI)
  })
  CI.mat <- as.matrix(t(CIs))
  plot.coverage(confidence.intervals = CI.mat,
                smooth=smooth,
                prob = level, ...)

}


plot.coverage <- function(confidence.intervals,dx,dn, # Matrix mit 2 Spalten, die die unteren
                          # und oberen Intervallgrenzen enthaelt und Zeilenzahl n+1 hat.
                          add = FALSE,  # soll der Plot hinzugefügt werden?
                          smooth = NULL, # falls nicht NULL, wird eine lokal gemittelte Überdeckungswkeit
                          # mit dem epsilon-Parameter smooth eingezeichnet
                          faktor = 20,  # legt die Zahl der p-Werte, für die das Coverage berechnet wird,
                          # als n * faktor - 1 fest
                          prob = 0.95,  # nominales Konfidenzniveau,
                          do.plot = FALSE,
                          ...
)
{
  # zu hohe Genauigkeiten (Wilson) vermeiden
  confidence.intervals = round(confidence.intervals, 12)
  n = nrow(confidence.intervals)-1
  x = 0:n
  # in welcher Spalte befindet sich die untere Grenze?
  lower.ind = which.min(confidence.intervals[1,])
  upper.ind = setdiff(c(1,2), lower.ind)
  # pvector berechnen : ohne 0 und 1
  pvector = seq(0, 1, length = faktor*n + 1)[-c(1, faktor*n+1)]
  gerade = !(length(pvector) %% 2)
  # pvector.halb = pvector[1:ceiling(length(pvector)/2)]
  # Überdeckungswahrscheinlichkeiten berechnen:
  # nur für eine Hälfte, da Symmetrie um p = 0.5
  coverage.probs = sapply(pvector,
                          function(p){
                            sum(dbinom(x, n, p)[which(p >= confidence.intervals[,lower.ind] & p <= confidence.intervals[,upper.ind])])
                          }
  )
  # if(gerade)
  #   coverage.probs = c(coverage.probs)#, rev(coverage.probs))
  # else
  #   coverage.probs = c(coverage.probs)[-1]#, rev(coverage.probs)[-1])
  # zeichnen
  if(do.plot){
  if(!add){
    plot(pvector, coverage.probs, type = "l", xlab = expression(paste("True ", pi)),
         ylab = "Coverage probability", col = "darkgrey", ...)
    abline(h = prob, lty = 2)
  }
  else
  {
    lines(pvector, coverage.probs, type = "l", col = "darkgrey", ...)
  }}
  # evtl. Smooth einzeichnen (vgl. Bayarri und Berger)
  if(!is.null(smooth)){
    # allgemeine a-Funktion
    a = function(p){
      if(p==0) NA
      else if (p==1) NA
      else
        if(p <= smooth)
          NA#(p*(1-p)*p^(-2) - 1)*p#1 - 2*smooth
      else if(p >= 1 - smooth)
        NA#(p*(1-p)*(1-p)^(-2) - 1)*p #1/smooth - 3 + 2*smooth
      else
        (p*(1-p)*smooth^(-2) - 1)*p
    }
    # Funktion zur Berechnung des local coverage
    local.coverage = function(p){
      ap = a(p)
      a1mp = a(1-p)
      alpha = ap + x
      beta = a1mp + n - x
      values.gamma = (lchoose(n, x)
                      + lgamma(ap + a1mp) - lgamma(ap) - lgamma(a1mp)
                      + lgamma(ap + x ) + lgamma(a1mp + n - x) - lgamma(ap + a1mp + n)
      )
      values.integral = log(
        pbeta(confidence.intervals[,upper.ind], alpha, beta) - pbeta(confidence.intervals[,lower.ind], alpha, beta)
      )
      sum(exp(values.gamma + values.integral))
    }

    # berechnen und einzeichnen:
    coverage.average = sapply(pvector, local.coverage)

    if(gerade)
      coverage.average = c(coverage.average)#, rev(coverage.average))
    else
      coverage.average = c(coverage.average)#, rev(coverage.average)[-1])

    if(do.plot) lines(pvector, coverage.average, lwd = 2)
    names(coverage.average) <- pvector
    return(coverage.average)
  }
}
