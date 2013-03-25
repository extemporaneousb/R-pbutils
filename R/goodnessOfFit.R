## 
## Some utility functions for doing goodness of fit tests. These
## binned procedures attempt to shelter goodness of fit tests from
## large outliers. The main issues here are determining the correct
## degrees of freedom.
##
## 

bGOF <- function(v, quantiles) {
  if (length(quantiles) != length(unique(quantiles))) {
    return(NA)
  }
  chisq.test(table(cut(v, quantiles, include.lowest = T)))$statistic
}

bGOFExp <- function(v, nbins = 20, df = nbins - 3, mu = 1/mean(v, na.rm = T)) {
  1-pchisq(bGOF(v, qexp(seq(0, 1, length = nbins), mu)), df)
}

bGOFNormal <- function(v, nbins = 20, df = nbins - 3.5, mu = mean(v, na.rm = T),
                       sigma = sd(v, na.rm = T)) {
  1 - pchisq(bGOF(v, quantiles = qnorm(seq(0, 1, length = nbins), mu, sigma)), df)
}


##
## Here, I calibrate the DF.
##
## hist(replicate(10000, {
##   bGOFExp(rexp(100, 1/rexp(1, 1/10)))
## }), breaks = 20)

## hist(replicate(10000, {
##   bGOFNormal(rnorm(100, rnorm(1), rchisq(1, 2)))
## }), breaks = 20)
