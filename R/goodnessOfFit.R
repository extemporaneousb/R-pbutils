## Copyright (c) 2010, Pacific Biosciences of California, Inc.

## All rights reserved.
 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted (subject to the limitations in the
## disclaimer below) provided that the following conditions are met:
 
##     * Redistributions of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
 
##     * Redistributions in binary form must reproduce the above
##        copyright notice, this list of conditions and the following
##        disclaimer in the documentation and/or other materials provided
##        with the distribution.
  
##     * Neither the name of Pacific Biosciences nor the names of its
##        contributors may be used to endorse or promote products derived
##        from this software without specific prior written permission.
 
## NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
## GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
## BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
## WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
## BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
## BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
## WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
## IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
