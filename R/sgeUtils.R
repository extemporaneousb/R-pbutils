## Copyright (c) 2010, Pacific Biosciences of California, Inc.

## All rights reserved.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:

##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.

##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.

##     * Neither the name of Pacific Biosciences nor the names of its
##       contributors may be used to endorse or promote products derived
##       from this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED

## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
## BE LIABLE FOR ANY

## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
## GOODS OR SERVICES;

## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
## LIABILITY, OR TORT

## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

.captureContext <- function(startFrame = 1, debug = FALSE) {
  nframes <- sys.nframe()
  res <- new.env(parent = emptyenv())
  
  for (i in startFrame:nframes) {
    p <- parent.frame(i)
    if (debug) {
      cat("In frame :", nframes - i, "\n")
      cat("Capturing:", paste(ls(p), collapse = ", "),  "\n")
    }
    for (v in ls(p)) {
      if (! exists(v, res)) {
        assign(v, get(v, p), res)
      }
    }
  }
  return(res)
}

sgeApply <- function(X, FUN, ..., basedir = ".sge_R", pe = "-pe smp 1", queue = "secondary",
                     debug = FALSE, distribute = TRUE, Rcmd = "R") {

  if (! distribute) {
    return(lapply(X, FUN, ...))
  }
  
  ## start with the calling frame.
  context <- .captureContext(2, debug = debug)
  
  sgePrepare <- function() {
    dir.create(basedir, showWarnings = FALSE)
    save(context, file = paste(basedir, 'context.rda', sep = '/'))
    session <- sessionInfo()
    if (length(session) >= 5) {
      paste(paste("require(", names(session[[5]]), ")\n"), collapse = "")
    }
  }
  
  sgeExecute <- function(exprString) {
    odir <- getwd()
    setwd(basedir)
    
    uid <- paste('obj', round(runif(1, 1e7, 1e8)), sep = "_")
    script <- paste(uid, ".R", sep = "")
    
    scriptTxt <-
      sprintf('load(\"context.rda\") \n attach(context) \n %s = { %s }\n\n save(%s, file = \"%s.rda\") \n', uid,
              exprString, uid, uid)
    scriptTxt <- paste(frontMatter, scriptTxt)
    cat(scriptTxt, file = script)
    submitCmd <- sprintf("qsub -N %s -V -cwd -q %s  %s -b y %s CMD BATCH --vanilla --no-save %s", uid, queue,
                         pe, Rcmd, script)

    print(submitCmd)
    system(submitCmd)
    setwd(odir)
    return(uid)
  }
  
  sgeWait <- function(uids) {
    while (TRUE) {
      if (any(sapply(uids, sgeJobExists))) {
        system("sleep 1")
      } else {
        break
      }
    }
    names(uids) <- uids
    
    lapply(uids, function(a) {
      file <- sprintf("%s/%s.rda", basedir, a)
      if (! file.exists(file)) {
        NULL
      } else {
      load(file)
      get(a)
    }
    })
  }
  
  sgeJobExists <- function(uid) {
    d <- system("qstat -xml", intern = T)
    if (length(grep(uid, d)) > 0) 
      TRUE
    else
      FALSE
  }
  
  sgeCleanup <- function(uid) {
    system(sprintf("rm %s/%s*", basedir, uid))
  }

  sgeLapply <- function(X, FUN, ...) {
    if (! is.atomic(X)) {
      stop("This procedure works only with atomic values of X, i.e., numbers, strings, etc.")
    }
    res <- lapply(X, function(elt) {
      force(elt)
      f <- match.fun(FUN)
      fdef <- sprintf("fff = %s", paste(deparse(f), collapse = "\n"))
      fcall <- sprintf("tryCatch(fff(%s), simpleError = function(se) { print(se) ; return(se) })",
                       paste(deparse(elt), collapse = "\n"))
      paste(fdef, fcall, sep = "\n")
    })
    res <- sgeEvalExprStrings(lapply(res, paste, collapse = "\n"))
    names(res) <- names(X)
    return(res)
  }
  
  sgeEvalExprStrings <- function(exprStrings) {
    uids <- sapply(exprStrings, function(e) sgeExecute(e))
    system("sleep 10")
    res <- sgeWait(uids)
    if (! debug) sapply(uids, sgeCleanup)
    errors <- Filter(function(z) inherits(z, "simpleError"), res)
    if (length(errors) > 0) {
      stop(paste(errors, collapse = '\n'))
    }
    return(res)
  }
  frontMatter <- sgePrepare()
  sgeLapply(X, FUN, ...)
}

