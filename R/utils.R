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

##
## This function was lifted directly from the IRangesOverview.pdf
##
plotIRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                        col = "black", sep = 0.5, ...) {
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins) * (height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom +
       height, col = col, ...)
  title(main)
  axis(1)
}

##
## This function was based on the IRangesOverview.pdf
## 
coverageDensity <- function(x, as.prob = FALSE) {
  cov <- coverage(x)
  cov <- as(cov, "vector")
  if (as.prob) cov/sum(cov) else cov
}  

##
## Add a character method for extracting by names.
## 
setMethod("[", "DNAStringSet", function(x, i, j, ..., drop) {
  if (!is.character(i))
    callNextMethod(x, i) # , j, ..., drop)
  else {
    idx <- do.call(c, lapply(i, function(a) match(a, names(x))))
    x[idx]
  }
})

makeMap <- function(dta, xname, ynames) {
  stopifnot(all(c(xname, ynames) %in% colnames(dta)))
  
  x <- dta[,xname]
  y <- dta[,ynames]
  w <- !duplicated(x)
  
  if (is.null(dim(y))) {
    y <- y[w]
    names(y) <- x[w]
  } else {
    y <- y[w, ]
    rownames(y) <- x[w]
  }
  return(y)
}

readM4 <- function(fname, parseReadNames = FALSE, ...) {
  cols <- c("queryId", "targetId", "smithWatermanEditScore", "hmmProbScore", "smithWatermanQVScore",
            "pctIdentity", "queryStrand", "queryStart", "queryEnd", "queryLength", "targetStrand",
            "targetStart", "targetEnd", "targetLength")
  tbl <- read.table(fname, stringsAsFactors = TRUE, ...)
  colnames(tbl) <- cols

  if (parseReadNames)
    tbl <- cbind(tbl, parseReadNames(tbl$queryId))
  
  return(tbl)
}

readM5 <- function(fname, parseReadNames = FALSE, ...) {
  cols <- c("queryId", "query", "read.start", "read.end",
            "read.strand", "ref", "ref.span", "ref.start", "ref.end",
            "ref.strand", "score", "match", "mismatch", "insert", "deletion",
            "read.seq", "aln", "ref.seq")
  tbl <- read.table(fname, stringsAsFactors = FALSE, ...)
  colnames(tbl) <- cols

  if (parseReadNames)
    tbl <- cbind(tbl, parseReadNames(tbl$queryId))

  return(tbl)
}

parseReadNames <- function(rn) {
  d <- as.data.frame(do.call(rbind, lapply(strsplit(rn, "_"), function(r) {
    r[1] <- gsub("x", "", r[1])
    r[2] <- gsub("y", "", r[2])
    r
  })))
  colnames(d) <- c("x", "y", "runid", "moviename", "timestamp", "machine", "pane",
                   "block_fragment_positioninoriginal")
  d$block <- sapply(strsplit(as.character(d[,8]), "\\."), "[", 1)
  d$fragment <- sapply(strsplit(as.character(d[,8]), "\\."), "[", 2)
  d <- d[,-8]
  return(d)
}


calcN50 <- function(v, pure = FALSE) {
  if (pure) {
    median(rep(v, v))
  } else {
    a <- sort(v)
    y <- cumsum(a)
    sum(y/sum(y) * a)
  }
}
 
cleanFASTA <- function(ifile, ofile = NULL, removeSpaces = NULL) {
  fsta <- readFASTA(ifile, strip.descs = TRUE)
  seqs <- sapply(fsta, "[[", "seq")
  desc <- sapply(fsta, "[[", "desc")

  if (!is.null(removeSpaces)) {
    stopifnot(is.character(removeSpaces))
    desc <- sapply(strsplit(desc, split = " "), paste, collapse = removeSpaces)
  }
  seqs <- as.character(seqs)
  names(seqs) <- desc
  xstring <- DNAStringSet(seqs)
  
  if (!is.null(ofile)) {
    write.XStringSet(xstring, file = ofile)
  }
  
  return(xstring)
}

plotDensity <- function(x, col = 1:length(x), legend = FALSE, xlim = NULL, ylim = NULL, log = NULL,
                        lwd = rep(1, length(x)),
                        lty = rep(1, length(x)), ...) {
  LOG <- FALSE
  
  if (! is.null(log)) {
    if (log == 'y')
      stop("Y access log not yet implemented.")
    else {
      LOG <- TRUE
      x   <- lapply(x, log10)
    }
  }

  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }
  
  a <- lapply(x, density, na.rm = TRUE)
  r <- sapply(a, function(b) c(range(b$x), range(b$y)))
  
  if (is.null(xlim)) xlim <- c(min(r[1,]), max(r[2,]))
  if (is.null(ylim)) ylim <- c(min(r[3,]), max(r[4,]))

  plot(NA, xlim = xlim, ylim = ylim, ylab = "density", xaxt = if (LOG) 'n' else 's', ...)
  
  mapply(function(d, col, lty, lwd) {
    points(d$x, d$y, type = 'l', col = col, lty = lty, lwd = lwd, ...)
  }, a, col, lty, lwd)

  if (LOG) {
    ## Steal a page out of ggplot.
    u <- par('usr')[1:2]
    p <- pretty(r, 8)
    e <- parse(text = sprintf("expression(%s)", paste(paste("10^", p, sep = ""), collapse = ",")))
    axis(1, p, labels = eval(e))
  }
  
  if (legend)
    legend("topright", fill = col, lty = lty, names(x))
}

##
## Produce a subsampled qqPlot with colorized points based on
## quantiles
##
qqPlot <- function(x, y, twoSided = TRUE, maxPoints = 5000, qtiles = c(.95, .99, .999),
                   pch = 16, colors = c("red", "violet", "orange"), legend = T, ...) {
  stopifnot(length(qtiles) == length(colors))

  m <- min(length(x), length(y))
  
  if (m > maxPoints) {
    m <- maxPoints
  }
  stopifnot(m <= length(x) && m <= length(y))
  
  x <- sort(sample(x, m))
  y <- sort(sample(y, m))

  if (twoSided) {
    qtiles <- c(rev((1 - qtiles)/2), (1 + qtiles)/2)
    colors <- c(rev(colors), "black", colors)
    qbins <- paste(paste(c(0, qtiles)*100, c(qtiles, 1)*100, sep = '-'), "%", sep = "")
    bincols <- colors
  } else {
    colors <- c("black", colors)
    qbins <- paste(paste(c(0, qtiles*100), c(qtiles, 1)*100, sep = "-"), "%", sep = "")
    bincols <- colors
  }
  colors <- colors[findInterval(seq(0, 1, length = m), qtiles) + 1]
  plot(x, y, col = colors, pch = pch, asp = 1, ...)
  abline(0, 1, lwd = 3, col = "grey", lty = 3)
  if (legend) legend("topleft", qbins, fill = bincols)
}

qqPairs <- function(x, twoSided = FALSE, ...) {
  if (class(x) == "list") {
    m <- min(sapply(x, length))
    m <- if (m > 1e4) 1e4 else m
    x <- sapply(x, quantile, prob = seq(0, 1, length = m), na.rm = TRUE)
  }
  u <- range(apply(x, 2, function(z) range(z[is.finite(z)])))

  if (ncol(x) > 2) {
    pairs(x, lower.panel = NULL, upper.panel = function(x, y) {
      par(new = TRUE)
      qqPlot(x, y, xlim = u, ylim = u, twoSided = twoSided, ...)
    })
  } else {
    qqPlot(x[,1], x[,2], xlim = u, ylim = u, twoSided = twoSided,
           xlab = colnames(x)[1], ylab = colnames(x)[2], ...)
  }
}

maPlot <- function(x, y = NULL, log = TRUE, ...) {
  if (class(x) == "matrix" && ncol(x) == 2) {
    y <- x[,2]
    x <- x[,1]
  }
  if (log) {
    x <- log(x)
    y <- log(y)
  }
  plot(xx <- (x+y)/2, yy <- y - x, ...)
  invisible(cbind(xx, yy))
}

interleave <- function(..., list = NULL) {
  if (!is.null(list))
    a <- list
  else 
    a <- list(...)
  
  stopifnot(length(unique(sapply(a, class))) == 1)
  
  doList <- function() {
    do.call(c, lapply(1:length(a[[1]]), function(i) {
      do.call(c, lapply(a, function(b) b[i]))
    }))
  }
  doMatrix <- function() {
    do.call(cbind, lapply(1:ncol(a[[1]]), function(i) {
      do.call(cbind, lapply(a, function(b) b[,i]))
    }))
  }
  x <- switch(class(a[[1]]), "matrix" = doMatrix(), "data.frame" = doMatrix(), "list" = doList())

  return(x)
}

plotDwellTimePDF <- function(x, col = 1:length(x), xlim = NULL, ylim = NULL, main = NULL,
                            xlab = NULL, ylab = NULL, lwd = 2, transform = function(a, b) (a * 10^b), ...) {

  ds <- lapply(x, function(b) {
    d <- density(log10(b), na.rm = TRUE)
    list(x = d$x, y = transform(d$y, d$x))
  })
  
  xlab <- if (is.null(xlab)) "Advance Time" else xlab
  ylab <- if (is.null(ylab)) "Advance Time * PDF" else  ylab
  main <- if (is.null(main)) "Advance Time Distribution" else main
  
  xrange <- range(sapply(ds, function(b) range(b$x, na.rm = TRUE, finite = TRUE)))
  yrange <- range(sapply(ds, function(b) range(b$y, na.rm = TRUE, finite = TRUE)))
  
  plot(NA, main = main, xlab = xlab, ylab = ylab,
       xlim = if(is.null(xlim)) xrange else xlim,
       ylim = if(is.null(ylim)) yrange else ylim,
       xaxt = 'n', yaxt = 'n')

  mapply(function(xy, color) {
    lines(supsmu(xy$x, xy$y), col = color, lwd = lwd, ...)
  }, ds, col)
  
  ## xaxis.
  p <- pretty(do.call(c, lapply(ds, function(a) a$x)), 8)
  e <- parse(text = sprintf("expression(%s)", paste(paste("10^", p, sep = ""), collapse = ",")))
  axis(1, p, labels = eval(e))
  
  ## yaxis
  p <- pretty(do.call(c, lapply(ds, function(a) a$y)), 8)
  axis(2, sqrt(p), p)

}

readGFF <- function(gffFile, attrClasses = c(), keepAttributes = FALSE, sep = "\t",
                  fill = TRUE, flush = TRUE) {
  ##
  ## Based on the description here:
  ## http://www.sanger.ac.uk/resources/software/gff/spec.html
  ## 
  d <- read.table(gffFile, comment.char = '#', stringsAsFactors = FALSE, sep = sep,
                  fill = fill, flush = flush)
  colnames(d) <- c('seqid','source','type','start','end','score','strand','frame','attributes')

  attrs <- strsplit(d$attributes, split = ";")
  keyValues <- lapply(attrs, function(a) strsplit(a, split = "="))
  attrNames <- unique(unlist(lapply(keyValues, function(a) sapply(a, "[[", 1))))
  flat <- lapply(keyValues, function(kv) {
    l <- vector("list", length(attrNames))
    l[1:length(l)] <- NA
    names(l) <- attrNames
    l[match(sapply(kv, "[[", 1), names(l))] <- sapply(kv, "[[", 2)
    return(l)
  })
  w <- if(keepAttributes) {
    w <- 1:ncol(d)
  } else {
    -which(colnames(d) == "attributes")
  }
  nd <- cbind(d[,w], do.call(rbind, lapply(flat, unlist)), stringsAsFactors = FALSE)

  for (i in seq.int(attrClasses)) {
    class(nd[,names(attrClasses)[i]]) <- attrClasses[i]
  }
  return(nd)
}

readVariantsGFF <- function(gffFile, keepAttributes = TRUE) {
  a <- readGFF(gffFile, attrClasses = c("coverage" = "integer", "confidence" = "numeric",
                          "length" = "integer"), keepAttributes = keepAttributes)
  a <- a[, !(colnames(a) %in% c("source", "frame"))]
  a$x <- (a$start + a$end)/2
  a$length[a$type == "SNV"] <- 1
  a$genotype[a$type == "deletion"] <- '-'
  return(a)
}

readAlignmentSummaryGFF <- function(gffFile, keepAttributes = TRUE) {
  a <- readGFF(gffFile, attrClasses = c("ins" = "integer", "del" = "numeric", "snv" = "integer"),
               keepAttributes = keepAttributes)
  a <- a[, !(colnames(a) %in% c("source", "frame"))]

  pCommaSep <- function(nm, colnames) {
    x <- as.data.frame(do.call(rbind, strsplit(a[,nm], ",")))
    x[,] <- as.numeric(as.matrix(x))
    colnames(x) <- colnames
    return(x)
  }
  a <- cbind(a, pCommaSep("cov2", c("v.mean", "v.sd")))
  a <- cbind(a, pCommaSep("cov", c("v.min", "v.med", "v.max")))
  a <- cbind(a, pCommaSep("gaps", c("v.ngaps", "v.gapsize")))
  a$x <- (a$start + a$end)/2
  return(a)
}

namedRange <- function(range, names) {
  names(range) <- names
  return(range)
}

collapse <- function(x, colname = NA) {
  stopifnot(is.list(x))

  x <- Filter(function(z) !is.null(z), x)
  
  if (! is.null(dim(x[[1]]))) {
    lf <- nrow
    cf <- rbind
    on <- colnames(x[[1]])
  } else {
    lf <- length
    cf <- c
    on <- "value"
  }
  nx <- rep(names(x), sapply(x, lf))
  dx <- as.data.frame(do.call(cf, x), row.names = NULL)

  if (is.na(colname)) {
    i <- 1
    while(paste("L", i, sep = "") %in% on) { 
      i <- i + 1
    }
    ncname <- paste("L", i, sep = "")
  } else {
    if (colname %in% on) warning("colname already exists, setting anyway.")
    ncname <- colname
  }
  dx[[ncname]] <- nx
  colnames(dx) <- c(on, ncname)
  rownames(dx) <- NULL
  return(dx)
}

toPhred <- function(ep, max) {
  -10 * log10(ifelse(ep == 0, 1/max, ep))
}

fromPhred <- function(qv, max) {
  1/(10^(qv/10))
}

readSam <- function(fileName, ...) {
  ## columns.
  ##
  ## 1	QNAME	Query template/pair NAME
  ## 2	FLAG	bitwise FLAG
  ## 3	RNAME	Reference sequence NAME
  ## 4	POS     1-based leftmost POSition/coordinate of clipped sequence
  ## 5	MAPQ	MAPping Quality (Phred-scaled)
  ## 6	CIAGR	extended CIGAR string
  ## 7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
  ## 8	MPOS	1-based Mate POSistion
  ## 9	TLEN	inferred Template LENgth (insert size)
  ## 10	SEQ     query SEQuence on the same strand as the reference
  ## 11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
  ## 12+  OPT     variable OPTional fields in the format TAG:VTYPE:VALUE

  ## bitmasks.
  ##
  ## 0x0001	p	the read is paired in sequencing
  ## 0x0002	P	the read is mapped in a proper pair
  ## 0x0004	u	the query sequence itself is unmapped
  ## 0x0008	U	the mate is unmapped
  ## 0x0010	r	strand of the query (1 for reverse)
  ## 0x0020	R	strand of the mate
  ## 0x0040	1	the read is the first read in a pair
  ## 0x0080	2	the read is the second read in a pair
  ## 0x0100	s	the alignment is not primary
  ## 0x0200	f	the read fails platform/vendor quality checks
  ## 0x0400	d	the read is either a PCR or an optical duplicate
  l <- length(grep("^@", readLines(fileName)))
  samFile <- read.table(fileName, skip = l, fill = T, stringsAsFactors = FALSE)

  colnames(samFile) <- c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIAGR', 'MRNM',
                         'MPOS', 'TLEN', 'SEQ', 'QUAL',
                         paste("OPT", 1:(ncol(samFile) - 11), sep = "_"))
  return(samFile)
}
