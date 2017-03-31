RNA-seq Generation
================

-   [How to install](#how-to-install)
-   [RNA Random expression](#rna-random-expression)
    -   [Normal](#normal)
    -   [Log Normal](#log-normal)
    -   [Real gene expression](#real-gene-expression)

How to install
==============

``` r
devtools::install_github('averissimo/gen.rna.seq.data')
```

RNA Random expression
=====================

Normal
------

``` r
my.random <- function(len, mean = runif(1, min = 0, max= 1), sd = runif(1, min = .5, max= 1.5)) { 
  rnorm(len, mean = mean, sd = sd)
}

len <- 100
dat.normal <- data.frame()
for (ix in 1:plot.num) {
  meanlog  <- runif(1, min = 0, max= 1)
  sdlog <- runif(1, min = .5, max= 1.5)
  
  new.meas <- data.frame(meas = my.random(len, meanlog, sdlog),
                         ix.me = array(sprintf('set-%03d', ix), dim = len),
                         fun = array('mean + sd', len))
  dat.normal <- rbind(dat.normal, new.meas)
}

ggplot(dat.normal, aes(meas, colour = ix.me)) + 
  theme_minimal() + theme(legend.position = 'none') +
  geom_freqpoly(bins = 200) + expand_limits(x = 0) +
  facet_wrap( ix.me~ fun, ncol = 3, scales = 'free',
              labeller = function(...) {label_value(..., multi_line = FALSE)}) +
  ggtitle('Generated gene expression (Gaussian)')
```

![](/home/averissimo/work/rpackages/gen.rna.seq/README_files/figure-markdown_github/gaussian-1.png)

Log Normal
----------

``` r
my.random <- function(len, 
                      meanlog = log(runif(1, min = 0, max= 3)), 
                      sdlog = log(1 + runif(1, min = 0.08, max= 1.5)),
                      transl = runif(1, min = 0, max= 1)) { 
  res <- rlnorm(len, meanlog = meanlog, sdlog = sdlog)
  if (meanlog <= log(1e-1)) {
    # higher change if low expressed genes
    res[runif(length(res)) >= 0.9] <- 0
    sdlog <- log(1) + sdlog / 10
  }
  #mean(res) * transl + res
  return(res)
}

len <- 100
dat.log <- data.frame()
for (ix in 1:plot.num) {
  new.meas <- data.frame(meas = my.random(len),
                         ix.me = array(sprintf('set-%03d', ix), dim = len),
                         fun = array('', len))
  dat.log <- rbind(dat.log, new.meas)
}

ggplot(dat.log, aes(meas, colour = ix.me)) + 
  theme_minimal() + theme(legend.position = 'none') +
  geom_freqpoly(bins = 200) + expand_limits(x = 0) +
  facet_wrap( ix.me~ fun, ncol = 2, scales = 'free',
              labeller = function(...) {label_value(..., multi_line = FALSE)}) +
  ggtitle('Generated gene expression (Log-normal)')
```

![](/home/averissimo/work/rpackages/gen.rna.seq/README_files/figure-markdown_github/log-normal-1.png)

Real gene expression
--------------------

Using breast cancer data

``` r
library(brca.data)
sample.data <- fpkm.per.tissue$solid.tissue.normal
sample.data.sd <- sapply(1:nrow(sample.data), function(ix) {sd(sample.data[ix,])})
sample.data <- sample.data[sample.data.sd  > mean(sample.data.sd),]
sample.data.subset <- sample.data[sample(seq(nrow(sample.data)),plot.num),]
sample.data.mel <- melt(sample.data.subset)
#
ggplot(sample.data.mel, aes(value, colour = Var1)) + 
  theme_minimal() + theme(legend.position = 'none') +
  geom_freqpoly(bins = 200) + expand_limits(x = 0) +
  facet_wrap(~ Var1,
             scales = 'free',
             ncol = 2, 
             labeller = function(...) {label_value(..., multi_line = FALSE)}) +
  ggtitle('Gene expression from BRCA')
```

![](/home/averissimo/work/rpackages/gen.rna.seq/README_files/figure-markdown_github/brca-1.png)
