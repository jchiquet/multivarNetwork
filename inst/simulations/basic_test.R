rm(list=ls())
working_directory <-
  switch(Sys.info()["user"],
    "jchiquet" = '~/git/multGGM/code/simulations'
## put your own directory here
  )
setwd(working_directory)

library(multivarNetwork)
library(pbmcapply)
library(tidyverse)
library(forcats)
library(viridis)
source("functions_utils.R")

p <- 20
n <- 20
hardeness <- 0.1
nedges <- p
K <- 4

G <- rgraph(p, nedges)
d <- rmultivar(n, K, G, "full", hardeness)

univar.merge <- multivarNetwork(Reduce("rbind", d$X), mc.cores=10)
univar.indep <- lapply(d$X, multivarNetwork, mc.cores=10)
multivar.net <- multivarNetwork(d$X, mc.cores=10)

print(perf.auc(perf.roc(multivar.net$networks, G)))
print(perf.auc(perf.roc(univar.merge$networks, G)))
print(mean(sapply(univar.indep, function(univ) perf.auc(perf.roc(univ$networks, G)))))

