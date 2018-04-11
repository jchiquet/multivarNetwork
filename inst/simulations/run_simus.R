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

## model parameters
p <- 40
hardeness <- .1
nedges <- p
seq_n <- list(expression(p/2),expression(p),expression(2*p))
seq_K <- c(2, 3, 4)
## simu parameters
cores <- 10
nsim <- 100

agreements <- c("anti", "full", "indep")

results <- as.tibble(
Reduce("rbind", lapply(agreements, function(agreement) {
  Reduce("rbind", lapply(seq_K, function(K) {
    cat("\n\t -agreement:", agreement, "K =",K, "\n\n")
    Reduce("rbind", pbmclapply(1:nsim, function(sim) {
      G <- rgraph(p, nedges)
      Reduce("rbind", lapply(seq_n, function(n) {
        cat("\n - n =", eval(n)," ")
        Xmult <- rmultivar(eval(n), K, G, agreement, hardeness)$X
        Xuniv <- Reduce("rbind", Xmult)
        univar.indep <- lapply(Xmult, multivarNetwork, select="none")
        univar.merge <- multivarNetwork(Xuniv, select="none")
        multivar     <- multivarNetwork(Xmult, select="none")
        data.frame(
          univar.indep = mean(sapply(univar.indep, function(univ) perf.auc(perf.roc(univ$networks, G)))) ,
          univar.merge = perf.auc(perf.roc(univar.merge$networks, G)),
          multivar     = perf.auc(perf.roc(multivar$networks, G)),
          scenario = paste0("n=",as.character(n)), K = K, sim = sim, agreement=agreement
        )
      }))
    }, mc.cores = cores))
  }))
  }))
) %>%
  gather(key = "method", value = "AUC", -scenario, -K, -sim, -agreement) %>%
  mutate(scenario = fct_inorder(factor(scenario)), K = factor(K), method = fct_inorder(factor(method)),
         agreement = fct_inorder(factor(agreement))) %>%
  mutate(agreement = fct_recode(agreement,
                    "same intra, same inter" = "full",
                    "same intra, no inter"   = "indep",
                    "no intra, same inter"   = "anti"),
         method = fct_recode(method,
                    "separate"       = "univar.indep",
                    "merge"          = "univar.merge",
                    "multiattribute" = "multivar") )

task_status <- c("2" = "K=2", "3" = "K=3", "4" = "K=4")

pAUC <- ggplot(results, aes(x=scenario, fill=method, y=AUC)) +
  geom_boxplot(notch = TRUE) + # geom_jitter(alpha=0.1, size=1) +
  facet_grid(agreement~K, scales = "free_y", labeller=labeller(K=task_status)) + scale_fill_viridis(discrete=TRUE, option="plasma")+ theme_bw(base_size = 15)
ggsave(plot = pAUC, filename ="../../chapter/figures/res_simu_new.pdf", width = 10, height=7)

save(results, file=paste0("simu_",Sys.Date(),".RData"))
