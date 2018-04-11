rm(list=ls())
library(tidyverse)
library(viridis)

load(file="nci60_prot+expr-2017-10-05.RData")

image.NA <- function(z,  zlim=c(0,1), col=c("white", "midnightblue"), na.color='red', outside.below.color='black', outside.above.color='white',...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
  
  par(mar=c(0.1,0.1,0.1,0.1))
  image(z[nrow(z):1,],  zlim=zlim, col=col, xaxt="n", yaxt="n") # we finally call image(...)
  box()
  par(mar=c(5.1,4.1,4.1,2.1))
}

pdf(file="../../chapter/figures/protNet_NCI60.pdf")
image.NA(as.matrix(protNet.1se$networks[[1]]))
dev.off()
pdf(file="../../chapter/figures/exprNet_NCI60.pdf")
image.NA(as.matrix(exprNet.1se$networks[[1]]))
dev.off()
pdf(file="../../chapter/figures/bivarNet_NCI60.pdf")
image.NA(as.matrix(bivarNet.1se$networks[[1]]))
dev.off()


nedges.prot <- sapply(protNet.all$networks, sum)
nedges.expr <- sapply(exprNet.all$networks, sum)
nedges.mult <- sapply(bivarNet.all$networks, sum)

max.edges <- 1500
J12 <-vector("numeric", max.edges)
J13 <-vector("numeric", max.edges)
J23 <-vector("numeric", max.edges)
for (nedges in 1:max.edges) {
  prot <- which(protNet.all$networks[[match(TRUE, nedges.prot >= 2*nedges)]] !=0)
  expr <- which(exprNet.all$networks[[match(TRUE, nedges.expr >= 2*nedges)]] !=0)
  mult <- which(bivarNet.all$networks[[match(TRUE, nedges.mult >= 2*nedges)]] !=0)

  J12[nedges] <- length(intersect(prot,expr)) / length(union(prot,expr))

  J13[nedges] <- length(intersect(prot,mult)) / length(union(prot,mult))

  J23[nedges] <- length(intersect(expr,mult)) / length(union(expr,mult))

}

dp <- data.frame(jaccard = c(J12, J13, J23),
           couple = rep(c("protein/expr.","protein/expr+protein","expr/expr+protein"), each=max.edges),
           edges = rep(1:max.edges, 3))
p <- ggplot(dp, aes(x=edges, y=jaccard, group=couple, colour=couple)) + geom_line() +
  theme_bw(base_size = 20) + scale_color_viridis(discrete=TRUE) + ylim(0.,0.4)
ggsave(plot = p, filename = "../../chapter/figures/jaccard_NCI60.pdf", width=10, height=7)

library(blockmodels)
mnet <- bivarNet.1se$networks[[1]]
## Adjust a SBM to the inferred CIT network
my_model <- BM_bernoulli("SBM_sym",as.matrix(mnet))
my_model$estimate(); SBM <- list()
SBM$cl <- apply(my_model$memberships[[which.max(my_model$ICL)]]$Z, 1, which.max)
SBM$pi <- my_model$model_parameters[[which.max(my_model$ICL)]]$pi

library(igraph)
library(RColorBrewer)
pal   <- brewer.pal(10, "Set3")

g <- graph_from_adjacency_matrix(mnet, mode = "undirected", weighted = NULL, diag = FALSE)
V(g)$class <- SBM$cl
V(g)$size <- 5
V(g)$frame.color <- "white"
V(g)$color <- pal[V(g)$class]
V(g)$label <- ""
E(g)$arrow.mode <- 0


