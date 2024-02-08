# minimal script to reproduce figure 1 B

library(tidyverse)
library(magrittr)
library(spatstat.geom)
source("alexander_experiments/my_compensate.r")

tfm <- function(x) asinh(x/5)
tb_real = read.csv("alexander_experiments/tb_real.csv")
tb_bead = read.csv("alexander_experiments/tb_bead.csv")
spillover_markers = c("Yb171Di", "Yb172Di", "Yb174Di", "Yb176Di")
target_marker = "Yb173Di"
runmed_k = 1
res = compensate(tb_real, tb_bead, target_marker, spillover_markers, runmed_k)

# reproduce figure 1 A
ggplot(tb_bead) + 
    geom_density(aes(x=tfm(Yb173Di), color=barcode))

# reproduce figure 1 B
ggplot(res$tb_spill_prob) +
    geom_line(aes(x=tfm(Yb173Di), y=spill_prob))

# TODO make sure all the spillover is actualy removed!