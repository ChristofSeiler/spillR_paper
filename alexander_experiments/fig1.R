# minimal script to reproduce figure 1 B

library(dplyr)
library(ggplot2)
library(magrittr)
library(spatstat.geom)
source("alexander_experiments/compensate.R")

tfm <- function(x) asinh(x/5)
tb_real = read.csv("alexander_experiments/tb_real.csv")
tb_bead = read.csv("alexander_experiments/tb_bead.csv")
spillover_markers = c("Yb171Di", "Yb172Di", "Yb174Di", "Yb176Di")
target_marker = "Yb173Di"
runmed_k = 3

# reproduce figure 1 A
ggplot(tb_bead) + 
    geom_density(aes(x=tfm(Yb173Di), color=barcode))

res = compensate(tb_real, tb_bead, target_marker, spillover_markers, runmed_k)

# reproduce figure 1 B
ggplot(res$tb_spill_prob) +
    geom_line(aes(x=tfm(Yb173Di), y=spill_prob), linewidth=0.8)

# TODO make sure all the spillover is actually removed!
# why does assigning the pi vs. uniform for missing counts make no difference??
# TODO the choice of target_p has suspiciously little effect


################################################################################
# try reimplementation of compensate

library(dplyr)
library(ggplot2)
library(magrittr)
library(spatstat.geom)
source("alexander_experiments/from_scratch_compensate.R")
tfm <- function(x) asinh(x/5)
print('testing compensate')
df = read.csv("alexander_experiments/tb_bead.csv")
df_real = read.csv("alexander_experiments/tb_real.csv")
target = 'Yb173Di'
pi_k = 0.9
n_iter = 10
# df = df |> filter(barcode != "Yb176Di")
res = from_scratch_compensate(df, df_real, target, pi_k, n_iter)

# TODO the barcode-colored densities look weird...

ggplot(res) +
    geom_density(aes(tfm(Yb173Di), color=barcode))

ggplot(res |> filter(barcode=="Yb173Di")) +
    geom_density(aes(tfm(Yb173Di)))

ggplot(df_real) +
    geom_density(aes(tfm(Yb173Di)))

# # plot the sanity-check compensation of the beads experiment itself
# ggplot(res) +
#     geom_density(aes(x=tfm(count), weight=compensated_number, color=as.factor(barcode)), alpha=0.5)
