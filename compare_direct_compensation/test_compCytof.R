### Test algorithm using example from Appendix A
setwd('compare_direct_compensation')
library(SingleCellExperiment)
source('compCytof.R')
source('compensate.R')


# define real experiment
nCells <- 6
counts <- matrix(c(3, 5, 17, 3, 17, 2), nrow = 1)
rownames(counts) <- c('protein1')
sce <- SingleCellExperiment(assays = list(counts = counts))
rowData(sce) <- DataFrame(channel_name = c('marker1'))


# TODO: haven't managed to assign bc_key correctly yet!!

# define beads experiment
bc_key <- c(2, 2, 2)
nCells_bead <- 3
counts_bead <- matrix(c(2, 3, 2), nrow = 1)
rownames(counts_bead) <- c("protein2")
sce_bead <- SingleCellExperiment(assays = list(counts = counts_bead))
rowData(sce_bead) <- DataFrame(channel_name = c('marker2'))
assignPrelim(sce, bc_key, assay = 'counts')
sm <- matrix(c(0., 0.05, 0.05, 0.), nrow = 2)
rownames(sm) <- c("marker1", "marker2")
colnames(sm) <- c("marker1", "marker2")
metadata(sce_bead) <- list(spillover_matrix = sm)

# translate markers into barcodes
marker_to_barc <- tibble(marker = c("marker1", "marker2"), barcode = c(1, 2))

# compute spillover compensation
sce_spillr <- compCytof(sce, sce_bead, marker_to_barc, overwrite = FALSE)

compensate



# setwd('compare_direct_compensation')
# library(CATALYST)
# library(tidyverse)
# library(magrittr)
# library(cowplot)
# library(ggplot2)
# source('compensate.R')
# source('compCytof.R')
# source('compDirect.R')
# source('plotDiagnostics.R')
# source('plotDiagnosticsDirect.R')

# # --------- experiment with beads ---------
# bc_key <- c(139, 141:156, 158:176)
# sce_bead_ <- prepData(ss_exp)
# sce_bead_ <- assignPrelim(sce_bead_, bc_key, verbose = FALSE)
# sce_bead_ <- applyCutoffs(estCutoffs(sce_bead_))
# sce_bead_ <- computeSpillmat(sce_bead_)

# # --------- experiment with real cells ---------
# data(mp_cells, package = "CATALYST")
# sce_ <- prepData(mp_cells)