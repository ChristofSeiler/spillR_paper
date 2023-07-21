setwd('compare_direct_compensation')
library(CATALYST)
library(tidyverse)
library(magrittr)
library(cowplot)
library(ggplot2)
source('compensate.R')
source('compCytof.R')
source('compDirect.R')

# --------- experiment with beads ---------
bc_key <- c(139, 141:156, 158:176)
sce_bead <- prepData(ss_exp)
sce_bead <- assignPrelim(sce_bead, bc_key, verbose = FALSE)
sce_bead <- applyCutoffs(estCutoffs(sce_bead))
sce_bead <- computeSpillmat(sce_bead)

# --------- experiment with real cells ---------
data(mp_cells, package = "CATALYST")
sce <- prepData(mp_cells)

# --------- table for mapping markers and barcode ---------
marker_to_barc <- 
  rowData(sce_bead)[,c("channel_name", "is_bc")] %>%
  as_tibble %>%
  dplyr::filter(is_bc == TRUE) %>%
  mutate(barcode = bc_key) %>%
  select(marker = channel_name, barcode)

# --------- call compensate from compCytof package ---------
sce_spillr <- compCytof(sce, sce_bead, marker_to_barc, overwrite = FALSE)
sce_direct <- compDirect(sce, sce_bead)


# --------- 2d histogram from spillR and direct ---------

ps <- list()
for (chs in list(c("Yb173Di", "Yb174Di"), c("Yb172Di", "Yb173Di"), 
                 c("Yb172Di", "Yb171Di"))) {
    p_orig <- plotScatter(sce_spillr, chs, assay = "exprs") +
        ylim(0, NA) + xlim(0, NA) +
        labs(title = "Uncorrected", x = chs[1], y = chs[2]) + 
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    p_spillR <- plotScatter(sce_spillr, chs, assay = "compexprs") +
        ylim(0, NA) + xlim(0, NA) +
        labs(title = "spillR", x = chs[1], y = chs[2]) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    p_direct <- plotScatter(sce_direct, chs, assay = 'compexprs') +
        ylim(0, NA) + xlim(0, NA) +
        labs(title = "compDirect", x = chs[1], y = chs[2]) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ps <- append(ps, list(p_orig, p_spillR, p_direct))
}

grid <- plot_grid(plotlist = ps, ncol = 3)
grid
ggsave("comparison.pdf", grid, device = "pdf")

# TODO need to make sure we really have the same assays in our return object!