#' Compute spillover probability and correct for spillover
#'
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import CATALYST
#' @import ggplot2
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param ch Character string specifying the channel to plot
#'
#' @return A list of \code{\link[ggplot2]{ggplot2}} plots
#'
#' @examples
#' sce <- prepData(mp_cells)    # real cells
#' sce_bead <- prepData(ss_exp) # beads
#' sce <- spillR::compCytof(sce, sce_bead, overwrite = FALSE)
#' plotDiagnostics(sce, "Yb173Di")
plotDiagnosticsDirect <- function(sce, ch) {
  
  tfm <- function(x) asinh(x/5) # TODO: this should not be hardcoded
  
  # before and after correction
  exprs <- sce %>% 
    assay("exprs") %>%
    t() %>%
    as_tibble()
  names(exprs) <- rowData(sce)$channel_name
  
  compexprs <- sce %>% 
    assay("compexprs") %>%
    t() %>%
    as_tibble()
  names(compexprs) <- rowData(sce)$channel_name
  
  before_after <- tibble(
    cell = 1:nrow(exprs),
    before = exprs %>% pull(ch),
    after = compexprs %>% pull(ch)
  )
  before_after %>% 
    pivot_longer(-cell, names_to = "correction") %>% 
    mutate(correction = factor(correction, levels = c('after', 'before'))) %>%
    ggplot(aes(value, color = correction, linetype = correction)) + 
    geom_freqpoly(alpha = 1.0, bins = 50, linewidth = 0.8) +
    scale_color_manual(values = c('#00BFC4', '#F8766D')) +
    scale_linetype_manual(values = c('solid', 'dashed')) +
    xlab(paste0("tfm(", ch, ")")) +
    ggtitle("Spillover Compensation on Real Cells")
  
}
