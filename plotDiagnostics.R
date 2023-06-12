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
plotDiagnostics <- function(sce, ch) {
  
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
  p_before_after <- before_after %>% 
    pivot_longer(-cell, names_to = "correction") %>% 
    mutate(correction = factor(correction, levels=c('before', 'after'))) %>%
    ggplot(aes(value, color = correction, linetype = correction)) + 
    geom_freqpoly(alpha = 1.0, bins = 50) +
    scale_color_manual(values = c('#00BFC4', '#F8766D')) +
    scale_linetype_manual(values = c('solid', 'dashed')) +
    xlab(paste0("tfm(", ch, ")"))
  
  # diagnostic plot for spillover estimate
  tb_bead <- metadata(sce)$beads_distr[[ch]]
  tb_spill_prob <- metadata(sce)$spillover_est[[ch]]
  p_spill <- tb_bead %>%
    ggplot(aes(tfm(.data[[ch]]), color = barcode)) +
    geom_density(adjust = 1) +
    geom_line(data = tb_spill_prob, 
              aes(tfm(.data[[ch]]), spill_prob_smooth), 
              color = "black", 
              linetype = "longdash") + 
    ylab("density")

  list(
    p_before_after = p_before_after, 
    p_spill = p_spill
    )
  
}
