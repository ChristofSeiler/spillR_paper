#' Compute spillover probability and correct for spillover
#'
#' @import magrittr
#' @import dplyr
#' @import spatstat.geom
#'
#' @param tb_real Data frame or tibble with proteins counts of real experiment
#' @param tb_bead Data frame or tibble with proteins counts of bead experiment
#' @param target_marker Marker name in real experiment
#' @param spillover_markers Marker names in bead experiment
#'
#' @return A list of class \code{spillr} containing
#'   \item{tb_compensate}{corrected real cells}
#'   \item{tb_spill_prob}{probability curve}
#'   \item{convergence}{covergence table of EM algorithm}
#'   \item{tb_real}{input real cells}
#'   \item{tb_bead}{input bead cells}
#'   \item{target_marker}{input marker in real experiment}
#'   \item{spillover_markers}{input markers in bead experiment}
#'
#' @examples
#' set.seed(23)
#' tb_real <- generate_real()
#' tb_bead <- generate_bead()
#' target_marker <- "A"
#' spillover_marker <- "B"
#' spillR::compensate(tb_real, tb_bead, "A", "B")
compensate <- function(tb_real, tb_bead, target_marker, spillover_markers) {
  
  tfm <- function(x) asinh(x/5)
  all_markers <- c(target_marker, spillover_markers)
  y_target <- tb_real %>% pull(.data[[target_marker]])
  y_min <- min(y_target)
  y_max <- max(y_target)
  #tb_bead %<>% filter(
  #  y_min <= .data[[target_marker]] & .data[[target_marker]] <= y_max
  #  )
  
  # support for target marker
  tb_beads_pmf <- tibble(y = y_min:y_max)
  
  # collect pmf from beads
  for(marker in spillover_markers) {
     n <- nrow(dplyr::filter(tb_bead, barcode == marker))
     if(n > 0) {
       Fn <- tb_bead %>% 
         dplyr::filter(barcode == marker) %>% 
         pull(all_of(target_marker)) %>%
         ecdf()
       tb <- tb_beads_pmf %>% 
         mutate(pmf = Fn(y) - Fn(y-1)) %>% 
         dplyr::select(y, pmf)
       names(tb) <- c("y", marker)
       tb_beads_pmf %<>% left_join(tb, by = "y")
     } else {
       tb <- tibble(y = y_min:y_max, pmf = 1/nrow(tb_beads_pmf))
       names(tb) <- c("y", marker)
       tb_beads_pmf %<>% left_join(tb, by = "y")
     }
  }
  
  # --------- step 1: initialize ---------
  
  # prior probability
  pi <- rep(1, length(all_markers))
  pi <- pi/length(pi)
  names(pi) <- all_markers
  
  # add pmf from real cells
  Fn <- tb_real %>% 
    pull(all_of(target_marker)) %>% 
    ecdf()
  tb_real_pmf <- tibble(y = y_min:y_max) %>% 
    mutate(pmf = Fn(y) - Fn(y-1)) %>%
    dplyr::select(y, pmf)
  names(tb_real_pmf) <- c("y", target_marker)
  
  # join beads and real
  tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
  
  # --------- step 2: iterate ---------
  
  n_iter <- 10
  convergence <- matrix(nrow = n_iter, ncol = length(pi)+1)
  colnames(convergence) <- c("iteration", names(pi))
  
  for(i in 1:n_iter) {
    # remove counts without any signal
    total <- rowSums(dplyr::select(tb_pmf, all_of(all_markers)))
    tb_pmf %<>% mutate(total = total)
    tb_pmf %<>% dplyr::filter(total > 0)
    tb_pmf %<>% dplyr::select(-total)
    
    # E-step
    
    # membership probabilities
    M <- tb_pmf %>% dplyr::select(all_of(all_markers)) %>% as.matrix()
    PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
    post_M <- PI * M
    post_M <- post_M/rowSums(post_M)
    
    # M-step
    
    # update prior probability
    y_pi <- bind_cols(y = tb_pmf$y, post_M)
    y_obsv <- tb_real %>% select(y = all_of(target_marker))
    y_obsv <- left_join(y_obsv, y_pi, by = "y")
    pi <- y_obsv %>% select(-y)
    pi <- colSums(pi)/nrow(pi)

    # new weighted empirical density estimate
    Fn <- ewcdf(pull(y_obsv, y), weights = pull(y_obsv, all_of(target_marker)))
    tb_real_pmf <- tibble(y = y_min:y_max) %>% mutate(pmf = Fn(y) - Fn(y-1))
    names(tb_real_pmf) <- c("y", target_marker)
    
    # update join
    tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
    
    # keep track
    convergence[i,] <- c(i, pi)
  }
  
  # --------- spillover probability curve ---------
  
  # calculate posterior spillover probability for each cell
  M <- tb_pmf %>% dplyr::select(all_of(all_markers)) %>% as.matrix()
  PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
  post_M <- PI * M
  post_M <- post_M/rowSums(post_M)
  spill_prob <- 1-post_M[,target_marker]
  tb_spill_prob <- tibble(y = y_min:y_max, spill_prob)
  names(tb_spill_prob) <- c(target_marker, "spill_prob")
  
  # count number of cells per y for weights
  tb_tally <- left_join(
    tb_real %>% group_by(.data[[target_marker]]) %>% tally(name = "cells"),
    tb_bead %>% group_by(.data[[target_marker]]) %>% tally(name = "beads"),
    by = target_marker
  ) 
  tb_tally[is.na(tb_tally)] <- 0
  tb_tally %<>% mutate(n = cells + beads) %>% select(target_marker, n)
  
  # post-process spillover probabilities
  tb_spill_prob %<>% mutate(y_tfm = tfm(.data[[target_marker]]))
  tb_spill_prob %<>% left_join(tb_tally, by = target_marker)
  
  # binomial smoothing
  # fit <- glm(spill_prob ~ y_tfm, family = binomial, weights = n,
  #            data = tb_spill_prob)
  # inverse_logit <- function(fit, x) {
  #   hat <- coef(fit)[1] + coef(fit)[2]*x
  #   1/(1+exp(-hat))
  # }
  # tb_spill_prob %<>% mutate(spill_prob_smooth = inverse_logit(fit, y_tfm))
  
  # moving average
  width <- 5
  window <- cbind(
    sapply(rev(seq(width)), function(i) lag(tb_spill_prob$spill_prob, n = i)),
    tb_spill_prob$spill_prob,
    sapply(seq(width), function(i) lead(tb_spill_prob$spill_prob, n = i))
  )
  weight <- cbind(
    sapply(rev(seq(width)), function(i) lag(tb_spill_prob$n, n = i)),
    tb_spill_prob$n,
    sapply(seq(width), function(i) lead(tb_spill_prob$n, n = i))
  )
  tb_spill_prob$spill_prob_smooth <-
    rowSums(window*weight, na.rm = T)/rowSums(weight, na.rm = T)
  
  # compensate
  tb_compensate <- tb_real
  tb_compensate %<>% left_join(tb_spill_prob, by = target_marker)
  tb_compensate %<>% mutate(
    spill = rbinom(n = nrow(tb_compensate), 
                   size = 1, 
                   prob = tb_compensate$spill_prob_smooth)
  )
  tb_compensate %<>% 
    mutate(corrected = ifelse(spill == 1, NA, .data[[target_marker]]))
  names(tb_compensate)[1] = "uncorrected"
  
  # return spillr object
  res <- NULL
  res$tb_compensate <- tb_compensate
  res$tb_spill_prob <- tb_spill_prob
  res$convergence <- convergence
  res$tb_real <- tb_real
  res$tb_bead <- tb_bead
  res$target_marker <- target_marker
  res$spillover_markers <- spillover_markers
  class(res) <- "spillr"
  res
  
}
