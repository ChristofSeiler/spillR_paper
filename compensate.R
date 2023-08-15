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
  
  # denoise <- function(y, y_min = min(y), y_max = max(y)) {
  #   
  #   # frequency table
  #   tb_obsv <- tibble(y) %>% group_by(y) %>% tally()
  #   y_min_obsv <- min(tb_obsv$y)
  #   y_max_obsv <- max(tb_obsv$y)
  #   tb_pred <- tibble(y = y_min:y_max)
  #   tb_pred %<>% left_join(tb_obsv, by = "y")
  #   
  #   # padding with zero outside of data support
  #   tb_pred %<>% mutate(n = ifelse(y > y_max_obsv | y < y_min_obsv, 0, n))
  #   
  #   # nonparametric density estimate
  #   tb_pred$n[is.na(tb_pred$n)] <- 0
  #   tb_pred %>% mutate(pmf = n/sum(tb_pred$n))
  #   
  # }
  
  # support for target marker
  tb_beads_pmf <- tibble(y = y_min:y_max)
  
  # collect pmf from beads
  for(marker in spillover_markers) {
    
    Fn <- tb_bead %>% 
      dplyr::filter(barcode == marker) %>% 
      pull(all_of(target_marker)) %>%
      ecdf()
    tb <- tb_beads_pmf %>% 
      mutate(pmf = Fn(y) - Fn(y-1)) %>% 
      dplyr::select(y, pmf)
    names(tb) <- c("y", marker)
    tb_beads_pmf %<>% left_join(tb, by = "y")
    
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
    
    # bug: pi <- colSums(post_M)/nrow(post_M)
    # bug fix
    y_pi <- bind_cols(y = tb_pmf$y, post_M)
    y_obsv <- tb_real %>% select(y = all_of(target_marker))
    y_obsv <- left_join(y_obsv, y_pi, by = "y")
    pi <- y_obsv %>% select(-y)
    pi <- colSums(pi)/nrow(pi)

    # # stochastic EM: 
    # # assigns each observation to a class with the highest posterior probability
    # #class <- apply(post_M, 1, function(Mi) sample(colnames(post_M), size = 1, prob = Mi))
    # # categorical EM:
    # # assigns each observation randomly based on posterior probabilities
    # # bug: class <- all_markers[apply(post_M, 1, which.max)]
    # # bug fix
    # class <- all_markers[apply(select(y_obsv, -y), 1, which.max)]
    #   
    # # update pmf from real cells
    # if(sum(class == target_marker) > 0) {
    #   
    #   # bug:
    #   # ys <- tb_pmf[class == target_marker, ] %>% pull(y)
    #   # tb_real_pmf <- tb_real %>% 
    #   #   dplyr::filter(.data[[target_marker]] %in% ys) %>%
    #   #   pull(all_of(target_marker)) %>% 
    #   #   denoise(y_min = y_min, y_max = y_max) %>% 
    #   #   dplyr::select(y, pmf)
    #   # bug fix
    #   Fn <- y_obsv %>% 
    #     mutate(masked = class != target_marker) %>%
    #     filter(!masked) %>% 
    #     pull(y) %>%
    #     ecdf()
    #   tb_real_pmf <- tibble(y = y_min:y_max) %>% mutate(pmf = Fn(y) - Fn(y-1))
    #   
    # } else {
    #   
    #   # if no signal, then use uniform distribution  
    #   tb_real_pmf <- tb_beads_pmf %>% 
    #     dplyr::select(y) %>% 
    #     mutate(pmf = 1/nrow(tb_beads_pmf))
    #   
    # }

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
  tb_spill_prob <- dplyr::select(tb_pmf, y) %>% 
    mutate(spill_prob = if_else(is.na(spill_prob), 0, spill_prob))
  names(tb_spill_prob) <- c(target_marker, "spill_prob")
  
  # compensate
  tb_compensate <- tb_real
  tb_compensate %<>% left_join(tb_spill_prob, by = target_marker)
  tb_compensate %<>% mutate(
    spill = rbinom(n = nrow(tb_compensate), 
                   size = 1, 
                   prob = tb_compensate$spill_prob)
  )
  tb_compensate %<>% mutate(corrected = ifelse(spill == 1, NA, .data[[target_marker]]))
  
  names(tb_compensate)[1] = "uncorrected"
  
  # postprocess spillover probabilities
  tb_spill_prob %<>% mutate(y_tfm = tfm(.data[[target_marker]]))
  fit <- glm(spill_prob ~ y_tfm, family = binomial, data = tb_spill_prob)
  inverse_logit <- function(fit, x) {
    hat <- coef(fit)[1] + coef(fit)[2]*x
    1/(1+exp(-hat))
  }
  tb_spill_prob %<>% mutate(spill_prob_smooth = inverse_logit(fit, y_tfm))
  
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
