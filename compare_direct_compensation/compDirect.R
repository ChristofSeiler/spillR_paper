library(CATALYST)
library(tidyverse)
library(magrittr)


compDirect <- function(sce, sce_bead) {

    # transformation
    tfm <- function(x) asinh(x / 5)

    # experiment with beads
    counts_bead <- t(assay(sce_bead, "counts"))
    counts_bead <- floor(counts_bead)
    counts_bead <- as_tibble(counts_bead)

    # experiment with real cells
    channel_names <- rowData(sce)[, "channel_name"]
    counts_real <- t(assay(sce, "counts"))
    counts_real <- floor(counts_real)
    colnames(counts_real) <- channel_names
    counts_real <- as_tibble(counts_real)

    # get count values and their probabilities for a given marker
    get_marker_probs <- function(marker, counts){
        return(table(counts[which(colnames(counts)==marker)]) / nrow(counts))
    }

    # get the probability of a given count for each marker in counts_table
    get_count_prob <- function(counts_table, count) {
        if (as.character(count) %in% names(counts_table)) {
            count_prob <- counts_table[match(as.character(count),
                                       names(counts_table))]
        } else {
            count_prob <- 0
        }
        return(count_prob)
    }

    # structure to store compensated counts
    compensated <- counts_real

    # compensate for spillover in all markerssy
    for (c in channel_names) {
        print(paste("compensating spillover in", c))
        # calculate probas for counts of a target marker (real & beads)
        real_probs <- get_marker_probs(c, counts_real)
        bead_probs <- lapply(setdiff(colnames(counts_bead),
                             c(c, "barcode")),
                             get_marker_probs, counts = counts_bead)

        # initialize counts that can be due to spillover
        overlap <- intersect(names(real_probs),
                                 Reduce("c", sapply(bead_probs, names)))

        # iteratively add counts that are more likely in the bead experiment
        done <- FALSE
        while (!done) {
            real_probs_order <- sort(as.vector(real_probs),
                                     index.return = TRUE, decreasing = TRUE)$ix
            for (i in real_probs_order){
                count <- names(real_probs[i])
                if ((count != "0") && (count %in% overlap)) {
                    max_prob <- max(sapply(bead_probs,
                                           get_count_prob,
                                           count = count))
                } else {
                    max_prob <- 0
                }
                if (real_probs[i] < max_prob){
                    print(paste(c, "removing count", count))
                    real_probs <- real_probs[-i]
                    real_probs <- real_probs / sum(real_probs)
                    overlap <- setdiff(overlap, count)
                    compensated[counts_real[, c] == as.integer(count), c] <- NA
                    # break for loop to re-start search
                    break
                }
            }
            if (i == real_probs_order[length(real_probs_order)]){
                done <- TRUE
            }
        }
    }
    # return assay with compensated counts
    assay(sce, "compcounts", FALSE) <- t(compensated)
    assay(sce, "compexprs", FALSE) <- t(tfm(compensated))
    return(sce)
}
