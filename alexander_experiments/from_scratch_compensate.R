# implement compensate function from scratch
if(interactive()){rm(list=ls())}
library(tidyr)
library(dplyr)


from_scratch_compensate = function(df, df_real, target, pi_k=0.9, n_iter=3){
    # choose max range of support
    x_min = min(min(df[[target]]), min(df_real[[target]]))
    x_max = max(max(df[[target]]), max(df_real[[target]]))
    markers = unique(df$barcode)
    d = length(markers)

    # initialize
    pdfs = list()
    data = data.frame('count'=c(x_min:x_max))
    for(i in markers){
        # get frequency of each count for barcode
        bc_counts = df[df$barcode==i, target]
        fit = density(bc_counts, from=x_min, to=x_max)
        pdfs = append(pdfs, approxfun(fit$x, fit$y))
        # create table of counts for each barcode
        bc_table = data.frame(table(bc_counts))
        names(bc_table) = c('count', i)
        bc_table$count = as.numeric(as.character(bc_table$count))
        data = left_join(data, bc_table, by='count')
    }
    names(pdfs) = markers
    pi = c(pi_k, vapply(2:d, function(x) (1-pi_k)/(d-1), numeric(1)))
    convergence <- matrix(nrow=n_iter, ncol=d)
    convergence[1,] = pi
    
    # EM algorithm
    for(i in 2:n_iter){
        # E-step
        pmf_pi = bind_cols(lapply(markers, function(x) pdfs[[x]](data[[x]])))
        pmf_pi = pmf_pi * pi
        pmf_pi = pmf_pi / colSums(pmf_pi, na.rm=TRUE)
        pmf_pi = pmf_pi / rowSums(pmf_pi, na.rm=TRUE)
        names(pmf_pi) = markers

        # M-step: update pi
        pi = colMeans(pmf_pi, na.rm=TRUE)
        pi = pi / sum(pi)

        # logging
        convergence[i,] = pi

        # # TODO have to be very careful here, not sure this is correct
        # pmf_pi[is.na(pmf_pi)] = 1/d

        # update density estimates
        for(m in markers){
            # TODO check if the weights really make sense
            obs = df[df$barcode==m, target]
            ws = pmf_pi[obs+1, m]  # need +1 because R is 1-indexed
            ws[is.na(ws)] = 0
            ws = ws / sum(ws)
            fit = density(obs, from=x_min, to=x_max, weights=ws)
            pdfs[[m]] = approxfun(fit$x, fit$y)
        }
    }
    print(convergence)

    # # perform correction on the bead data itself - just as a sanity check
    # # (basically perform the E-step once more)
    # pmf_pi = bind_cols(lapply(markers, function(x) pdfs[[x]](data[[x]])))
    # pmf_pi = pmf_pi / colSums(pmf_pi, na.rm=TRUE)
    # names(pmf_pi) = markers
    # pmf_pi = pmf_pi / rowSums(pmf_pi, na.rm=TRUE)
    # pmf_pi[is.na(pmf_pi)] = 0
    # compensated = pmf_pi * rowSums(data[markers], na.rm=TRUE)
    # compensated$count = data$count
    # # get back into original format!
    # res = data.frame('count'=c(), 'compensated_number'=c(), 'barcode'=c())
    # for(m in markers){
    #     compensated_m = compensated[compensated[,m]>0, c('count', m)]
    #     names(compensated_m) = c('count', 'compensated_number')
    #     compensated_m$barcode = m
    #     res = bind_rows(res, compensated_m)
    # }

    # Compensate spillover on real data (basically E-step)
    r_pmf_pi = bind_cols(lapply(markers, 
                                function(x) pdfs[[x]](df_real[[target]])))
    r_pmf_pi = r_pmf_pi * pi
    r_pmf_pi = r_pmf_pi / colSums(r_pmf_pi, na.rm=TRUE)
    r_pmf_pi / rowSums(r_pmf_pi, na.rm=TRUE)
    names(r_pmf_pi) = markers
    # find the index of the largest value of each row
    df_real$barcode = markers[apply(r_pmf_pi, 1, function(x)which.max(x))]
    
    # TODO: instead of removing counts, we could multiply the target marker probability with the actual count? (or just set the counts which we believe to be spillover to zero?)

    return(df_real)
}


if(interactive()){
    print('testing from_scratch_compensate')
    df = read.csv("alexander_experiments/tb_bead.csv")
    df_real = read.csv("alexander_experiments/tb_real.csv")
    target = 'Yb173Di'
    pi_k = 0.9
    n_iter = 10
    from_scratch_compensate(df, df_real, target, pi_k, n_iter)
}
