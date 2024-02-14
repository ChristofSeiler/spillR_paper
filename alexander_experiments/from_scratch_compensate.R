# implement compensate function from scratch
library(tidyr)
library(dplyr)


compensate = function(df, target, pi_k=0.9, n_iter=3){
    x_min = min(df[[target]])
    x_max = max(df[[target]])
    markers = unique(df$barcode)
    d = length(markers)

    # initialize
    pdfs = list()
    data = data.frame('count'=c(x_min:x_max))
    # TODO I'm not sure I'm estimating the pdfs correctly!
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
        names(pmf_pi) = markers
        pmf_pi = pmf_pi / colSums(pmf_pi, na.rm=TRUE)
        # pmf_pi[is.na(pmf_pi)] = 0
        ka_given_x = pmf_pi / rowSums(pmf_pi, na.rm=TRUE)

        # M-step
        ka = colMeans(ka_given_x, na.rm=TRUE)
        ka = ka / sum(ka)

        # logging
        convergence[i,] = ka

        # # TODO have to be very careful here, not sure this is correct
        # ka_given_x[is.na(ka_given_x)] = 1/d

        # update density estimates
        for(m in markers){
            # TODO check if the weights really make sense
            obs = df[df$barcode==m, target]
            ws = ka_given_x[[m]][obs]
            ws[is.na(ws)] = 0
            ws = ws / sum(ws)
            fit = density(obs, from=x_min, to=x_max, weights=ws)
            pdfs[[m]] = approxfun(fit$x, fit$y)
        }
    }

    print(convergence)
}


if(interactive()){
    print('testing compensate')
    df = read.csv("alexander_experiments/tb_bead.csv")
    target = 'Yb173Di'
    pi_k = 0.9
    n_iter = 3
    compensate(df, target, pi_k, n_iter)
}
