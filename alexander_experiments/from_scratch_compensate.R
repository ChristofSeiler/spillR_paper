# implement compensate function from scratch
library(tidyr)
library(dplyr)


compensate = function(df, target, pi_k=0.9){
    min_x = min(df$X)
    max_x = max(df$X)
    d = ncol(df)
    markers = unique(df$barcode)

    # data
    pdfs = list()
    data = data.frame('count'=c(min_x:max_x))
    for(i in markers){
        # get frequency of each count for barcode
        bc_counts = df |>
                    filter(barcode==i) |>
                    transmute('count'=.data[["X"]], !!i:=.data[[target]])
        data = left_join(data, bc_counts, by='count')
        # compute ecdf for barcode
        fit = density(rep(bc_counts$count, times=bc_counts[[i]]))
        pdfs = append(pdfs, approxfun(fit$x, fit$y))
    }
    names(pdfs) = markers
    pi = c(pi_k, vapply(2:d, function(x) (1-pi_k)/(d-1), numeric(1)))

    # initialize
    pmf_data = lapply(markers, function(x) pdfs[[x]](data[[x]])) |> bind_cols()
    pmf_data[is.na(pmf_data)] = 0
    names(pmf_data) = markers

    # TODO we only ever seem to have one marker for each count?? not clear why...

    # E-step


    return(0)
}