#' @title rra
#'
#' @description The main implementation of reduced-rank approach
#'              by Dashan Huang et al. (MS, 2022)
#'
#' @param R The t*(1+n) dataframe containing the return data with date.
#'          Make sure that the date column is the first column.
#' @param G The t*(1+l) dataframe containing the factor proxies with date.
#'          Make sure that the date column is the first column.
#' @param K The number of factors you want to extract.
#' @param price_err Bool flag of including price error
#' @param alpha The annual pricing error e.g. alpha=1 if 1% per year
#'
#' @return The K extracted factors factors($Phi*G$)
#' @importFrom POET "POET"
#' @export rra


rra <- function(
    R, G, K,
    price_err=FALSE, annual_error=1,
    weight=FALSE
) {
    if (nrow(R)!=nrow(G)){
        stop("Time periods do not match between two dataframe.")
    }

    R_m <- as.matrix(R[, 2:ncol(R)])
    G_m <- as.matrix(G[, 2:ncol(G)])

    if(price_err){
        std_err <- sapply(data.frame(R_m), sd)
        alpha <- std_err*(annual_error/12)
        rra <- rra_pere(R, G, K, alpha)
        return(rra)
    }

    if(weight){
        S2 <- POET(t(cbind(rep(1, nrow(G_m)),G_m)),K,.5,'soft','vad')$SigmaY
        rra <- rra_w(R, G, K, S2)
        return(rra)
    }

    rra <- rra_main(R, G, K)
    return(rra)
}
