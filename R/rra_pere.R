#' @title rra_pere
#'
#' @description he implementation of reduced-rank approach with pricing error
#' by Dashan Huang et al. (MS, 2022)
#'
#' @param R The t*n matrix containing the return data.
#' @param G The t*l matrix containing the factor proxies.
#' @param K The number of factors you want to extract.
#' @param alpha The matrix of alpha*std_error of each basis portfolio.
#'
#' @return The K extracted factors($Phi*G$)

rra_pere <- function(R, G, K, alpha) {
    # using R=a+BG+u to extract K factors
    # R: t*n, G: t*l
    # output: T*K factors with RRA

    t = nrow(G); # the number of months
    l = ncol(G); # the number of latent factors
    n = ncol(R); # the number of excess returns of asset
    z = cbind(matrix(data=1, nrow=t, ncol=1), G)
    m = l+1;
    R = R - kronecker(matrix(1,t,1), t(alpha))

    x = matrix(data=1, nrow=t, ncol=1)

    w1 = diag(n) # weighting matrix 1
    w2 = diag(m) # weighting matrix 2
    p0 = z%*%w2%*%t(z)
    p = p0 # - p0%*%x%*%solve(t(x)%*%p0%*%x)%*%t(x)%*%p0
    q = t(G)%*%p%*%G
    sqrt_q = sqrtm(solve(q))
    a_part = t(G)%*%p%*%R/(t^2)
    a = t(sqrt_q)%*%(a_part)%*%w1%*%t(a_part)%*%sqrt_q
    E = eigen(a, TRUE)$vectors
    v = eigen(a, TRUE)$values
    phi = sqrt_q %*% E[, 1:K]
    Gstar = G %*% phi
    beta = solve(t(Gstar)%*%p%*%Gstar) %*% t(Gstar)%*%p%*%R
    theta = phi %*% beta
    alpha = solve(t(x)%*%p0%*%x)%*%t(x)%*%p0%*%(R-G%*%theta)
    alpha = t(alpha)
    u = R - x%*%t(alpha) - G%*%theta
    s1 = diag(diag(t(u)%*%u/t))
    s2 = diag(diag(t(z)%*%z/t))

    w1 = solve(s1)
    w2 = solve(s2)
    p0 = z%*%w2%*%t(z)
    p = p0 # - p0%*%x%*%solve(t(x)%*%p0%*%x)%*%t(x)%*%p0
    q = (t(G)%*%p%*%G)/t^2
    sqrt_q = sqrtm(solve(q))
    a_part = t(G)%*%p%*%R/(t^2)
    a = t(sqrt_q)%*%(a_part)%*%w1%*%t(a_part)%*%sqrt_q
    E = eigen(a, TRUE)$vectors
    v = eigen(a, TRUE)$values
    phi = sqrt_q %*% E[, 1:K]
    Gstar = G %*% phi
    Gstar = Re(Gstar)
    return(Gstar)

}
