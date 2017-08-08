#' @importFrom stats pnorm 
#' @importFrom stats rbinom 
#' @importFrom stats rmultinom 
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom stats dnorm 
#' @importFrom stats dpois 
#' @importFrom stats rpois 
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom abind abind
NULL

#' Bernoulli likelihood
#'
#' Returns the log-likelihood for a vector of parameter values, binary responses and design matrix
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of count responses, length equals number of patients
#' @param x design matrix
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates, 
#' \code{z}: covariate matrix
#' @return log-likelihood value
#'
#' @export

llik_bernoulli = function(theta, y, x, ...) {
  more = list(...)
  if (!is.null(more$beta)){
  beta = more$beta
  z = more$beta
  add = t(z)%*%beta
  } else add = 0
  if (!all(y %in% c(0,1))) stop('non-binary response, y takes 0 or 1 values')
  c = c(t(x)%*%theta + add)
  out = c(mapply(function(y, c) log(y * (1 - pnorm(0, c, 1)) +
                                      (1 - y) * pnorm(0, c, 1)), y, c))
  out = sum(out)
  if (out == -Inf) out = -1e6
  return(out)
}

#' Poisson likelihood
#'
#' Returns the log-likelihood for a vector of parameter values, count responses and design matrix
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of count responses, length equals number of patients
#' @param x design matrix
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates, 
#' \code{z}: covariate matrix
#' @return log-likelihood value
#'
#' @export

llik_poisson = function(theta, y, x, ...) {
  more = list(...)
  if (!is.null(more$beta)){
  beta = more$beta
  z = more$beta
  add = t(z)%*%beta
  } else add = 0
  if (!all(y >= 0)) stop('negative response, y takes non-negative')
  if (all(y %in% c(0,1))) warning('response values are 0 or 1; response.type = binary may be a better option')
  c = exp(c(t(x)%*%theta + add))
  out = c(mapply(function(y, c) dpois(y, c, log = T), y, c))
  out = sum(out)
  if (out == -Inf) out = -1e6
  return(out)
}

#' Prior
#'
#' Returns log-prior for a vector of parameter values - this is the default normal prior function in case the user does not want to
#' bother to decide on a prior distribution
#'
#' @param par vector of parameter values
#' @param ... the normal prior hyper-parameters (optional), can be specified as vectors of the same size as \code{par} or
#' scalars. The default is \code{pmean = 0} and \code{psd = 10}.
#' @return log-prior value
#'
#' @export

lp_default = function(par, ...) {
  args = list(...)
  pmean = if(is.null(args$pmean)) 0 else args$pmean
  psd = if(is.null(args$psd)) 10 else args$psd
  out = if(is.null(par)) 0 else sum(mapply(dnorm, x = par, mean = pmean, sd = psd, MoreArgs = list(log = T)))
  return(out)
}

#' Posterior
#'
#' Returns log-posterior for a vector of parameter values, count responses and design matrix
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of responses, length equals number of patients
#' @param x design matrix
#' @param response.type response type, either 'binary' or 'count'
#' @param lp_theta a function returning log-prior for effect sizes; the default is lp_default which is a diffuse normal prior,
#' but can be replaced with a user-specified functin
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates, 
#' \code{z}: covariate matrix, \code{lp_beta}: a function returning log-prior for coefficients; the default is lp_default which is a diffuse normal prior,
#' but can be replaced with a user-specified functin
#' @return log-posterior value
#'
#' @export

lpost = function(theta, y, x, response.type, lp_theta = lp_default, ...) {
  more = list(...)
  if (!is.null(more$beta)) {
    beta = more$beta
    z = more$z
    lp_beta = more$lp_beta
    bprior = lp_beta(beta)
  } else bprior = 0
  if (response.type == 'binary') ll = llik_bernoulli(theta, y, x, ...)
  if (response.type == 'count') ll = llik_poisson(theta, y, x, ...)
  lp = lp_theta(theta) + bprior
  return(ll + lp)
}

#' MCMC step
#'
#' This function perfomrs the nested MCMC steps required in SMC
#'
#' @param theta vector of effect sizes, length equals number of arms
#' @param y vector of responses, length equals number of patients
#' @param x design matrix
#' @param qt proposal standard deviation
#' @param response.type response type, either 'binary' or 'count'
#' @param ng number of MCMC steps for each particle. It can be as small as 1, but need to be larger if batch updates are made
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates, 
#' \code{z}: covariate matrix, \code{lp_beta}: a function returning log-prior for coefficients; the default is lp_default which is a diffuse normal prior,
#' but can be replaced with a user-specified functin, \code{qb}: proposal variance for \code{beta}
#' @return vector of updated parameters
#'
#' @export

Gibbs = function(theta, y, x, qt, response.type, ng, ...) {
  more = list(...)
  nt = length(qt)
  new = theta + qt*rnorm(nt)
  for (j  in 1:ng) {
    for (i in 1:nt) {
      new_th = theta
      new_th[i] = new[i]
      lp = lpost(theta, y, x, response.type = response.type, ...)
      lpnew = lpost(new_th, y, x, response.type = response.type, ...)
      p = min(1, exp(lpnew - lp))
      if (runif(1) <= p) theta = new_th
    }
    if (!is.null(more$beta)) {
      beta = more$beta
      z = more$z
      qb = more$qb
      if (!is.null(more$lp_beta)) lp_beta = more$lp_beta else lp_beta = lp_default
      nb = length(qb)
      newb = beta + qb*rnorm(nb)
      for (i in 1:nb) {
        new_b = beta
        new_b[i] = newb[i]
        lp = lpost(theta, y, x, response.type = response.type, beta = beta, z = z, lp_beta = lp_beta)
        lpnew = lpost(theta, y, x, response.type = response.type, beta = newb, z = z, lp_beta = lp_beta)
        p = min(1, exp(lpnew - lp))
        if (runif(1) <= p) beta = new_b
      }
    } else beta = NULL
  }
  out = list(theta = theta, beta = beta)
  return(out)
}

#' SMC update
#'
#' This function updates a posterior sample of effect sizes according to newly observed data
#'
#' @param theta matrix of posterior samples for effect sizes, columns: arms, rows: samples
#' @param y vector of recent responses
#' @param x design matrix for the recent observations
#' @param N number of particles
#' @param response.type response type, either 'binary' or 'count'
#' @param ... covariate related arguments, \code{beta}: vector of coefficients, length equals number of covariates, 
#' \code{z}: covariate matrix
#' @return updated posterior samples
#'
#' @export

SMC_step = function(theta, y, x, N, response.type = 'binary', ...) {
  more = list(...)
  smA = split(theta, row(theta))
  if (!is.null(more$beta)) {
    beta = more$beta
    z = more$z
    smB = split(beta, row(beta)) 
    if (response.type == 'binary') ll = mapply(llik_bernoulli, smA, beta = smB, MoreArgs = list(x = x, y = y, z = z))
    if (response.type == 'count') ll = mapply(llik_poisson, smA, beta = smB,  MoreArgs = list(x = x, y = y, z = z))
    }
  else {
    if (response.type == 'binary') ll = apply(theta, 1, llik_bernoulli, x = x, y = y)
    if (response.type == 'count') ll = apply(theta, 1, llik_poisson, x = x, y = y)
    beta = NULL
    }
  w = exp(ll)/sum(exp(ll))
  index = sample(1:N, N, prob = w, replace = T)
  theta = theta[index,]
  out = list(theta = theta, beta = beta)
  return(out)
}

#' Probabilities of superiority
#'
#' This function estimates the probability of seperiority for each treatment arm from a posterior sample of effect sizes
#'
#' @param theta matrix of posterior samples for effect sizes, columns: arms, rows: samples
#' @return vector of probabilities of superiority
#'
#' @export

sup_check = function(theta) {
  n = length(theta)
  check = c()
  for (i in 1:n) check[i] = (theta[i] == max(theta))
  return(check)
}

#' RAR simulator
#'
#' This function simulates a RAR trial with SMC updates
#'
#' @param nt Number of treatment arms
#' @param theta0 vector of effect sizes
#' @param nb batch size for each update
#' @param maxN maximum number of patients that the trial can run for; the trial is stopped when \code{maxN} is achieved
#' @param N Number of particles
#' @param upper upper threshold for probability of superiority; the trial is stopped when one of the trial arms reaches
#' this value
#' @param lower lower threshold for probability of superiority; an arm is stopped when its corresponding probability of
#' superiority falls below this value
#' @param burn the size of pre-adaptation sample; the default is 10 times number of treatments
#' @param response.type response type, either 'binary' or 'count'
#' @return updated posterior samples
#'
#' @export

RAR_sim = function(nt, theta0, nb = 1, maxN = 500, N = 1000, upper = 0.95, lower = .05,
                   burn = 10*nt, response.type = 'binary') {
  ng = nb
  j = 0
  x = array(0, dim = c(nt, 1))
  y = NULL
  theta = array(rnorm(N*nt, 0, 10), dim = c(N, nt, 1))
  psup = array(rep(1/nt, nt), dim = c(nt, 1))
  repeat {
    j = j + 1
    xb = rmultinom(nb, 1, prob = sqrt(psup[,j]))
    if (response.type == 'binary') {
      yb = apply(t(xb)%*%theta0, 1, function(z) rbinom(1, 1, prob = 1/(1+exp(-z))))
    }
    if (response.type == 'count') {
      yb = apply(t(xb)%*%theta0, 1, function(z) rpois(1, exp(z)))
    }
    x = abind(x, xb, along = 2)
    y = c(y, yb)
    smc_out = SMC_step(theta[,,j], x = xb, y = yb, N = N, response.type)
    theta_new = smc_out$theta
    qt = apply(theta_new, 2, sd)/3
    Gibbs_out = apply(theta_new, 1, Gibbs, y = y, x = x[,-1], qt = qt,
                        response.type = response.type, ng)
    theta_new = t(sapply(Gibbs_out, function(l) return(l$theta)))
    theta = abind(theta, theta_new, along = 3)
    check = t(apply(theta[,,j+1], 1, sup_check))
    if (length(y) < burn) {
      psup = abind(psup, psup[,j], along = 2)
    } else psup = abind(psup, apply(check, 2, mean), along = 2)
    psup[psup[,j+1] < lower | psup[,j] ==0,j+1] = 0
    if (max(psup[,j+1]) > upper | length(y) > maxN) break
  }
  out = list(psup = psup, theta = theta, y = y, x = x[,-1])
  class(out) = 'trial'
  return(out)
}





