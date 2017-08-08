#' Probability of superiority plot
#'
#' @param trial An object of class 'trial'
#'
#' @return Probabilities of superiority for each arm evolving through the course of the trial
#' @export


psup_plot = function(trial) {
  psup = trial$psup
  df = melt(psup)
  names(df) = c('treatment', 'interim_look', 'p.best')
  df$treatment = as.factor(df$treatment)
  p = ggplot(df, aes(x = interim_look, y = p.best, color = treatment)) +
    geom_line() + geom_hline(yintercept = 0.95, color = 'darkgrey') +
    scale_color_brewer(palette = 'Set1')
  return(p)
}

#' Posterior density plot
#'
#' @param trial An object of class 'trial'
#'
#' @return 5 snapshots of the posterior density plots for the effect sizes through teh trial
#' @export

post_plot = function(trial) {
  theta = trial$theta
  np = dim(theta)[3]
  select = floor(seq(1, np, length = 5))
  df = melt(theta)[,-1]
  names(df) = c('treatment', 'interim_look', 'theta')
  df0 = df[df$interim_look %in% select,]
  df0$interim_look = as.factor(df0$interim_look)
  #df1 = data.frame(theta0 = theta0, treatment = as.factor(sort(unique(df0$treatment))))
  p = ggplot(df0, aes(x = theta, fill = interim_look)) +
    geom_density(alpha = .5, color = 'grey') +
    #geom_vline(data = df1, aes(xintercept = theta0)) +
    facet_grid(treatment ~ ., labeller = label_both) +
    scale_fill_brewer() + xlim(-3, 3)
  return(p)
}

#' Data plot
#'
#' @param trial An object of class 'trial'
#'
#' @return A visual summary of the trial data
#' @export

data_plot = function(trial) {
  y = trial$y
  x = trial$x
  nt = dim(trial$theta)[2]
  treat = 1:nt
  df = data.frame(cbind(y, t(x)%*%treat, 1:length(y)))
  names(df) = c('response', 'treatment', 'patient')
  df$treatment = as.factor(df$treatment)
  if (all(y %in% c(0,1))) df$response = as.factor(df$response)
  p = ggplot(df, aes(x = patient, y = response, color = response)) +
    geom_point(size = 3) + facet_grid(treatment ~ ., labeller = label_both)
  return(p)
}



 
