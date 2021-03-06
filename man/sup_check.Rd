% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAR_sim.R
\name{sup_check}
\alias{sup_check}
\title{Probabilities of superiority}
\usage{
sup_check(theta)
}
\arguments{
\item{theta}{matrix of posterior samples for effect sizes, columns: arms, rows: samples}
}
\value{
vector of probabilities of superiority
}
\description{
This function estimates the probability of seperiority for each treatment arm from a posterior sample of effect sizes
}
