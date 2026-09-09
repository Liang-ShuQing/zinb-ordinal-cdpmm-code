\name{fit_zinb_ordinal}
\alias{fit_zinb_ordinal}
\title{Fit joint ZINB-ordinal model with CDPMM random effects}
\description{
Thin wrapper around \code{fit_joint_model}.
}
\usage{
fit_zinb_ordinal(data, chain = 2000, burn = 1000, thin = 5, G = 8, ...)
}
\arguments{
\item{data}{A list / \code{zinbcdpmm_data} from \code{simulate_zinb_ordinal}
or with the same fields as \code{generate_data} output.}
\item{chain}{MCMC iterations.}
\item{burn}{Burn-in.}
\item{thin}{Thinning interval.}
\item{G}{Truncation level of the stick-breaking CDPMM (default 8).}
\item{...}{Additional arguments passed to \code{fit_joint_model}.}
}
\value{
An object of class \code{zinbcdpmm_fit}.
}

