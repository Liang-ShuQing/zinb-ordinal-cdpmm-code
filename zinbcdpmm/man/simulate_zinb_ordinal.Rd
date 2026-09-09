\name{simulate_zinb_ordinal}
\alias{simulate_zinb_ordinal}
\title{Simulate longitudinal ZINB + ordinal data}
\description{
Thin wrapper around \code{generate_data}. Returns an object of class
\code{zinbcdpmm_data}.
}
\usage{
simulate_zinb_ordinal(
  n = 100,
  nis = 10,
  seed = 2026,
  random_nis = FALSE,
  nis_range = 1:20,
  scenario = 1,
  C = 5,
  re_dist = "normal",
  mvt_df = 4,
  rho = NULL,
  zi_intercept = NULL,
  ...
)
}
\arguments{
\item{n}{Number of subjects.}
\item{nis}{Visits per subject (scalar or length-n vector).}
\item{seed}{RNG seed.}
\item{random_nis}{If TRUE, sample visit counts from \code{nis_range}.}
\item{nis_range}{Integer range for random visit counts.}
\item{scenario}{Simulation scenario (1 or 2).}
\item{C}{Number of ordinal categories.}
\item{re_dist}{Random-effects truth: \"normal\", \"mixture\", or \"mvt\".}
\item{mvt_df}{Degrees of freedom when \code{re_dist = \"mvt\"}.}
\item{rho}{Common RE correlation (default 0.5).}
\item{zi_intercept}{Optional override of zero-inflation intercept.}
\item{...}{Passed to \code{generate_data}.}
}
\value{
A list with class \code{zinbcdpmm_data}.
}

