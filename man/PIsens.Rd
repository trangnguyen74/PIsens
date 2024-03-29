% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/C-sens.R
\name{PIsens}
\alias{PIsens}
\title{Sensitivity analysis for violation of principal ignorability (PI)}
\usage{
PIsens(
  data,
  estimator,
  targeted,
  Z.form,
  C.form,
  Y.form,
  Y.model,
  Y.bounds,
  Y.FUN = NULL,
  sens.type,
  sens.range = NULL,
  sens.step = NULL,
  boot.num = 999,
  boot.seed = 12345,
  boot.cont.wt = TRUE,
  double.boot = FALSE,
  BCa.CI = TRUE
)
}
\arguments{
\item{data}{A data frame that has been processed by function \code{declare_data()}.}

\item{estimator}{One of the estimators listed in the paper: "epi", "pimu" or "emu" (for the simple estimator \eqn{e\pi}, \eqn{\pi\mu} or \eqn{e\mu}), or "MS", "IF" or "IFH" (for the multi-step, IF-as-estimating-function, or Hajek-ized estimator).}

\item{targeted}{If \code{TRUE}, the principal score model and the Y1 models are fit to inverse propensity score weighted treatment arm data; and the Y0 models are fit to control arm data that are inverse probability of control weighted and principal score weighted. If \code{FALSE}, such weighting is not used.}

\item{Z.form}{A character string specifying the formula for the propensity score model.}

\item{C.form}{A character string specifying the formula for the principal score model.}

\item{Y.form}{Either a character string specifying one formula to be used for all outcome model(s), or a character vector of length 3 specifying formulas for the outcome in compliers in the treatment arm (first element), in noncompliers in the treatment arm (second element), and in controls (third element.).}

\item{Y.model}{The type of model to be used for the outcome. Can be (1) a character string or (2) an object of class family. In case (1), if "linear", "gaussian" or "identity", will fit a linear model; if "binomial", "quasibinomial", "logit", "logistic" or "odds ratio", will fit a quasibinomial model with logit link; if "poisson",  "quasipoisson", "count" or "rate ratio", will fit a quasipoisson model with log link; if "Gamma" or "gamma", will fit a quasi model with log link and variance proportional to mean squared.}

\item{Y.bounds}{A two element vector for the lower and upper bounds of the outcome. Needed when modeling a bounded outcome by a quasibinomial model and/or when conducting a GOR-based sensitivity analysis.}

\item{Y.FUN}{Optional function for customized estimation of outcome mean functions. For an example, see function \code{jobsII_earn3_mus()} in this package, and see how it is used in the vignette. If use \code{Y.FUN}, arguments for outcome model forms and family do not need to be specified, and are ignored.}

\item{sens.type}{Name of the sensitivity parameter. Options are: "OR", "GOR", "MR", "MD", "beta.quant".}

\item{sens.range}{A numeric vector of length 2, holding the minimum and maximum values for the sensitivity parameter. If not specified, defaults to (1/3, 3) for \code{sens.type} being OR, GOR or MR; (-1, 1) for \code{sens.type} being MD; and (-.05, .05) for \code{sens.type} being "beta.quant".}

\item{sens.step}{The width of the step by which to increase the sensitivity parameter. If not specified, defaults to 0.05 for mean-based sensitivity analyses and 0.05 for distribution-based sensitivity analyses.}

\item{boot.num}{Number of bootstrap samples. Default is 999.}

\item{boot.seed}{Seed for reproducibility. Default is 12345.}

\item{boot.cont.wt}{If \code{TRUE} (default), implement continuous weights bootstrap (with weights based on the Dirichlet distribution). If \code{FALSE}, implement classic bootstrap (with weights based on the multinomial distribution). If}

\item{double.boot}{Whether to implement bias correction of point estimate by double bootstrap. If FALSE, implements bias correction of point estimate by single bootstrap only; this substantially reduces computational burden for the beta-quantile sensitivity analysis. Default is FALSE.}

\item{BCa.CI}{Whether to compute BCa confidence intervals. If FALSE, compute percentile intervals only, which is not recommended. However, BCa CIs are computationally burdensome for the beta-quantile sensitivity analysis. Default is TRUE.}
}
\description{
Sensitivity analysis for violation of principal ignorability (PI)
}
