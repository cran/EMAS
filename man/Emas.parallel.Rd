% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Emas.parallel.R
\name{Emas.parallel}
\alias{Emas.parallel}
\title{Epigenome-Wide Mediation Analysis Study: Parallel multiple mediation model}
\usage{
Emas.parallel(
  data,
  M.matrix,
  id = "",
  x = "",
  y = "",
  x.cov = c(),
  y.cov = c(),
  m.cov = c(),
  m.cor = TRUE,
  boot = FALSE,
  lavaan = FALSE,
  ...
)
}
\arguments{
\item{data}{A data.frame included id, x, y, x.cov, y.cov, m.cov.}

\item{M.matrix}{A matrix with the CpG information screened from EMAS., maybe a M-value matrix or a beta value matrix.}

\item{id}{Variable name of the id.}

\item{x, y}{Variable name of exposure(x) and outcome(y).}

\item{x.cov}{Variable names of covariates related to exposure(x).}

\item{y.cov}{Variable names of covariates related to outcome(y).}

\item{m.cov}{Variable names of covariates related to mediator(m).}

\item{m.cor}{A logical value. If 'TRUE', the mediators in the parallel multiple mediation model are set to correlate with each other.}

\item{boot}{A logical value or a numeric value. If a numeric value, the number for bootstrap.}

\item{lavaan}{A logical value. If 'TRUE', a lavaan object will be given.}

\item{...}{Other arguments passed to \code{\link{sem}} from \code{\link{lavaan}} package.}
}
\value{
\code{Emas.parallel} returns a data.frame with the average mediation effects(AME), 
   average direct effects(ADE), and total effects(Tot). If \code{lavaan} is 'TRUE', a lavaan object will be given.
}
\description{
This function can perform the parallel multiple mediation model after the Epigenome-Wide Mediation Analysis Study (EMAS).
}
\details{
This function can perform the parallel multiple mediation model after the Epigenome-Wide Mediation Analysis Study (EMAS) to further explore the potential parallel mediating CpG sites of exposure variables affecting outcome variables.
}
\examples{
\donttest{data(data.m)
data(Mvalue)
EP.result <- Emas.parallel(data.m, Mvalue, 
                           id = "ID", x = "x", y = "y", 
                           x.cov = c("age", "gender"), 
                           y.cov = c("age", "gender"), 
                           m.cov = c("age", "gender", "CD8T", "CD4T"), 
                           m.cor = TRUE, boot = FALSE, lavaan = FALSE)}
}
\references{
Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling.
\emph{Journal of Statistical Software}, 48(2), 1–36. \doi{doi:10.18637/jss.v048.i02}.
}
\author{
Xiuquan Nie, niexiuquan1995@foxmail.com
}
