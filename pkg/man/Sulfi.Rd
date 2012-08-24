\name{Sulfi}
\alias{Sulfi}
\non_function{}
\title{Pharmacokinetics of sulfisoxazole}
\description{
  The \code{Sulfi} data frame has 12 rows and 2 columns from an experiment
  on the pharmacokinetics of sulfisoxazole.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{
      a numeric vector giving the time since drug administration (min).
    }
    \item{conc}{
      a numeric vector giving the observed concentration of
      sulfisoxazole (\eqn{\mu}{u}g/ml).
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
fm1 <- nls(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2), Sulfi)
coef(summary(fm1))
plotfit(fm1,
     xlab = "Time since drug administration (min)",
     ylab = expression("Sulfisoxazole concentration ("*mu*"g/ml)"))
}
\keyword{datasets}
