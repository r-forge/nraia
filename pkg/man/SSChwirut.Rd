\name{SSChwirut}
\title{Chiwrut model for ultrasonic response}
\usage{
SSChwirut(input, lrc, b0, b1)
}
\alias{SSChwirut}
\arguments{
 \item{input}{a numeric vector of values at which to evaluate the model.}
 \item{lrc}{a numeric parameter representing the logarithm of the rate constant.}
 \item{b0}{a numeric parameter representing the intercept in linear
   denominator expression.}
 \item{b1}{a numeric parameter representing the slope in linear
   denominator expression.}
}
\description{
  This \code{\link{selfStart}} model evaluates the Chwirut model,
  for an ultrasonic response.
  It has an \code{initial} attribute that creates initial estimates of
  the parameters \code{lrc}, \code{b0} and \code{b1}.  The model
  function also returns a \code{"gradient"} attribute.
}
\value{
  a numeric vector of the same length as \code{input}.  It is the value of
  the expression \code{exp(-exp(lrc)*input)/(b0+b1*input)}.
  If all of the arguments \code{lrc}, \code{b0} and \code{b1}
  are names of objects, the gradient matrix with respect to
  these names is attached as an attribute named \code{gradient}.
}
\seealso{\code{\link{nls}}, \code{\link{selfStart}}
}
\examples{
if (require(NISTnls)) {
   coef(summary(nls(y ~ SSChwirut(x, lrc, b0, b1), Chwirut1)))
   coef(summary(nls(y ~ SSChwirut(x, lrc, b0, b1), Chwirut2)))
}
}
\keyword{models}
