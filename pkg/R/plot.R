plot.nls <- function(x, ...) {
    ff <- fitted(x)
    nn <- length(ff)
    rr <- residuals(x)
    df <- data.frame(x=c(fitted=ff,qnorm(ppoints(nn))), y=c(rr, sort(rr)),
                     f=factor(rep.int(c("Fitted values", "Normal quantiles"),
                     c(nn,nn))))
    if ("ylab" %in% names(list(...)))
        xyplot(y ~ x|f, df, type=c("g","p"), scales=list(x=list(relation="free")), xlab=NULL, ...)
    else
        xyplot(y ~ x|f, df, type=c("g","p"), scales=list(x=list(relation="free")),
               xlab=NULL, ylab="Residuals", ...)
}

## check if the results of confint.nls have a class.  If so check the class or use a null S3 default and a method.
cigrid <- function(ci, fun, npts=200) {
  vals <- as.data.frame(
    apply(ci, 1, function(v) seq(v[1],v[2],length.out=npts)))
  ans <- array(apply(do.call(expand.grid,vals),1,fun),
               c(npts,npts))
  attr(ans, "vals") <- vals
  class(ans) <- "cigrid"
  ans
}

contourplot.cigrid <- function(x, data, ...) {
  dots <- list(...)
  dnms <- names(dots)
  vals <- attr(x, "vals")
  vnms <- names(vals)
  if (!("xlab" %in% dnms)) dots$xlab <- vnms[1]
  if (!("ylab" %in% dnms)) dots$ylab <- vnms[2]
  dots$x <- unclass(x)
  dots$row.values <- vals[[1]]
  dots$column.values <- vals[[2]]
  dots$xlim <- range(vals[[1]])
  dots$ylim <- range(vals[[2]])
  do.call(contourplot, dots)
}
