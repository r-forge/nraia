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
