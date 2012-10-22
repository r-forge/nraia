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
cigrid <- function(x, fun, npts=200, ...) UseMethod("cigrid")

cigrid.matrix <- function(x, fun, npts=200, ...) {
    ## Create a data frame of the parameter values on the grid
    vals <- as.data.frame(
        apply(x, 1, function(v) seq(v[1], v[2], length.out=npts)))

    ## Evaluate fun on the grid points and convert to an array
    ans <- array(apply(do.call(expand.grid,vals),1,fun),
                 rep.int(npts,length(vals)))
    attr(ans, "vals") <- vals
    class(ans) <- "cigrid"
    ans
}

cigrid.profile.nls <- function(x, fun, npts=200, maxconf=0.99, ...) {
    nu   <- df.residual(attr(x, "original.fit"))
    P    <- length(x)
    cis  <- confint(x, level=pf(P*qf(maxconf, P, nu), 1, nu))
    if (any(is.na(cis))) stop("write code to fix NA's in confints")
    ans  <- cigrid(cis, fun, npts, ...)
    class(ans) <- c("cigrid.profile", class(ans))
    attr(ans, "profile") <- as.data.frame(x)
    attr(ans, "fdf")     <- c(P, nu)
    attr(ans, "deviance")<- deviance(attr(x, "original.fit"))
    attr(ans, "coef")    <- coef(attr(x, "original.fit"))
    ans
}

contourplot.cigrid.profile <- function(x, data, levs=c(0.5,0.8,0.9,0.95,0.99), ...) {
    cc   <- attr(x, "coef")
    dev  <- attr(x, "deviance")
    fdf  <- attr(x, "fdf")
    frm  <- attr(x, "profile")
    frms <- split(frm, frm$.pnm)
    vals <- attr(x, "vals")
    forw <- attr(frm, "forward")
    bakw <- attr(frm, "backward")
    pnms <- names(vals)
    P    <- length(vals)
    Pm1  <- P - 1L
    ww   <- 1./Pm1                    # widths of viewports in plotViewport
    cens <- (1:Pm1)*ww - 1./(2*Pm1)   # centers of subviewports
    grid.newpage()
    pushViewport(plotViewport(c(5,4,2,2), name="canvas"))
    for (j in 1:Pm1) {
        for (i in (j+1):P) {
            nms <- pnms[c(j,i)]
            pushViewport(viewport(cens[j], cens[P + 1 - i], ww, ww))
            pushViewport(dataViewport(vals[[j]], vals[[i]], extension=0, name="plotRegion"))
            if (j==1L) {
                grid.yaxis()
                grid.text(pnms[i], x=unit(-3,"lines"), rot=90)
            }
            if (i==P) {
                grid.xaxis()
                grid.text(pnms[j], y=unit(-3,"lines"))
            }
            if (i==(j+1)) {
                grid.xaxis(label=FALSE,main=FALSE)
                grid.yaxis(label=FALSE,main=FALSE)
            }
            grid.rect()
            lapply(contourLines(vals[[j]], vals[[i]], apply(x, c(j,i), min),
                                levels = dev*(1+fdf[1]*qf(levs,fdf[1],fdf[2])/fdf[2])),
                   function(cl) grid.lines(cl$x, cl$y, default.units="native"))
            grid.points(cc[j], cc[i], pch="+", gp=gpar(cex=1.2, col='black'))
            pushViewport(dataViewport(vals[[j]], vals[[i]], extension=0, name="plotRegion", clip="on"))
            frmi <- frms[[pnms[i]]]
            frmj <- frms[[pnms[j]]]
            grid.lines(frmi[[j]], frmi[[i]], gp=gpar(lty=5), default.units="native")
            grid.lines(frmj[[j]], frmj[[i]], gp=gpar(lty=6), default.units="native")
            popViewport(3)
        }
    }
}
