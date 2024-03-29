plotfit <- function(fm, newpars=NULL, ...) UseMethod("plotfit")

## Determine the names of the x and y variables in a fitted model
vnm <- function(fm)
{
    mm <- fm$m
    cc <- fm$call
    pnms <- names(mm$getPars())
    form <- cc$formula
    rhsnms <- all.vars(form[[3]])
    vnms <- rhsnms[!(rhsnms %in% pnms)]
    if (length(vnms) > 1)
        stop("plotfit not yet implemented for >1 covariate")
    list(x = as.name(vnms), y = form[[2]])
}

## Create the predictor function from a fitted model
pfun <- function(fm, newpars=NULL)
{
    fm1 <- fm
    if (!is.null(newpars)) {
        stopifnot(is.numeric(newpars),
                  all(names(newpars) == names(pp <- fm$m$getPars())))
        if (fm1$m$setPars(newpars))
            warning("singular gradient matrix at newpars")
    }
    vnmx <- as.character(vnm(fm1)$x)
    function(x) {
        ll <- list(x)
        names(ll) <- vnmx
        predict(fm1, ll)
    }
}
    
plotfit.nls <- function(fm, newpars=NULL, ...)
{
    predfun <- pfun(fm, newpars)
    xyplot(eval(substitute(y ~ x, vnm(fm))), fm$m$getEnv(),
           panel = function(x, y, ...) {
               panel.grid(h = -1, v = -1)
               panel.points(x, y, ...)
               panel.curve(predfun, ...)
           }, ...)
}

plotfit.list <- function(fm, newpars=NULL, ...)
{
    if (!all(unlist(lapply(fm, inherits, "nls"))))
        stop("plotfit of a list must be a list of nls models")
    nms <- names(fm)
    pfuns <- lapply(fm, pfun)
    xyplot(eval(substitute(y ~ x, vnm(fm[[1]]))), fm[[1]]$m$getEnv(),
           panel = function(x, y, ...) {
               panel.grid(h = -1, v = -1)
               panel.points(x, y, ...)
               dots <- list(...)
               lims <- current.panel.limits()$x
               if (!is.null(dots$from)) lims[1] <- as.numeric(dots$from)[1]
               if (!is.null(dots$to)) lims[2] <- as.numeric(dots$to)[1]
               n <- 101
               if (!is.null(dots$n)) n <- as.integer(max(2, dots$n[1]))
               xv <- seq(lims[1], lims[2], len = n)
               ln <- trellis.par.get("superpose.line")
               for (i in seq_along(pfuns))
                   llines(xv, pfuns[[i]](xv), col = ln$col[i],
                          lty = ln$lty[i], lwd = ln$lwd[i])
           }, ...)
}
