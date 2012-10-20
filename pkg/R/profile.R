## functions and methods for profile.nls objects

as.data.frame.profile.nls <-
    function(x, row.names = NULL, optional = FALSE, ..., tauCoords=FALSE)
{
    fr   <- do.call("rbind", lapply(x, "[[", "par.vals"))
    tau  <- lapply(x, "[[", "tau")
    xnms <- names(x)
    pnms <- factor(rep.int(xnms, sapply(tau, length)))
    fr   <- as.data.frame(fr)
    tnm  <- ifelse("tau" %in% names(fr), ".tau", "tau")
    fr[[tnm]] <- unlist(tau)
    fr$.pnm <- pnms
    ss   <- split(fr, fr$.pnm)
    fwd  <- list()
    for (nm in xnms) {
        form <- eval(substitute(a ~ b,
                                list(a=as.name(tnm), b=as.name(nm))))
        fwd[[nm]] <- interpSpline(form, ss[[nm]])
    }
    attr(fr,"forward") <- fwd
    attr(fr,"backward") <- lapply(fwd, backSpline)
    class(fr) <- c("profile.frame", "data.frame")
    if (!tauCoords) return(fr)
    fr[[tnm]] <- NULL
    for (nm in xnms)
        fr[[nm]] <- zapsmall(predy(fwd[[nm]], fr[[nm]]),digits=12)
    class(fr) <- c("tau.frame", class(fr))
    fr
}

taugrid <- function(pr, fun, level=0.99, npts=200) {
    stopifnot(inherits(pr, "profile.nls"),
              0 < (level <- as.numeric(level[1])),
              level < 1,
              5 < (npts <- as.integer(npts[1])))
    orig <- attr(pr, "original.fit")
    bd   <- sqrt(2*qf(level, length(coef(orig)), df.residual(orig)))
    pnms <- names(coef(orig))
    names(pnms) <- pnms
    args <- c(lapply(pnms,function(nm) seq(-bd, bd, length.out=npts)),
              KEEP.OUT.ATTRS=FALSE)
    gr   <- do.call(expand.grid, args)
    pmat <- matrix(0., nrow=nrow(gr), ncol=length(pnms))
    colnames(pmat) <- unname(pnms)
    fr   <- as.data.frame(pr, tauCoords=TRUE)
    bac  <- attr(fr, "backward")
    for (nm in pnms) pmat[,nm] <- predy(bac[[nm]], gr[[nm]])
    gr$pmat <- pmat
    gr$RSS  <- apply(pmat, 1, fun)
    attr(gr, "tau.frame") <- fr
    gr
}

## extract only the y component from a prediction
predy <- function(sp, vv) predict(sp, vv)$y
## A lattice-based plot method for profile.nls objects
plot.profile.nls <-
    function (x, levels = sqrt(qf(pmax.int(0, pmin.int(1, conf)), 1, df[2])),
              conf = c(50, 80, 90, 95, 99)/100, 
              absVal = TRUE, ...) 
{
    df <- attr(x, "summary")$df
    levels <- sort(levels[is.finite(levels) && levels > 0])
    spl <- lapply(x, function(x)
                  splines::interpSpline(x$par.vals[, attr(x, "parameters")$par],
                                        x$tau))
    bspl <- lapply(spl, splines::backSpline)
    tau <- c(-rev(levels), 0, levels)
    fr <- data.frame(tau = rep.int(tau, length(x)),
                     pval = unlist(lapply(bspl, predy, tau)),
                     pnm = gl(length(x), length(tau), labels = names(x)))
    ylab <- expression(tau)
    if (absVal) {
        fr$tau <- abs(fr$tau)
        ylab <- expression("|" * tau * "|")
    }
    xyplot(tau ~ pval | pnm, fr,
           scales = list(x = list(relation = 'free')),
           ylab = ylab, xlab = "", panel = function(x, y, ...)
       {
           pfun <- function(x) predy(spl[[panel.number()]], x)
           panel.grid(h = -1, v = -1)
           lsegments(x, y, x, 0, ...)
           if (absVal) {
               lsegments(x, y, rev(x), y)
               pfun <- function(x) abs(predy(spl[[panel.number()]], x))
           } else {
               panel.abline(h = 0, ...)
           }
           panel.curve(pfun, ...)
       }, ...)
}

## convert the x-cosine and y-cosine to an average and difference,
## ensuring that the difference is positive by flipping signs if
## necessary
ad <- function(xc, yc)
{
    a <- (xc + yc)/2
    d <- (xc - yc)
    cbind(ifelse(d > 0, a, -a), abs(d))
}

## convert d versus a (as an xyVector) and level to a matrix of taui and tauj
tauij <- function(xy, lev) lev * cos(xy$x + outer(xy$y/2, c(-1, 1)))

## safe arc-cosine
sacos <- function(x) acos(pmax.int(-0.999, pmin.int(0.999, x)))

cont <- function(sij, sji, levels, nseg = 101)
{
    ada <- array(0, c(length(levels), 2, 4))
    ada[, , 1] <- ad(0, sacos(predy(sij,  levels)/levels))
    ada[, , 2] <- ad(sacos(predy(sji, levels)/levels), 0)
    ada[, , 3] <- ad(pi, sacos(predy(sij, -levels)/levels))
    ada[, , 4] <- ad(sacos(predy(sji, -levels)/levels), pi)
    pts <- array(0, c(length(levels), nseg + 1, 2))
    for (i in seq_along(levels))
        pts[i, ,] <- tauij(predict(periodicSpline(ada[i, 1, ], ada[i, 2, ]),
                                   nseg = nseg), levels[i])
    levs <- c(-rev(levels), 0, levels)
    list(tki = predict(sij, levs), tkj = predict(sji, levs), pts = pts)
}

splom.profile.nls <-
    function (x, data, ## unused - only for compatibility with generic
              levels = sqrt(2 *
              qf(pmax.int(0, pmin.int(1, conf)), 2, df[2])),
              conf = c(50, 80, 90, 95, 99)/100, ...)
{
    df <- attr(x, "summary")$df
    levels <- sort(levels[is.finite(levels) && levels > 0])
    mlev <- max(levels)
    spl <- lapply(x, function(x)
                  interpSpline(x$par.vals[, attr(x, "parameters")$par], x$tau))
    frange <- sapply(spl, function(x) range(x$knots))
    bsp <- lapply(spl, backSpline)
    brange <- sapply(bsp, function(x) range(x$knots))
    pfr <- do.call(cbind, lapply(bsp, predy, c(-mlev, mlev)))
    pfr[1, ] <- pmax.int(pfr[1, ], frange[1, ], na.rm = TRUE)
    pfr[2, ] <- pmin.int(pfr[2, ], frange[2, ], na.rm = TRUE)
    nms <- names(spl)
    np <- length(nms)

    ## Create data frame fr of par. vals in tau coordinates
    fr <- as.data.frame(x)
    for (nm in nms) fr[[nm]] <- predy(spl[[nm]], fr[[nm]])
    ## create a list of lists with the names of the parameters
    traces <- lapply(x, function(el) lapply(x, function(el1) list()))
    for (j in seq_along(nms)[-1]) {
        for (i in seq_len(j - 1)) {
            fri <- subset(fr, .pnm == nms[i])
            sij <- interpSpline(fri[ , i], fri[ , j])
            frj <- subset(fr, .pnm == nms[j])
            sji <- interpSpline(frj[ , j], frj[ , i])
            ll <- cont(sij, sji, levels)
            traces[[j]][[i]] <- list(sij = sij, sji = sji, ll = ll)
        }
    }
    ## panel function for lower triangle
    lp <- function(x, y, groups, subscripts, i, j, ...) {
        tr <- traces[[j]][[i]]
        grid::pushViewport(viewport(xscale = c(-1.07, 1.07) * mlev,
                              yscale = c(-1.07, 1.07) * mlev))
        dd <- sapply(current.panel.limits(), diff)/50
        psij <- predict(tr$sij)
        ll <- tr$ll
        ## now do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(psij$y, psij$x, ...)
        llines(predict(tr$sji), ...)
        with(ll$tki, lsegments(y - dd[1], x, y + dd[1], x, ...))
        with(ll$tkj, lsegments(x, y - dd[2], x, y + dd[2], ...))
        for (k in seq_along(levels)) llines(ll$pts[k, , ], ...)
        grid::popViewport(1)
    }
    ## panel function for upper triangle
    up <- function(x, y, groups, subscripts, i, j, ...) {
        ## panels are transposed so reverse i and j
        jj <- i
        ii <- j
        tr <- traces[[jj]][[ii]]
        ll <- tr$ll
        pts <- ll$pts
        limits <- current.panel.limits()
        psij <- predict(tr$sij)
        psji <- predict(tr$sji)
        ## do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(predy(bsp[[ii]], psij$y), predy(bsp[[jj]], psij$x), ...)
        llines(predy(bsp[[ii]], psji$x), predy(bsp[[jj]], psji$y), ...)
        for (k in seq_along(levels))
            llines(predy(bsp[[ii]], pts[k, , 1]),
                   predy(bsp[[jj]], pts[k, , 2]), ...)
    }
    dp <- function(x = NULL,            # diagonal panel
                   varname = NULL, limits, at = NULL, lab = NULL,
                   draw = TRUE,

                   varname.col = add.text$col,
                   varname.cex = add.text$cex,
                   varname.lineheight = add.text$lineheight,
                   varname.font = add.text$font,
                   varname.fontfamily = add.text$fontfamily,
                   varname.fontface = add.text$fontface,

                   axis.text.col = axis.text$col,
                   axis.text.alpha = axis.text$alpha,
                   axis.text.cex = axis.text$cex,
                   axis.text.font = axis.text$font,
                   axis.text.fontfamily = axis.text$fontfamily,
                   axis.text.fontface = axis.text$fontface,

                   axis.line.col = axis.line$col,
                   axis.line.alpha = axis.line$alpha,
                   axis.line.lty = axis.line$lty,
                   axis.line.lwd = axis.line$lwd,
                   i, j, 
                   ...)
    {
        n.var <- eval.parent(expression(n.var))
        add.text <- trellis.par.get("add.text")
        axis.line <- trellis.par.get("axis.line")
        axis.text <- trellis.par.get("axis.text")
        if (!is.null(varname))
            grid::grid.text(varname,
                            gp =
                            gpar(col = varname.col,
                                 cex = varname.cex,
                                 lineheight = varname.lineheight,
                                 fontface = lattice:::chooseFace(varname.fontface,
                                 varname.font),
                                 fontfamily = varname.fontfamily))
        if (draw) {
            at <- pretty(limits)
            sides <- c("left", "top")
            if (j == 1) sides <- "top"
            if (j == n.var) sides <- "left"
            for (side in sides)
                panel.axis(side = side,
                           at = at,
                           labels = format(at, trim = TRUE),
                           ticks = TRUE,
                           check.overlap = TRUE,
                           half = side == "top" && j > 1,

                           tck = 1, rot = 0,

                           text.col = axis.text.col,
                           text.alpha = axis.text.alpha,
                           text.cex = axis.text.cex,
                           text.font = axis.text.font,
                           text.fontfamily = axis.text.fontfamily,
                           text.fontface = axis.text.fontface,

                           line.col = axis.line.col,
                           line.alpha = axis.line.alpha,
                           line.lty = axis.line.lty,
                           line.lwd = axis.line.lwd)
            lims <- c(-1.07, 1.07) * mlev
            grid::pushViewport(viewport(xscale = lims, yscale = lims))
            side <- ifelse(j == 1, "right", "bottom")
            which.half <- ifelse(j == 1, "lower", "upper")
            at <- pretty(lims)
            panel.axis(side = side, at = at, labels = format(at, trim = TRUE),
                       ticks = TRUE, half = TRUE, which.half = which.half,
                       tck = 1, rot = 0,

                       text.col = axis.text.col,
                       text.alpha = axis.text.alpha,
                       text.cex = axis.text.cex,
                       text.font = axis.text.font,
                       text.fontfamily = axis.text.fontfamily,
                       text.fontface = axis.text.fontface,

                       line.col = axis.line.col,
                       line.alpha = axis.line.alpha,
                       line.lty = axis.line.lty,
                       line.lwd = axis.line.lwd)
            grid::popViewport(1)
        }
    }

    splom(~ pfr, lower.panel = lp, upper.panel = up, diag.panel = dp, ...)
}
