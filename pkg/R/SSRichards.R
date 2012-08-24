## self-starting model function for the Richards growth model
SSRichards <-
    selfStart(~ Asym * (1+exp((xmid - input)/scal))^(-exp(-lpow)),
              function(mCall, data, LHS)
          {
              linit <- unname(attr(SSlogis, "initial")(mCall, data, LHS))
              xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
              if (nrow(xy) < 5) {
                  stop("too few distinct input values to fit a logistic model")
              }
              pars <- unname(coef(nls(y ~ (1+exp((xmid - x)/scal))^(-exp(-lpow)),
                                      xy, c(xmid = linit[2], scal = linit[3],
                                            lpow = 0.001), alg = "plinear")))
              value <- pars[c(4,1:3)]
              names(value) <- mCall[c("Asym", "xmid", "scal", "lpow")]
              value
          },
              c("Asym", "xmid", "scal", "lpow"),
              function(input, Asym, xmid, scal, lpow) {})

SSChwirut <-
    selfStart(~ exp(-exp(lrc)*input)/(b0 + b1*input),
              function(mCall, data, LHS)
          {
              xy   <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
              if (nrow(xy) < 4) {
                  stop("too few distinct input values to fit the Chwirut model")
              }
              rc1  <- -coef(lm(log(y) ~ x, xy))[2]
              pars <- coef(nls(y ~ exp(-exp(lrc)*x)/(1+p3*x), xy,
                               c(lrc=log(unname(rc1)), p3=1/mean(xy$x)), alg="plinear"))
              value <- c(pars[1], c(1, pars[2])/pars[3])
              names(value) <- mCall[c("lrc", "b0", "b1")]
              value
          },
              c("lrc", "b0", "b1"),
              c("input", "lrc", "b0", "b1"))
