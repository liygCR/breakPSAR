# hypothesis test of a singe structral break by the likelihood ratio method.
#
# This function provides An quick hypothesis test of a singe structral break
# by the likelihood ratio method proposed by Davis et al. (1995). Typcally for
# linear model lm(formula=y1~x1-1). Confidence level is 0.05.
#
#
#
# @param y1 an autoregressive time series.
# @param x1 design matrix without the intercept.
# @return Returns an object with
# @return \item{cp}{the estimate of chang-point if there is one and 0 otherwise.}
# @export css
# @references {Davis et al. (1995), Testing for a change in the parameter values and
#  order of an autoregressive model, \emph{The Annals of Statistics} 23(1), 282-304}
#

#rm(list=ls())

css <- function(y1, x1) {
  #
  n <- dim(x1)[1]
  p <- dim(x1)[2]
  bn <- (2 * log(log(n)) + (p + 1)/2 * log(log(log(n))) - log(gamma((p +
                                                                       1)/2)))^2/(2 * log(log(n)))
  an <- sqrt(bn/(2 * log(log(n))))

  # Recursive Calculation of the likelihood ratio statistic

  S <- t(y1[1:(n)]) %*% x1[1:(n), ]
  C1 <- solve(t(x1[1:(n), ]) %*% x1[1:(n), ], t(S))
  Q <- S %*% C1
  sigma <- (sum(y1^2) - Q)/n
  Dic <- c(rep(0, n))
  C0 <- solve(t(x1[1:(2 * p), ]) %*% x1[1:(2 * p), ])
  C1 <- solve(t(x1[(2 * p + 1):n, ]) %*% x1[(2 * p + 1):n, ])
  S0 <- t(y1[1:(2 * p)]) %*% x1[1:(2 * p), ]
  S1 <- t(y1[(2 * p + 1):(n)]) %*% x1[(2 * p + 1):(n), ]
  for (i in (2 * p + 1):(n - (2 * p + 1))) {
    xi <- matrix(x1[i, ], 1, p)
    C0 <- C0 - C0 %*% t(xi) %*% xi %*% C0/c(1 + xi %*% C0 %*% t(xi))
    C1 <- C1 + C1 %*% t(xi) %*% xi %*% C1/c(1 - xi %*% C1 %*% t(xi))
    S0 <- S0 + y1[i] * xi
    S1 <- S1 - y1[i] * xi
    Q1 <- S0 %*% C0 %*% t(S0)
    Q2 <- S1 %*% C1 %*% t(S1)
    Dic[i] <- Q1 + Q2
  }

  cp <- which(Dic == max(Dic))
  d <- (Dic[cp] - Q)/sigma
  if (((d - bn)/an) >= 2 * log(-2/(log(1 - 0.05))))
    return(cp) else return(0)

}

#' Use ASCAD or AMCP to estimate the chang-points.
#'
#' This function provides a novel and fast methodology for simultaneous multiple
#' structural break estimation and variable selection for nonstationary time series
#' models by using ASCAD or AMCP to estimate the chang-points by Jin, Shi and Wu (2011).
#'
#'
#' @param Y an autoregressive time series.
#' @param method the method to be used by mcp or scad.
#' See \code{\link[plus]{plus}} in R packages \pkg{plus} for details.
#' @param p an upper bound of autoregressive orders.
#' @param lam.d an parameter of the refining procedure. Suggest lam.d=0.05 if
#' the length of Y is smaller than 5000 and lam.d=0.01 otherwise.
#' @return pluscpv returns an object of class "pluscpv".
#' An object of class "pluscpv" is a list containing the following components:
#' @return \item{variable.selection}{a matrix of variable selection.
#' Element of the	matrix is 1 if the variable is selected and 0	otherwise.}
#' @return \item{change.points}{estimators of change points.}
#' @export pluscpv
#' @import plus
#' @seealso plus
#' @references {Jin, B., Shi, X., and Wu, Y. (2013). A novel and fast methodology for
#'  simultaneous multiple structural break estimation and variable selection for
#'  nonstationary time series models. \emph{Statistics and Computing}, 23(2), 221-231.}
#'

# Use ASCAD or AMCP to estimate the chang-points by Jin,Shi and Wu (2011).
# "A novel and fast methodology for simultaneous multiple structural break
# estimation and variable selection for nonstationary time series models"


pluscpv <- function(Y, method = c("mcp", "scad"), p, lam.d) {



  library("plus")  #call the plus package build in R by
  # ZHANG(2010, Nearly unbiased variable selection under minimax concave
  # penalty. The Annals of Statistics 38(2), 894-942)
  n = length(Y)
  m = 4 * ceiling(sqrt(n - p))
  q = floor((n - p)/m)
  K_temp = matrix(0, nrow = q, ncol = q, byrow = TRUE)
  ########### transform###################

  X_temp = matrix(1, n - p, (p + 1))
  for (i in 1:(n - p)) {
    for (j in 2:(p + 1)) {
      X_temp[i, j] = Y[p + i - (j - 1)]
    }
  }
  Y_temp = c(Y[(p + 1):n])
  for (i in 1:q) K_temp[i, 1:i] = rep(1, i)




  x = NULL
  y = NULL
  x[[1]] = X_temp[1:((n - p - (q - 1) * m)), ]
  y[[1]] = Y_temp[1:((n - p - (q - 1) * m))]

  for (i in 2:q) {
    x[[i]] = X_temp[(n - p - (q - i + 1) * m + 1):((n - p - (q - i) *
                                                      m)), ]
    y[[i]] = Y_temp[(n - p - (q - i + 1) * m + 1):((n - p - (q - i) *
                                                      m))]
  }
  X_temp1 <- sapply(1:length(x), function(j, mat, list) kronecker(K_temp[j,
                                                                         , drop = FALSE], x[[j]]), mat = K_temp, list = x, simplify = F)
  X = do.call("rbind", X_temp1)  #X is a design matrix after transformation




  #################### mcp#########################
  if (method == "mcp") {
    object <- plus(X, Y_temp, method = "mc+", gamma = 2.4, intercept = F,
                   normalize = F, eps = 1e-50)

    bic = log(dim(X)[1]) * object$dim + dim(X)[1] * log(as.vector((1 -
                                                                     object$r.square) * sum(Y_temp^2))/length(Y_temp))
    step.bic = which.min(bic)
    mcp.coef <- coef(object, lam = object$lam[step.bic])

    mcp.coef.s = sum(abs(mcp.coef))
    if (mcp.coef.s == 0) {
      mcp.coef.v = mcp.coef
    }
    if (mcp.coef.s > 0) {
      object <- plus(diag(length(mcp.coef)), mcp.coef, method = "mc+",
                     gamma = 2.4, intercept = F)
      mcp.coef.v = coef(object, lam = lam.d)
    }

    mcp.coef.ar = c(rep(0, q * (p + 1)))
    mcp.coef.ar[which(mcp.coef.v[1:(p + 1)] != 0)] = mcp.coef[which(mcp.coef.v[1:(p +
                                                                                    1)] != 0)]


    mcp.coef.v.m = abs(matrix(c(mcp.coef.v), q, (p + 1), byrow = T))
    mcp.coef.m = c(apply(mcp.coef.v.m, 1, max))




    mcp.cp1 = which(mcp.coef.m != 0)

    mcp.cp1 = mcp.cp1[mcp.cp1 > 1]

    d1 = length(mcp.cp1)

    if (d1 == 0) {
      mcpcss.cp = 0
      mcp.v.c = mcp.coef.v.m[1, ]
      mcp.v.c[which(mcp.v.c != 0)] = 1
    }
    if (d1 >= 1) {
      mcpcss.cp = NULL
      mcp.v.i = 1
      j = 0
      for (i in 1:d1) {
        if (j == 1) {
          if (mcp.cp1[i] - mcp.cp1[i - 1] == 1) {
            j = 0
            next
          }
        }
        y1 = c(y[[mcp.cp1[i] - 1]], y[[mcp.cp1[i]]])
        x1 = rbind(x[[mcp.cp1[i] - 1]], x[[mcp.cp1[i]]])
        cp = css(y1, x1)  #call css function for single change point test
        if (cp == 0)
          next
        if (cp != 0) {
          if (mcp.cp1[i] == 2)
            mcpcss.cp = c(mcpcss.cp, cp)
          j = 1
          if (mcp.cp1[i] > 2)
            mcpcss.cp = c(mcpcss.cp, n - (q - mcp.cp1[i] + 2) *
                            m + cp)
          j = 1
          mcp.v.i = c(mcp.v.i, mcp.cp1[i])
          s.i1 = c((mcp.v.i[length(mcp.v.i)] * (p + 1) - p):(mcp.v.i[length(mcp.v.i)] *
                                                               (p + 1)))
          s.i0 = c((mcp.v.i[length(mcp.v.i) - 1] * (p + 1) - p):(mcp.v.i[length(mcp.v.i) -
                                                                           1] * (p + 1)))
          mcp.coef.si1 = mcp.coef.v[s.i1]
          mcp.coef.si0 = mcp.coef[s.i1]
          mcp.coef.si1[which(mcp.coef.si1 != 0)] = mcp.coef.si0[which(mcp.coef.si1 !=
                                                                        0)]
          mcp.coef.ar[s.i1] = mcp.coef.ar[s.i0] + mcp.coef.si1
        }
      }
      # do variable selection
      mcp.v = mcp.coef.v.m[mcp.v.i, ]
      mcp.v[which(mcp.v != 0)] = 1
      object <- plus(diag(length(mcp.coef.ar)), mcp.coef.ar, method = "mc+",
                     gamma = 2.4, intercept = F)
      mcp.coef.v.c = coef(object, lam = lam.d)
      mcp.coef.v.m = abs(matrix(c(mcp.coef.v.c), q, (p + 1), byrow = T))
      mcp.coef.m = c(apply(mcp.coef.v.m, 1, max))
      mcp.cp1 = which(mcp.coef.m != 0)
      mcp.v.c = mcp.coef.v.m[mcp.cp1, ]
      mcp.v.c[which(mcp.v.c != 0)] = 1
    }
    mcp.v.c = as.matrix(mcp.v.c)
    if (length(mcpcss.cp) == 1) {
      if (mcpcss.cp == 0)
        mcp.v.c = t(mcp.v.c)
    }
    dimc = dim(mcp.v.c)
    if (dimc[1] == 1)
      row_names = "model_1" else {
        row_names = c(rep(0, dimc[1]))
        for (i in 1:dimc[1]) row_names[i] = paste("model", i, sep = "_")
      }
    if (dimc[2] == 1)
      col_names = "intercept" else {
        col_names = c(rep(0, dimc[2]))
        col_names[1] = "intercept"
        for (i in 2:dimc[2]) col_names[i] = paste("lag", i - 1, sep = "_")
      }
    dimnames(mcp.v.c) = list(row_names, col_names)

    return(list(variable.selection = mcp.v.c, change.points = mcpcss.cp))
  }
  #################### scad#########################

  if (method == "scad") {
    object <- plus(X, Y_temp, method = "scad", gamma = 2.4, intercept = F,
                   normalize = F, eps = 1e-50)

    bic = log(dim(X)[1]) * object$dim + dim(X)[1] * log(as.vector((1 -
                                                                     object$r.square) * sum(Y_temp^2))/length(Y_temp))
    step.bic = which.min(bic)
    scad.coef <- coef(object, lam = object$lam[step.bic])

    scad.coef.s = sum(abs(scad.coef))
    if (scad.coef.s == 0) {
      scad.coef.v = scad.coef
    }
    if (scad.coef.s > 0) {
      object <- plus(diag(length(scad.coef)), scad.coef, method = "scad",
                     gamma = 2.4, intercept = F)
      scad.coef.v = coef(object, lam = lam.d)
    }

    scad.coef.ar = c(rep(0, q * (p + 1)))
    scad.coef.ar[which(scad.coef.v[1:(p + 1)] != 0)] = scad.coef[which(scad.coef.v[1:(p +
                                                                                        1)] != 0)]

    scad.coef.v.m = abs(matrix(c(scad.coef.v), q, (p + 1), byrow = T))
    scad.coef.m = c(apply(scad.coef.v.m, 1, max))




    scad.cp1 = which(scad.coef.m != 0)

    scad.cp1 = scad.cp1[scad.cp1 > 1]

    d1 = length(scad.cp1)

    if (d1 == 0) {
      scadcss.cp = 0
      scad.v.c = scad.coef.v.m[1, ]
      scad.v.c[which(scad.v.c != 0)] = 1
    }
    if (d1 >= 1) {
      scadcss.cp = NULL
      scad.v.i = 1
      j = 0
      for (i in 1:d1) {
        if (j == 1) {
          if (scad.cp1[i] - scad.cp1[i - 1] == 1) {
            j = 0
            next
          }
        }
        y1 = c(y[[scad.cp1[i] - 1]], y[[scad.cp1[i]]])
        x1 = rbind(x[[scad.cp1[i] - 1]], x[[scad.cp1[i]]])
        cp = css(y1, x1)  #call css function for single change point test
        if (cp == 0)
          next
        if (cp != 0) {
          if (scad.cp1[i] == 2)
            scadcss.cp = c(scadcss.cp, cp)
          j = 1
          if (scad.cp1[i] > 2)
            scadcss.cp = c(scadcss.cp, n - (q - scad.cp1[i] + 2) *
                             m + cp)
          j = 1
          scad.v.i = c(scad.v.i, scad.cp1[i])
          s.i1 = c((scad.v.i[length(scad.v.i)] * (p + 1) - p):(scad.v.i[length(scad.v.i)] *
                                                                 (p + 1)))
          s.i0 = c((scad.v.i[length(scad.v.i) - 1] * (p + 1) -
                      p):(scad.v.i[length(scad.v.i) - 1] * (p + 1)))
          scad.coef.si1 = scad.coef.v[s.i1]
          scad.coef.si0 = scad.coef[s.i1]
          scad.coef.si1[which(scad.coef.si1 != 0)] = scad.coef.si0[which(scad.coef.si1 !=
                                                                           0)]
          scad.coef.ar[s.i1] = scad.coef.ar[s.i0] + scad.coef.si1
        }
      }
      # do variable selection
      scad.v = scad.coef.v.m[scad.v.i, ]
      scad.v[which(scad.v != 0)] = 1
      object <- plus(diag(length(scad.coef.ar)), scad.coef.ar, method = "scad",
                     gamma = 2.4, intercept = F)
      scad.coef.v.c = coef(object, lam = lam.d)
      scad.coef.v.m = abs(matrix(c(scad.coef.v.c), q, (p + 1), byrow = T))
      scad.coef.m = c(apply(scad.coef.v.m, 1, max))
      scad.cp1 = which(scad.coef.m != 0)
      scad.v.c = scad.coef.v.m[scad.cp1, ]
      scad.v.c[which(scad.v.c != 0)] = 1
    }
    scad.v.c = as.matrix(scad.v.c)
    if (length(scadcss.cp) == 1) {
      if (scadcss.cp == 0)
        scad.v.c = t(scad.v.c)
    }
    dimc = dim(scad.v.c)
    if (dimc[1] == 1)
      row_names = "model_1" else {
        row_names = c(rep(0, dimc[1]))
        for (i in 1:dimc[1]) row_names[i] = paste("model", i, sep = "_")
      }
    if (dimc[2] == 1)
      col_names = "intercept" else {
        col_names = c(rep(0, dimc[2]))
        col_names[1] = "intercept"
        for (i in 2:dimc[2]) col_names[i] = paste("lag", i - 1, sep = "_")
      }
    dimnames(scad.v.c) = list(row_names, col_names)

    return(list(variable.selection = scad.v.c, change.points = scadcss.cp))
  }
}

