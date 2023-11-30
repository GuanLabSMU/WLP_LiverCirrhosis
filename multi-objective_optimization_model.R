# Description
# The number of optimal components is obtained by using multi-objective optimization algorithm, 
# and then the components are selected from high to low according to the degree of C-T network.

mogptab <- function (factors, values, conditions) {
  levels <- max(conditions$priority)
  goals <- nrow(conditions)
  variables <- ncol(factors)
  objectives <- nrow(factors)
  nonbasics <- variables + objectives
  if (nrow(factors) != objectives) 
    stop("factors and conditions do not have the same number of objectives ")
  if (length(values) != objectives) 
    stop("conditions and values do not have the same number of objectives")
  te <- cbind(factors, diag(-1, objectives))
  tw <- matrix(0, nrow = levels, ncol = nonbasics)
  tu <- matrix(0, nrow = objectives, ncol = levels)
  for (goal in 1:goals) {
    i <- conditions[goal, "objective"]
    k <- conditions[goal, "priority"]
    w <- conditions[goal, "p"]
    u <- conditions[goal, "n"]
    s <- variables + i
    tw[k, s] <- w
    tu[i, k] <- u
  }
  col.headings <- c(paste("X", 1:variables, sep = ""), paste("P", 
                                                             1:objectives, sep = ""))
  row.headings <- paste("N", 1:objectives, sep = "")
  tb <- values
  ta <- rep(0, levels)
  ti <- matrix(0, nrow = levels, ncol = nonbasics)
  tab <- list(iter = 0, variables = variables, levels = levels, 
              objectives = objectives, nonbasics = nonbasics, level = 0, 
              te = te, tb = tb, tw = tw, tu = tu, ti = ti, ta = ta, 
              row.headings = row.headings, col.headings = col.headings)
  class(tab) <- "mogptab"
  return(tab)
}

check.tb <- function (tab) {
  for (i in 1:tab$objectives) {
    if (tab$tb[i] < 0) {
      for (j in 1 <- 1:tab$variables) {
        tab$te[i, j] <- -tab$te[i, j]
      }
      tab <- swp.headings(tab, i, i + tab$variables)
      tab <- swp.vec(tab, i, i + tab$variables)
    }
  }
  return(tab)
}

calc.ti.k <- function (tab, k) {
  tab$ti <- t(tab$tu) %*% tab$te - tab$tw
  return(tab)
}

calc.ta.k <- function (tab, k) {
  tab$ta <- t(tab$tu) %*% tab$tb
  return(tab)
}

swp.headings <- function (tab, nr, nc) {
  temp <- tab$col.headings[nc]
  tab$col.headings[nc] <- tab$row.headings[nr]
  tab$row.headings[nr] <- temp
  return(tab)
}

swp.vec <- function (tab, nr, nc) {
  for (k in 1:tab$levels) {
    temp <- tab$tu[nr, k]
    tab$tu[nr, k] <- tab$tw[k, nc]
    tab$tw[k, nc] <- temp
  }
  return(tab)
}

ev.mogp <- function (tab, k) {
  sp <- 0
  vsp <- 0
  for (s in 1:tab$nonbasics) {
    if (tab$ti[k, s] > 0) {
      if ((k == 1) || (neg.ind.rows(tab, k, s) == 0)) {
        if (tab$ti[k, s] > vsp) {
          sp <- s
          vsp <- tab$ti[k, s]
        }
      }
    }
  }
  return(sp)
}

mogpout <- function (tab, factors, values) {
  x.headings <- paste("X", 1:tab$variables, sep = "")
  neg.headings <- paste("N", 1:tab$objectives, sep = "")
  pos.headings <- paste("P", 1:tab$objectives, sep = "")
  x <- rep(0, tab$variables)
  n <- rep(0, tab$objectives)
  p <- rep(0, tab$objectives)
  for (i in 1:tab$objectives) {
    for (j in 1:tab$variables) {
      if (tab$row.headings[i] == x.headings[j]) {
        x[j] <- tab$tb[i]
      }
    }
    for (j in 1:tab$objectives) {
      if (tab$row.headings[i] == neg.headings[j]) {
        n[j] <- tab$tb[i]
      }
    }
    for (j in 1:tab$objectives) {
      if (tab$row.headings[i] == pos.headings[j]) {
        p[j] <- tab$tb[i]
      }
    }
  }
  f <- factors %*% x
  output <- list(x = x, n = n, p = p, f = f, a = tab$ta, b = values)
  class(output) <- "mogpout"
  return(output)
}

dv.mogp <- function (tab, sp) {
  ip <- 0
  for (i in 1:tab$objectives) {
    if (tab$te[i, sp] > 0) {
      v <- tab$tb[i]/tab$te[i, sp]
      if (ip == 0) {
        vip <- v
        ip <- i
      }
      else if (v < vip) {
        vip <- v
        ip <- i
      }
      else if (v == vip) {
        ip <- dv.tie(tab, i, ip)
      }
    }
  }
  return(ip)
}

piv.mogp <- function (tab, nevc, ndvr, verbose) {
  objectives <- tab$objectives
  nonbasics <- tab$nonbasics
  if (verbose) {
    cat("\n")
    cat(paste("Iteration:", tab$iter, ", Entering Variable:", 
              tab$col.headings[nevc], ", Departing Variable:", 
              tab$row.headings[ndvr]))
    cat("\n")
  }
  tab <- swp.headings(tab, ndvr, nevc)
  tab <- swp.vec(tab, ndvr, nevc)
  piv <- tab$te[ndvr, nevc]
  pib <- tab$tb[ndvr]
  for (i in 1:objectives) {
    if (i != ndvr) {
      pix <- tab$te[i, nevc]/piv
      tab$tb[i] <- (tab$tb[i] - pix * pib)
      for (s in 1:nonbasics) {
        if (s != nevc) {
          tab$te[i, s] <- (tab$te[i, s] - tab$te[ndvr, 
                                                 s] * pix)
        }
      }
    }
  }
  for (s in 1:nonbasics) {
    tab$te[ndvr, s] <- (tab$te[ndvr, s]/piv)
  }
  for (i in 1:objectives) {
    tab$te[i, nevc] <- (-tab$te[i, nevc]/piv)
  }
  tab$tb[ndvr] <- tab$tb[ndvr]/piv
  tab$te[ndvr, nevc] <- (1/piv)
  return(tab)
}

mogp <- function (factors, values, conditions, maxiter = 1000, verbose = FALSE){
  tab <- mogptab(factors, values, conditions)
  prnt <- 0
  tab$iter <- 0
  check.tb(tab)
  for (k in 1:tab$levels) {
    tab$level <- k
    tab <- calc.ti.k(tab, k)
    tab <- calc.ta.k(tab, k)
    sp <- ev.mogp(tab, k)
    while (sp != 0) {
      tab$iter <- tab$iter + 1
      if (tab$iter >= maxiter) {
        prnt <- prnt + 1
        cat(paste("Algorithm did not finish", tab$iter, 
                  "iterations at level", k))
        cat("\n")
        print(tab)
        out <- mogpout(tab, factors, values)
        result <- list(tab = tab, out = out, converged = FALSE)
        class(result) <- "mogp"
        return(result)
      }
      ip <- dv.mogp(tab, sp)
      if (ip == 0) {
        cat(paste("Failed pivot computation at level", 
                  k))
        cat("\n")
        prnt <- prnt + 1
        print(tab)
        out <- mogpout(tab, factors, values)
        result <- list(tab = tab, out = out, converged = FALSE)
        class(result) <- "mogp"
        return(result)
      }
      tab <- piv.mogp(tab, sp, ip, verbose)
      tab <- calc.ti.k(tab, k)
      tab <- calc.ta.k(tab, k)
      sp <- ev.mogp(tab, k)
      if (verbose) 
        print(tab)
    }
  }
  out <- mogpout(tab, factors, values)
  result <- list(tab = tab, out = out, converged = TRUE)
  class(result) <- "mogp"
  return(result)
}
factors=matrix(c(1,0,1,1,1,0),3)
values=c(1329*0.7,137*0.8,134*0.4)
conditions=data.frame(objective=1:3,priority=c(1,2,3),p=c(0,0,1),n=c(1,1,0))
soln=mogp(factors,values,conditions) 
soln$converged
soln$out$x[1]

# -------------------------------------------------------------------------

# Arguments
# factors: A coefficient matrix of constraint variables (excluding deviation variables).

# values: A vector of target values corresponding to the coefficient matrix.

# conditions: A data frame with the deviation variables for each objective together with the priority level.
# It consists of four vectors: objective, priority, p and n.
# Objective and priority are positive integers, respectively indicating which pair of deviation variables and the priority level of the deviation variable.
# P and n represent the weight coefficients of d+ (positive deviation variable) and D (negative deviation variable) respectively.

# 1329: The number of targets   
# 137: The intersection between targets and diseases genes   
# 134: Components with a target number above the median