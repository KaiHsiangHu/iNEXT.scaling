TD.m.est = function(x, rho, m, q){
  x = x[x > 0]
  n = sum(x)
  N = ceiling(n / rho)
  
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0 = ifelse(f2 > 0, f1^2 / (n/(n - 1) * 2 * f2 + rho/(1 - rho) * f1), f1*(f1 - 1) / (n/(n - 1) * 2 + rho/(1 - rho) * f1))
  if (is.nan(f0)) f0 = 0
  
  fk.hat <- function(x, m) {
    
    sapply(1:m, function(j) sum(dhyper(j, x, n - x, m)))
  }
  
  D0.hat <- function(x, mie) {
    
    Sub <- function(m) {
      if (m <= n) {
        
        ifi <- table(x)
        ifi <- cbind(i = as.numeric(names(ifi)), fi = ifi)
        RD(ifi, n, m, 0)
        
      } else {
        
        ms = m - n
        N0 = (1/rho - 1) * f1 / f0
        N0 = (N - n + 1) * f1 / (n * f0 + f1)
        if (N0 == "NaN") N0 = 0
        
        length(x) + f0 * (1 - (1 - ms/(N - n))^N0)
      }
    }
    
    int.m = c(floor(mie[mie <= n]), ceiling(mie[mie <= n])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > n] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= n) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
    
  }
  
  D1.hat <- function(x, mie) {
    
    Sub <- function(m) {
      
      if (m <= n) {
        
        ifi <- table(x)
        ifi <- cbind(i = as.numeric(names(ifi)), fi = ifi)
        RD(ifi, n, m, 1)
        
      } else {
        
        N0 = (1/rho - 1) * f1 / f0
        N0 = (N - n + 1) * f1 / (n * f0 + f1)
        if (N0 == "NaN") N0 = 0
        
        ms = m - n
        
        obs + (asy - obs) * (1 - (1 - ms / (N - n))^N0)
      }
    }
    
    asy = Diversity_profile(x, rho, q = 1)
    obs = Diversity_profile_MLE(x, rho, q = 1)
    asy[asy < obs] = obs
    
    int.m = c(floor(mie[mie <= n]), ceiling(mie[mie <= n])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > n] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= n) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
  }
  
  D2.hat <- function(x, mie) {
    
    Sub <- function(m) {
      
      if (m <= n) {
        
        1 / ((1 - m/N) / m + (m + m/N - 1)/m * p2)
        
      } else {
        
        1 / ((1 - m/N) / m + (m + m/N - 1)/m * p2)
      }
    }
    
    p2 = (sum(x * (x - 1)) + n * rho) / (n^2 - n + n * rho)
    
    int.m = c(floor(mie[mie <= n]), ceiling(mie[mie <= n])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > n] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= n) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
  }
  
  Dq.hat <- function(x, mie, q) {
    
    asy = Diversity_profile(x, rho, q)
    obs = Diversity_profile_MLE(x, rho, q)
    
    Sub <- function(m) {
      
      if (m <= n) {
        
        k <- 1:m
        sum((k / m)^q * fk.hat(x, m))^(1 / (1 - q))
        
      } else {
        
        N0 = (1/rho - 1) * f1 / f0
        N0 = (N - n + 1) * f1 / (n * f0 + f1)
        if (N0 == "NaN") N0 = 0
        
        ms = m - n
        
        obs + (asy - obs) * (1 - (1 - ms / (N - n))^N0)
      }
    }
    
    asy = Diversity_profile(x, rho, q = q)
    obs = Diversity_profile_MLE(x, rho, q = q)
    asy[asy < obs] = obs
    
    int.m = c(floor(mie[mie <= n]), ceiling(mie[mie <= n])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > n] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= n) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
  }
  
  iNEXT.func <- function(x, q, mie) {
    
    if (q == 0) 
      
      D0.hat(x, mie)
    
    else if (q == 1) 
      
      D1.hat(x, mie)
    
    else if (q == 2) 
      
      D2.hat(x, mie)
    
    else Dq.hat(x, mie, q)
  }
  
  sapply(q, function(i) iNEXT.func(x, i, m)) %>% as.vector
}


TD.m.est_inc <- function(x, rho, t_, q){
  
  nt = x[1]
  nT = ceiling(nt / rho)
  x = x[-1]
  x = x[x > 0]
  
  Q1 = sum(x == 1)
  Q2 = sum(x == 2)
  Q0 = ifelse(Q2 > 0, Q1^2 / (nt/(nt - 1) * 2 * Q2 + rho/(1 - rho) * Q1), Q1*(Q1 - 1) / (nt/(nt - 1) * 2 + rho/(1 - rho) * Q1))
  if (is.nan(Q0)) Q0 = 0
  
  Qk.hat <- function(x, m) {
    
    sapply(1:m, function(j) sum(dhyper(j, x, nt - x, m)))
  }
  
  D0.hat <- function(x, mie) {
    
    Sub <- function(m) {
      
      if (m <= nt) {
        
        iQi <- table(x)
        iQi <- cbind(i = as.numeric(names(iQi)), Qi = iQi)
        RD_inc(iQi, nt, m, 0)
        
      } else {
        
        ms = m - nt
        M0 = (1/rho - 1) * Q1 / Q0
        M0 = (nT - nt + 1) * Q1 / (nt * Q0 + Q1)
        if (M0 == "NaN") M0 = 0
        
        length(x) + Q0 * (1 - (1 - ms/(nT - nt))^M0)
      }
    }
    
    int.m = c(floor(mie[mie <= nt]), ceiling(mie[mie <= nt])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > nt] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= nt) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
    
  }
  
  D1.hat <- function(x, mie) {
    
    Sub <- function(m) {
      
      if (m <= nt) {
        
        iQi <- table(x)
        iQi <- cbind(i = as.numeric(names(iQi)), Qi = iQi)
        RD_inc(iQi, nt, m, 1)
        
      } else {
        
        M0 = (1/rho - 1) * Q1 / Q0
        M0 = (nT - nt + 1) * Q1 / (nt * Q0 + Q1)
        if (M0 == "NaN") M0 = 0
        
        ms = m - nt
        
        obs + (asy - obs) * (1 - (1 - ms / (nT - nt))^M0)
      }
    }
    
    asy = Diversity_profile.inc(c(nt, x), rho, q = 1)
    obs = Diversity_profile_MLE.inc(c(nt, x), rho, q = 1)
    asy[asy < obs] = obs
    
    int.m = c(floor(mie[mie <= nt]), ceiling(mie[mie <= nt])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > nt] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= nt) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
    
  }
  
  D2.hat <- function(x, mie) {
    
    Sub <- function(m) {
      
      if (m <= nt) {
        
        (m * sum(x) / nt)^2 / ((1 - m/nT) * m * sum(x) / nt  + (m^2 - (1 - m / nT) * m) * p2)
        
      } else {
        
        (m * sum(x) / nt)^2 / ((1 - m/nT) * m * sum(x) / nt  + (m^2 - (1 - m / nT) * m) * p2)
      }
    }
    
    p2 = (sum(x * (x - 1)) + sum(x) * rho) / (nt^2 - nt + nt * rho)
    
    int.m = c(floor(mie[mie <= nt]), ceiling(mie[mie <= nt])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > nt] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= nt) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
    
  }
  
  Dq.hat <- function(x, mie, q) {
    
    Sub <- function(m) {
      
      if (m <= nt) {
        
        Ut = sum(x) * m / nt
        k <- 1:m
        sum((k/Ut)^q * Qk.hat(x, m))^(1 / (1 - q))
        
      } else {
        
        M0 = (1/rho - 1) * Q1 / Q0
        M0 = (nT - nt + 1) * Q1 / (nt * Q0 + Q1)
        if (M0 == "NaN") M0 = 0
        
        ms = m - nt
        
        obs + (asy - obs) * (1 - (1 - ms / (nT - nt))^M0)
      }
    }
    
    asy = Diversity_profile.inc(c(nt, x), rho, q = q)
    obs = Diversity_profile_MLE.inc(c(nt, x), rho, q = q)
    asy[asy < obs] = obs
    
    int.m = c(floor(mie[mie <= nt]), ceiling(mie[mie <= nt])) %>% unique %>% sort
    
    mRTD = sapply(int.m, function(m) if (m == 0) 0 else Sub(m))
    
    ext.m = mie[mie > nt] %>% unique
    if (length(ext.m) != 0) mETD = sapply(ext.m, function(m) Sub(m))
    
    sapply(mie, function(m) {
      
      if (m <= nt) {
        
        if (m == round(m)) mRTD[int.m == m] else 
          (ceiling(m) - m) * mRTD[int.m == floor(m)] + (m - floor(m)) * mRTD[int.m == ceiling(m)]
        
      } else mETD[ext.m == m]
      
    })
    
  }
  
  iNEXT.func <- function(x, q, mie) {
    
    if (q == 0) 
      
      D0.hat(x, mie)
    
    else if (q == 1) 
      
      D1.hat(x, mie)
    
    else if (q == 2) 
      
      D2.hat(x, mie)
    
    else Dq.hat(x, mie, q)
  }
  
  sapply(q, function(i) iNEXT.func(x, i, t_)) %>% as.vector
}


iNEXT.Ind <- function(Spec, rho, q=0, m=NULL, knots=40, nboot=200, conf=0.95){
  qtile <- qnorm(1 - (1 - conf)/2)
  n <- sum(Spec)
  
  Dq.hat <- TD.m.est(Spec, rho, m, q)
  C.hat <- Coverage(Spec, rho, 'abundance', m)
  
  goalSC <- unique(C.hat)
  Dq.hat_unc <- unique(invChat.Ind(x = Spec, rho, q = q, C = goalSC))
  
  refC <- Coverage(Spec, rho, 'abundance', n)
  Dq.hat_unc$Method[Dq.hat_unc$SC == refC] = "Observed"
  
  if(nboot > 1 & length(Spec) > 1) {
    
    Abun.Mat <- EstiBootComm.Ind(Spec, rho, nboot)
    
    ses_m <- apply(matrix(apply(Abun.Mat,2 ,function(x) TD.m.est(x, rho, m, q)), nrow = length(Dq.hat)), 
                   1, sd, na.rm = TRUE)
    
    ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(x) Coverage(x, rho, 'abundance', m)), nrow = length(m)),
                        1, sd, na.rm = TRUE)
    
    ses_C <- apply(matrix(apply(Abun.Mat, 2, function(x) invChat.Ind(x, rho, q, unique(Dq.hat_unc$SC))$qTD), nrow = nrow(Dq.hat_unc)), 
                   1, sd, na.rm = TRUE)
    
  } else {
    
    ses_m <- rep(NA, length(Dq.hat))
    
    ses_C_on_m <- rep(NA, length(m))
    
    ses_C <- rep(NA, nrow(Dq.hat_unc))
  }
  
  out_m <- data.frame(m = rep(m, length(q)), 
                      qTD = Dq.hat, 
                      qTD.LCL = Dq.hat - qtile * ses_m,
                      qTD.UCL = Dq.hat + qtile * ses_m,
                      SC = rep(C.hat, length(q)), 
                      SC.LCL = C.hat - qtile * ses_C_on_m, 
                      SC.UCL = C.hat + qtile * ses_C_on_m)
  
  out_m$Method <- ifelse(out_m$m < n, "Rarefaction", ifelse(out_m$m == n, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(m))
  
  id_m <- match(c("Order.q", "m", "Method", "qTD", "qTD.LCL", "qTD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  
  out_m <- out_m[, id_m]
  out_m$qTD.LCL[out_m$qTD.LCL < 0] <- 0
  out_m$SC.LCL[out_m$SC.LCL < 0] <- 0
  out_m$SC.UCL[out_m$SC.UCL > 1] <- 1
  
  
  out_C <- data.frame(Dq.hat_unc, 
                      qTD.LCL = Dq.hat_unc$qTD - qtile * ses_C,
                      qTD.UCL = Dq.hat_unc$qTD + qtile * ses_C) 
  id_C <- match(c("Order.q", "SC", "m", "Method", "qTD", "qTD.LCL", "qTD.UCL"), names(out_C), nomatch = 0)
  
  out_C <- out_C[, id_C]
  out_C$qTD.LCL[out_C$qTD.LCL < 0] <- 0
  
  return(list(size_based = out_m, coverage_based = out_C))
  }


iNEXT.Sam <- function(Spec, rho, t=NULL, q=0, knots=40, nboot=200, conf=0.95){
  qtile <- qnorm(1 - (1 - conf) / 2)
  
  nT <- ncol(Spec)
  
  Dq.hat <- TD.m.est_inc(as.incfreq(Spec), rho, t, q)
  C.hat <- Coverage(Spec, rho, "incidence_raw", t)
  
  goalSC <- unique(C.hat)
  Dq.hat_unc <- unique(invChat.Sam(x = as.incfreq(Spec), rho, q = q, C = goalSC))
  
  refC <- Coverage(Spec, rho, "incidence_raw", nT)
  Dq.hat_unc$Method[Dq.hat_unc$SC == refC] = "Observed"
  
  
  if(nboot > 1 & length(Spec) > 2) {
    
    Abun.Mat <- EstiBootComm.Sam(Spec, rho, nboot)
    
    tmp <- which(colSums(Abun.Mat) == nT)
    
    if (length(tmp) > 0) Abun.Mat <- Abun.Mat[, -tmp]
    
    if (ncol(Abun.Mat) == 0) {
      
      out <- cbind("t" = t, "qTD" = Dq.hat, "SC" = C.hat)
      
      warning("Insufficient data to compute bootstrap s.e.")
      
    } else {		
      
      ses_m <- apply(matrix(apply(Abun.Mat, 2, function(y) TD.m.est_inc(y, rho, t, q)), nrow = length(Dq.hat)),
                     1, sd, na.rm = TRUE)
      
      ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(y) Coverage(y, rho, "incidence_freq", t)),nrow = length(t)),
                          1, sd, na.rm = TRUE)
      
      ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(y) invChat.Sam(y, rho, q,unique(Dq.hat_unc$SC))$qTD), nrow = nrow(Dq.hat_unc)),
                     1, sd, na.rm = TRUE)
      
    }
    
  } else {
    
    ses_m <- rep(NA, length(Dq.hat))
    
    ses_C_on_m <- rep(NA, length(t))
    
    ses_C <- rep(NA, nrow(Dq.hat_unc))
    
  }
  
  out_m <- data.frame(nT = rep(t,length(q)), 
                      qTD = Dq.hat, 
                      qTD.LCL = Dq.hat - qtile * ses_m,
                      qTD.UCL = Dq.hat + qtile * ses_m,
                      SC = rep(C.hat, length(q)), 
                      SC.LCL = C.hat - qtile * ses_C_on_m, 
                      SC.UCL = C.hat + qtile * ses_C_on_m)
  
  out_m$Method <- ifelse(out_m$nT < nT, "Rarefaction", ifelse(out_m$nT == nT, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(t))
  
  id_m <- match(c("Order.q", "nT", "Method", "qTD", "qTD.LCL", "qTD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  
  out_m <- out_m[, id_m]
  out_m$qTD.LCL[out_m$qTD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  
  out_C <- data.frame(Dq.hat_unc,
                      qTD.LCL = Dq.hat_unc$qTD - qtile * ses_C,
                      qTD.UCL = Dq.hat_unc$qTD + qtile * ses_C) 
  id_C <- match(c("Order.q", "SC", "nT", "Method", "qTD", "qTD.LCL", "qTD.UCL"), names(out_C), nomatch = 0)
  
  out_C <- out_C[, id_C]
  out_C$qTD.LCL[out_C$qTD.LCL < 0] <- 0
  
  return(list(size_based = out_m, coverage_based = out_C))
  }


invChat <- function(x, rho, q, datatype = "abundance", C = NULL,nboot = 0, conf = NULL) {
  
  qtile <- qnorm(1 - (1 - conf) / 2)
  
  if (datatype == "abundance") {
    
    out <- lapply(1:length(x), function(i) {
      
      est <- invChat.Ind(x[[i]], rho[i], q, C)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Ind(x[[i]], rho[i], nboot)
        
        ses <- apply(matrix(apply(Abun.Mat, 2, function(a) invChat.Ind(a, rho[i], q,C)$qTD), nrow = length(q) * length(C)), 1, sd)
        
      } else {
        
        ses <- rep(0,nrow(est))
      }
      
      cbind(est, s.e. = ses, qTD.LCL = est$qTD - qtile * ses, qTD.UCL = est$qTD + qtile * ses)
    })
    
    out <- do.call(rbind,out)
    out$Assemblage <- rep(names(x), each = length(q) * length(C))
    
    out <- out[, c(ncol(out), seq(1, (ncol(out) - 4)), (ncol(out) - 2),(ncol(out) - 1),(ncol(out) - 3))]
    rownames(out) <- NULL
    
    out = out %>% select(c('Assemblage', 'Order.q', 'SC', 'm', 'Method', 'qTD', 's.e.', 'qTD.LCL', 'qTD.UCL'))
    
  } else if (datatype == "incidence_raw") {
    
    out <- lapply(1:length(x), function(i) {
      
      est <- invChat.Sam(as.incfreq(x[[i]]), rho[i], q, C)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Sam(x[[i]], rho[i], nboot)
        
        tmp <- which(colSums(Abun.Mat) == ncol(x[[i]]))
        
        if(length(tmp) > 0) Abun.Mat <- Abun.Mat[,-tmp]
        
        if (ncol(Abun.Mat) == 0) warning("Insufficient data to compute bootstrap s.e.")
        
        ses <- apply(matrix(apply(Abun.Mat, 2, function(a) invChat.Sam(a, rho[i], q, C)$qTD), nrow = length(q)* length(C)), 1, sd)
        
      } else {
        
        ses <- rep(0,nrow(est))
      }
      
      cbind(est, s.e. = ses, qTD.LCL = est$qTD - qtile * ses, qTD.UCL = est$qTD + qtile * ses)
    })
    
    out <- do.call(rbind,out)
    out$Assemblage <- rep(names(x),each = length(q)*length(C))
    
    out <- out[, c(ncol(out), seq(1, (ncol(out) - 4)), (ncol(out) - 2), (ncol(out) - 1), (ncol(out) - 3))]
    rownames(out) <- NULL
    
    out = out %>% select(c('Assemblage', 'Order.q', 'SC', 'nT', 'Method', 'qTD', 's.e.', 'qTD.LCL', 'qTD.UCL'))
    
  }
  
  out$qTD.LCL[out$qTD.LCL < 0] <- 0
  out
}


invChat.Ind <- function (x, rho, q, C) {
  
  n = sum(x)
  N = ceiling(n / rho)
  x = x[x > 0]
  refC = Coverage(x, rho, "abundance", sum(x))
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  
  f <- function(m, C) abs(Coverage(x, rho, "abundance", m) - C)
  
  mm <- sapply(C, function(cvrg) {
    
    if (refC == cvrg) {
      
      n
      
    } else if (refC > cvrg) {
      
      opt <- optimize(f, C = cvrg, lower = 0, upper = sum(x))
      
      opt$minimum
      
    } else if (refC < cvrg) {
      
      N1 = (2 * (N - n + 2) * f2 + (n - 1) * f1) / ((n - 1) * f1 + 2 * f2)
      if (N1 == "NaN" | N1 == Inf) N1 = 0
      
      if (N1 == 0 | f1 == 0) {
        
        ms = 0
        
      } else ms <- (N - n) * (1 - exp(log((1 - cvrg) / (1 - rho) * n / f1) / N1) )
      
      n + ms
    }
  })
  
  mm[mm > N] = N
  mm
  
  out <- TD.m.est(x = x, rho, m = mm, q = q)
  method <- ifelse(mm > n, 'Extrapolation', ifelse(mm < n, 'Rarefaction', 'Observed'))
  method <- rep(method, length(q))
  
  m <- rep(mm, length(q))
  order <- rep(q, each = length(mm))
  SC <- rep(C, length(q))
  
  data.frame(m = m,
             Method = method,
             Order.q = order,
             SC = SC,
             qTD = out)
}


invChat.Sam <- function (x, rho, q, C) {
  
  nt = x[1]
  nT = ceiling(nt / rho)
  
  x = x[-1]
  x = x[x > 0]
  refC = Coverage(c(nt, x), rho, "incidence_freq", nt)
  Q1 = sum(x == 1)
  Q2 = sum(x == 2)
  
  
  f <- function(m, C) abs(Coverage(c(nt, x), rho, "incidence_freq", m) - C)
  
  mm <- sapply(C, function(cvrg) {
    
    if (refC == cvrg) {
      
      nt
      
    } else if (refC > cvrg) {
      
      opt <- optimize(f, C = cvrg, lower = 0, upper = nt)
      
      opt$minimum
      
    } else if (refC < cvrg) {
      
      M1 = (2 * (nT - nt + 2) * Q2 + (nt - 1) * Q1) / ((nt - 1) * Q1 + 2 * Q2)
      if (M1 == "NaN" | M1 == Inf) M1 = 0
      
      
      if (M1 == 0 | Q1 == 0) {
        
        ms = 0
        
      } else ms <- (nT - nt) * (1 - exp(log((1 - cvrg) / (1 - rho) * sum(x) / Q1) / M1) )
      
      nt + ms
    }
  })
  
  mm[mm > nT] = nT
  
  out <- TD.m.est_inc(c(nt, x), rho, t_ = mm, q = q)
  
  method <- ifelse(mm > nt, 'Extrapolation', ifelse(mm < nt, 'Rarefaction', 'Observed'))
  method <- rep(method, length(q))
  m <- rep(mm, length(q))
  
  order <- rep(q, each = length(mm))
  SC <- rep(C, length(q))
  
  data.frame(nT = m,
             Method = method,
             Order.q = order,
             SC = SC,
             qTD = out)
  
}


invSize <- function(x, rho, q, datatype="abundance", size=NULL, nboot=0, conf=NULL){
  
  qtile <- qnorm(1 - (1 - conf) / 2)
  
  if (datatype == "abundance") {
    
    out <- lapply(1:length(x), function(i){
      
      est <- invSize.Ind(x[[i]], rho[i], q, size)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Ind(x[[i]], rho[i], nboot)
        
        ses <- apply(matrix(apply(Abun.Mat, 2, function(a) invSize.Ind(a, rho[i], q,size)$qTD), nrow = length(q) * length(size)), 1, sd)
        
      } else {
        
        ses <- rep(0, nrow(est))
      }
      
      cbind(est,
            s.e. = ses,
            qTD.LCL = est$qTD - qtile * ses,
            qTD.UCL = est$qTD + qtile * ses)
      
    })
    
    out <- do.call(rbind,out)
    out$Assemblage <- rep(names(x), each = length(q) * length(size))
    
    out <- out[, c(ncol(out), seq(1, (ncol(out) - 1)))]
    rownames(out) <- NULL
    
  } else if (datatype == "incidence_raw") {
    
    out <- lapply(1:length(x), function(i){
      
      est <- invSize.Sam(as.incfreq(x[[i]]), rho[i], q, size)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Sam(x[[i]], rho[i], nboot)
        
        tmp <- which(colSums(Abun.Mat) == ncol(x[[i]]))
        
        if (length(tmp) > 0) Abun.Mat <- Abun.Mat[, -tmp]
        
        if (ncol(Abun.Mat) == 0) warning("Insufficient data to compute bootstrap s.e.")
        
        ses <- apply(matrix(apply(Abun.Mat, 2, function(a) invSize.Sam(a, rho[i], q, size)$qTD), nrow = length(q) * length(size)), 1, sd)
        
      } else {
        
        ses <- rep(0,nrow(est))
      }
      
      cbind(est,
            s.e. = ses,
            qTD.LCL = est$qTD - qtile*ses,
            qTD.UCL = est$qTD + qtile*ses)
    })
    
    out <- do.call(rbind,out)
    out$Assemblage <- rep(names(x), each = length(q) * length(size))
    
    out <- out[, c(ncol(out), seq(1, (ncol(out) - 1)))]
    rownames(out) <- NULL
    
  }
  
  out
}


invSize.Ind <- function(x, rho, q, size){
  
  n <- sum(x)
  out <- TD.m.est(x, rho, m = size, q = q)
  
  SC <- Coverage(x, rho, 'abundance', size)
  
  method <- ifelse(size > n, 'Extrapolation', ifelse(size < n, 'Rarefaction', 'Observed'))
  method <- rep(method, length(q))
  
  m <- rep(size, length(q))
  order <- rep(q, each = length(size))
  SC <- rep(SC,length(q))
  
  data.frame(Order.q = order,
             m = m,
             Method = method,
             SC = SC,
             qTD = out)
  }


invSize.Sam <- function(x, rho, q, size){
  
  n <- ncol(x)
  out <- TD.m.est_inc(x, rho, t_ = size, q)
  
  SC <- Coverage(x, rho,"incidence_freq", size)
  
  method <- ifelse(size > n, 'Extrapolation', ifelse(size < n, 'Rarefaction', 'Observed'))
  method <- rep(method, length(q))
  
  m <- rep(size, length(q))
  order <- rep(q, each = length(size))
  SC <- rep(SC, length(q))
  
  data.frame(Order.q = order,
             nT = m,
             Method = method,
             SC = SC,
             qTD = out)
  }



Diversity_profile <- function(x, rho, q){
  
  x = x[x > 0]
  n = sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  N = ceiling(n / rho)
  p1 = ifelse(f2 > 0 | f1 > 0, ((1 - rho) * 2 * f2 + rho * f1) / ((n-1) * f1 + 2 * f2), 0)
  
  Hill <- function(q) {
    
    if (q == 0) {
      
      lbd = ifelse(f2 > 0, f1^2 / (n/(n - 1) * 2 * f2 + rho/(1 - rho) * f1), f1 * (f1 - 1) / (n/(n - 1) * 2 + rho/(1 - rho) * f1))
      sum(x > 0) + ifelse(is.nan(lbd), 0, lbd)
      
    } else if (q == 1) {
      
      A = (1 - rho) * sum(tab * sortx / n * (digamma(n) - digamma(sortx)))
      
      B <- D1_2nd(n, f1, f2, rho)
      
      if (B == "NaN" | B == Inf) B = 0
      
      MLEpart = -sum((x / n) * log(x / n))
      
      exp(A + B + rho * MLEpart)
      
    } else if (q == 2) {
      
      1 / ((sum(x * (x - 1)) + n * rho) / (n^2 - n + n * rho) )
      
    } else {
      
      AB = ans[which(q_part2 == q)]
      
      MLEpart = sum( (x / n)^q )
      
      (AB + rho * MLEpart)^(1 / (1 - q))
      
    }
  }
  
  sortx = sort(unique(x))
  tab = table(x)
  
  q_part2 <- q[!q %in% c(0, 1, 2)]
  if (length(q_part2) > 0) ans <- Dq(ifi = cbind(i = sortx, fi = tab), n = n, qs = q_part2, f1 = f1, A = p1, rho = rho)
  
  sapply(q, Hill)
}


Diversity_profile.inc <- function(x, rho, q){
  
  x = x[x > 0]
  nt = x[1]
  x = x[-1]
  
  Q1 = sum(x == 1)
  Q2 = sum(x == 2)
  nT = ceiling(nt / rho)
  p1 = ifelse(Q2 > 0 | Q1 > 0, ((1 - rho) * 2 * Q2 + rho * Q1)/((nt - 1) * Q1 + 2 * Q2), 0)
  
  Hill <- function(q) {
    
    if(q == 0){
      
      lbd = ifelse(Q2 > 0, Q1^2 / (nT/(nT - 1) * 2 * Q2 + rho/(1 - rho) * Q1), Q1 * (Q1 - 1) / (nT/(nT - 1) * 2 + rho/(1 - rho) * Q1))
      sum(x > 0) + ifelse(is.nan(lbd), 0, lbd)
      
    } else if(q == 1) {
      
      r <- 1:(nt - 1)
      
      A <- sum(sapply(r, function(k) {
        
        1 / k * (1 - rho) * sum(tab * sortx / nt * exp(lchoose(nt - sortx, k) - lchoose(nt - 1, k)))
        
      }))
      
      B <- D1_2nd(nt, Q1, Q2, rho)
      
      if (B == "NaN" | B == Inf) B = 0
      
      
      MLEpart = -sum((x / nt) * log(x / nt))
      
      exp(nt / sum(x) * (A + B + rho * MLEpart) + log(sum(x) / nt))
      
    } else if (q == 2) {
      
      (sum(x) / nt)^2 / ((sum(x * (x - 1)) + sum(x) * rho) / (nt * (nt - 1) + nt * rho))
      
    } else {
      
      AB = ans[which(q_part2 == q)]
      
      MLEpart = sum( (x / nt)^q )
      
      (sum(x) / nt)^(q / (q - 1)) * (AB + rho * MLEpart)^(1 / (1 - q))
      
    }
  }
  
  sortx = sort(unique(x))
  tab = table(x)
  
  q_part2 <- q[!q %in% c(0, 1, 2)]
  if (length(q_part2) > 0) ans <- Dq(ifi = cbind(i = sortx, fi = tab), n = nt, qs = q_part2, f1 = Q1, A = p1, rho = rho)
  
  sapply(q, Hill)
  }


Diversity_profile_MLE <- function(x, rho, q){
  x <- x[x > 0]
  x = x / rho
  N = sum(x)
  
  Hill <- function(q){
    
    if (q == 1) {
      
      exp(-sum(x/N * log(x/N)))
      
    } else {
      
      sum((x / N)^q)^(1 / (1 - q))
      
    }
  }
  
  sapply(q, Hill)
  }


Diversity_profile_MLE.inc <- function(x, rho, q){
  x = x[-1]
  x <- x[x > 0] / rho
  
  Hill <- function(q){
    
    if (q == 1) {
      
      exp(-sum(x / sum(x) * log(x / sum(x))))
      
    } else {
      
      sum((x / sum(x))^q)^(1 / (1 - q))
      
    }
  }
  
  sapply(q, Hill)
}


EstiBootComm.Ind <- function(x, rho, B){
  
  x <- x[x != 0]
  n <- sum(x)
  N = ceiling(n / rho)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0 = ceiling( ifelse(f2 > 0, f1^2 / (n/(n - 1) * 2 * f2 + rho/(1 - rho) * f1), f1 * (f1 - 1) / (n/(n - 1) * 2 + rho/(1 - rho) * f1)) )
  if (f0 == "NaN") f0 = 0
  
  Chat = 1 - (1 - rho) * f1/n
  
  lamda_hat = (1 - Chat) / sum((x / n) * (1 - rho)^(x / rho)) 
  
  if (lamda_hat == "NaN") lamda_hat = 0
  
  Ni_det = (x / rho)*(1 - lamda_hat * (1 - rho)^(x / rho)) 
  Ni_undet = N * (1 - Chat) / f0
  
  if (Ni_undet == "NaN") Ni_undet = 0 else if (Ni_undet < 1 & Ni_undet > 0) Ni_undet = 1
  
  N_hat = round( c(Ni_det, rep(Ni_undet, f0)), 0)
  
  ex.sample = unlist(lapply(1:length(N_hat), function(i) rep(i, N_hat[i])))
  sample_result <- sapply(1:B, function(i) sample(ex.sample, n, replace = FALSE))
  random = lapply(1:B, function(j) as.numeric(table(sample_result[,j])))
  
  random <- do.call(cbind.data.frame,
                    lapply(lapply(random, unlist), `length<-`, max(lengths(random))))
  
  random[is.na(random)] <- 0
  colnames(random) = NULL
  
  return(random)
}


EstiBootComm.Sam <- function(x, rho, B){
  
  if (rho != 1) {
    
    raw = x[rowSums(x) > 0,]
    nt <- ncol(x)
    nT = ceiling(nt / rho)
    x = rowSums(x)
    x <- x[x != 0]
    
    Q1 = sum(x == 1)
    Q2 = sum(x == 2)
    Q0 = ceiling( ifelse(Q2 > 0, Q1^2 / (nt/(nt - 1) * 2 * Q2 + rho/(1 - rho) * Q1), Q1 * (Q1 - 1) / (nt/(nt - 1) * 2 + rho/(1 - rho) * Q1)) )
    if (Q0 == "NaN") Q0 = 0
    
    Chat = 1 - (1 - rho) * Q1/sum(x)
    
    lamda_hat = (1 - Chat) / sum((x / nt) * (1 - rho)^(x / rho)) 
    
    if (lamda_hat == "NaN") lamda_hat = 0
    
    Ni_det = (x / rho)*(1 - lamda_hat * (1 - rho)^(x / rho)) 
    Ni_undet = nT * (1 - Chat) / Q0
    
    if (Ni_undet == "NaN") Ni_undet = 0 else if (Ni_undet < 1 & Ni_undet > 0) Ni_undet = 1
    
    N_hat = round( c(Ni_det, rep(Ni_undet, Q0)), 0)
    
    random = sapply(1:B, function(j) {
      
      raw.undet = matrix(0, nrow = Q0, ncol = ncol(raw))
      N_hat[1:length(x)] = N_hat[1:length(x)] - x
      
      quad.undet = sapply(N_hat, function(k) {
        candidate = rep(0, nT - nt)
        candidate[sample(1:(nT - nt), k, replace = FALSE)] = 1
        candidate
      }) %>% t
      
      boot.popu = cbind(rbind(raw, raw.undet), quad.undet)
      c(nt, boot.popu[, sample(1:ncol(boot.popu), nt, replace = FALSE)] %>% rowSums)
    })
    
  } else if (rho == 1) random = sapply(1:B, function(j) c(ncol(x), rowSums(x)))
  
  
  return(random)
}



asyTD = function(data, rho, datatype, q, nboot, conf) {
  
  if (datatype == "abundance") {
    
    out <- lapply(1:length(data), function(i) {
      
      dq <- Diversity_profile(data[[i]], rho[i],q)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Ind(data[[i]], rho[i], nboot)
        
        mt = apply(Abun.Mat, 2, function(xb) Diversity_profile(xb, rho[i], q))
        
        if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
        
        error <- qnorm(1 - (1 - conf) / 2) * apply(mt, 1, sd, na.rm = TRUE)
        
      } else error = NA
      
      out <- data.frame(Assemblage = names(data)[i], 
                        Order.q = q, 
                        qTD = dq, 
                        s.e. = error / qnorm(1 - (1 - conf) / 2),
                        qTD.LCL = dq - error, 
                        qTD.UCL = dq + error, 
                        Method = "Asymptotic")
      
      out$qTD.LCL[out$qTD.LCL < 0] <- 0
      
      out
    })
    
    out <- do.call(rbind, out)
    
  } else if (datatype == "incidence_raw") {
    
    out <- lapply(1:length(data), function(i) {
      
      dq <- Diversity_profile.inc(as.incfreq(data[[i]]), rho[i], q)
      
      names(dq) = NULL
      
      if(nboot > 1){
        
        Abun.Mat <- EstiBootComm.Sam(data[[i]], rho[i], nboot)
        
        tmp <- which(colSums(Abun.Mat) == ncol(data[[i]]))
        
        if (length(tmp) > 0) Abun.Mat <- Abun.Mat[, -tmp]
        
        if (ncol(Abun.Mat) == 0) {
          
          error = 0
          warning("Insufficient data to compute bootstrap s.e.")
          
        } else {
          
          mt = apply(Abun.Mat, 2, function(yb) Diversity_profile.inc(yb, rho[i], q))
          
          if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
          
          error <- qnorm(1-(1-conf)/2) * apply(mt, 1, sd, na.rm = TRUE)
        }
        
      } else error = NA
      
      out <- data.frame(Assemblage = names(data)[i], 
                        Order.q = q, 
                        qTD = dq, 
                        s.e. = error / qnorm(1 - (1 - conf) / 2),
                        qTD.LCL = dq - error, 
                        qTD.UCL = dq + error, 
                        Method = "Asymptotic")
      
      out$qTD.LCL[out$qTD.LCL < 0] <- 0
      
      out
    })
    
    out <- do.call(rbind,out)
  }
  
  return(out)
}



obsTD = function(data, rho, datatype, q, nboot, conf) {
 
  if (datatype == "abundance") {
    
    out <- lapply(1:length(data), function(i) {
      
      dq <- Diversity_profile_MLE(data[[i]], rho[i],q)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Ind(data[[i]], rho[i], nboot)
        
        mt = apply(Abun.Mat, 2, function(xb) Diversity_profile_MLE(xb, rho[i], q))
        
        if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
        
        error <- qnorm(1 - (1 - conf) / 2) * apply(mt, 1, sd, na.rm = TRUE)
        
      } else error = NA
      
      out <- data.frame(Assemblage = names(data)[i], 
                        Order.q = q, 
                        qTD = dq, 
                        s.e. = error / qnorm(1 - (1 - conf) / 2),
                        qTD.LCL = dq - error, 
                        qTD.UCL = dq + error, 
                        Method = "Observed")
      
      out$qTD.LCL[out$qTD.LCL < 0] <- 0
      
      out
    })
    
    out <- do.call(rbind,out)
    
  } else if (datatype == "incidence_raw") {
    
    out <- lapply(1:length(data),function(i){
      
      dq <- Diversity_profile_MLE.inc(as.incfreq(data[[i]]), rho[i], q)
      
      if (nboot > 1) {
        
        Abun.Mat <- EstiBootComm.Sam(data[[i]], rho[i], nboot)
        
        tmp <- which(colSums(Abun.Mat) == ncol(data[[i]]))
        
        if (length(tmp) > 0) Abun.Mat <- Abun.Mat[, -tmp]
        
        if (ncol(Abun.Mat) == 0) {
          
          error = 0
          
          warning("Insufficient data to compute bootstrap s.e.")
          
        } else {	
          
          mt = apply(Abun.Mat, 2, function(yb) Diversity_profile_MLE.inc(yb, rho[i], q))
          
          if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
          
          error <- qnorm(1 - (1 - conf) / 2) * apply(mt, 1, sd, na.rm = TRUE)
        }
        
      } else error = NA
      
      out <- data.frame(Assemblage = names(data)[i],
                        Order.q = q, 
                        qTD = dq, 
                        s.e. = error / qnorm(1 - (1 - conf) / 2),
                        qTD.LCL = dq - error, 
                        qTD.UCL = dq + error, 
                        Method = "Observed")
      
      out$qTD.LCL[out$qTD.LCL < 0] <- 0
      
      out
    })
    
    out <- do.call(rbind,out)
  }
  
  return(out)
}



TDinfo = function(data, rho, datatype) {
  
  Fun.abun <- function(x, rho){
    
    n <- sum(x)
    N = n / rho
    
    fk <- sapply(1:5, function(k) sum(x == k))
    Sobs <- sum(x > 0)
    
    Chat <- Coverage(x, rho, datatype, sum(x))
    Chat2n <- Coverage(x, rho, datatype, 2 * sum(x))
    
    c(n, round(N), rho, Sobs, Chat, Chat2n, fk)
  }
  
  Fun.inci <- function(x, rho) {
    
    nT <- x[1]
    x <- x[-1]
    Ts = nT / rho
    U <- sum(x)
    
    Qk <- sapply(1:5, function(k) sum(x == k))
    Sobs <- sum(x > 0)
    
    Chat <- Coverage(c(nT, x), rho, datatype, nT)
    Chat2T <- Coverage(c(nT, x), rho, datatype, 2 * nT)
    
    c(nT, Ts, rho, U, Sobs, Chat, Chat2T, Qk)
  }
  
  if (datatype == "abundance") {
    
    out <- do.call("rbind", lapply(1:length(data), function(i) Fun.abun(data[[i]], rho[i])))
    
    out <- data.frame(site = names(data), out)
    
    colnames(out) <- c("Assemblage", "n", "N", "rho", "S.obs", "SC(n)", "SC(2n)", paste("f", 1:5, sep = ""))
    rownames(out) <- NULL
    
  } else if (datatype == "incidence_raw") {
    
    out <- do.call("rbind", lapply(1:length(data), function(i) Fun.inci(as.incfreq(data[[i]]), rho[i])))
    
    out <- data.frame(site = names(data), out)
    
    colnames(out) <- c("Assemblage","T", "T*", "rho", "U", "S.obs", "SC(T)", "SC(2T)", paste("Q", 1:5, sep = ""))
    rownames(out) <- NULL
    
  }
  
  out
  }



