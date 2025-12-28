as.incfreq <- function(data, nT = NULL) {
  
  if (inherits(data, c("data.frame", "matrix"))) {
    
    if (sum(data > 1) != 0) stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
    
    if (is.null(nT)) nT = ncol(data)
    
    if (inherits(nT, 'data.frame')) nT = unlist(nT)
    mydata = list()
    
    if (ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
    
    ntmp <- 0
    
    for (i in 1:length(nT)) {
      
      mydata[[i]] <- data[, (ntmp + 1):(ntmp + nT[i]), drop=FALSE]
      ntmp <- ntmp+nT[i]
    }
    
    if (is.null(names(nT))) {
      
      names(mydata) <- paste0("Assemblage",1:length(nT))
      
    } else {
      
      names(mydata) = names(nT)
    }
    
    data = lapply(mydata, function(i) {
      
      out = c('nT' = ncol(i), rowSums(i))
      return(out)
    })
    
  } else if (inherits(data, "list")) 
    
    data <- lapply(data, function(i) {
      
      if (sum(i > 1) != 0) stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
      c('nT' = ncol(i), rowSums(i))
      
      }) else if (inherits(data, "integer")) data <- list("Assemblage_1" = c('nT' = length(data), sum(data)))
        
  if (length(data) == 1) data = data[[1]]
  
  return(data)
}

Coverage = function(x, rho, datatype, m){
  
  if (!(datatype %in% c('abundance', 'incidence_freq', 'incidence_raw'))) stop("Invalid Coverage datatype", call. = FALSE)
  
  if (datatype == 'incidence_raw') {x = as.incfreq(x); datatype = 'incidence_freq'}
  
  x <- x[x>0]
  
  Sub <- function(m) {
    n = sum(x)
    N = ceiling(n / rho)
    m[m > N] = N
    
    x = x[x > 0]
    f1 = sum(x == 1)
    f2 = sum(x == 2)
    
    sapply(m, function(k) {
      
      if (k < n) {
        
        Chat = 1 - (N - k)/N * sum(x/n * exp(lgamma(n - x + 1) - lgamma(n - x - k + 1) - lgamma(n) + lgamma(n - k)))
        
      } else {
        
        ms = k - n
        N1 = (2 * (N - n + 2) * f2 + (n - 1) * f1) / ((n - 1) * f1 + 2 * f2)
        if (N1 == "NaN" | N1 == Inf) N1 = 0
        
        Chat = 1 - (1 - rho) * f1/n * (1 - ms / (N - n))^N1
        if (rho == 1) Chat = 1
        if (k >= N) Chat = 1
        
      }
      
      Chat
    })
  }
  
  Sub2 <- function(m) {
    nt = x[1]
    nT = ceiling(nt / rho)
    m[m > nT] = nT
    
    x = x[-1]
    x = x[x > 0]
    Q1 = sum(x == 1)
    Q2 = sum(x == 2)
    
    sapply(m, function(k) {
      
      if (k < nt) {
        
        Chat = 1 - (nT - k)/nT * sum(x/sum(x) * exp(lgamma(nt - x + 1) - lgamma(nt - x - k + 1) - lgamma(nt) + lgamma(nt - k)))
        
      } else {
        
        ms = k - nt
        M1 = (2 * (nT - nt + 2) * Q2 + (nt - 1) * Q1) / ((nt - 1) * Q1 + 2 * Q2)
        if (M1 == "NaN" | M1 == Inf) M1 = 0
        
        Chat = 1 - (nT - nt)/nT * Q1/sum(x) * (1 - ms / (nT - nt))^M1
        if (rho == 1) Chat = 1
        if (k >= nT) Chat = 1
        
      }
      
      Chat
    })
  }
  
  sapply(m, function(i) ifelse(datatype != 'abundance', Sub2(i), Sub(i) ))
}


#' Printing iNEXTscaling object
#' 
#' \code{print.iNEXTscaling}: Print method for objects inheriting from class "iNEXTscaling"
#' @param x an \code{iNEXTscaling} object computed by \code{iNEXTscaling}.
#' @param ... additional arguments.
#' @return a list of three objects (see \code{iNEXTscaling} for more details) with simplified outputs and notes.
#' @export
print.iNEXTscaling <- function(x, ...) {
  site.n <- nrow(x[[1]])
  order.n <- paste(unique(x[[2]]$size_based$Order.q), collapse = ", ")
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep = "")
  cat("$class: iNEXTscaling\n\n")
  cat(names(x)[1], ": basic data information\n", sep = "")
  print(x[[1]])
  cat("\n")
  cat(names(x)[2],": diversity estimates with rarefied and extrapolated samples.\n", sep = "")
  cat("$size_based (LCL and UCL are obtained for fixed size.)\n")
  cat("\n")
  
  res <- lapply(x[[2]], function(y) {
    Assemblages <- unique(x[[2]]$size_based$Assemblage)
    tmp <- lapply(1:length(Assemblages),function(i) {
      y_each <- y[y$Assemblage==Assemblages[i],]
      m <- quantile(unlist(y_each[,3]), type = 1)
      y_each[unlist(y_each[,3]) %in% m,]
    })
    do.call(rbind,tmp)
  })
  
  print(data.frame(res[[1]]))
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage in each order q; call iNEXTscaling.object$", names(x)[2],
      "$size_based to view complete output.\n", sep = "")
  cat("\n")
  cat("$coverage_based (LCL and UCL are obtained for fixed coverage; interval length is wider due to varying size in bootstraps.)\n")
  cat("\n")
  print(data.frame(res[[2]]))
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage in each order q; call iNEXTscaling.object$", names(x[2]), 
      "$coverage_based to view complete output.\n", sep = "")
  cat("\n")
  cat(names(x)[3], ": asymptotic diversity estimates along with related statistics.\n", sep = "")
  print(x[[3]])
  return(invisible())
}

check.datatype <- function(data, datatype, nT = nT, empirical = FALSE) {
  
  DATATYPE <- c("abundance", "incidence_raw")
  
  if (is.na(pmatch(datatype, DATATYPE))) stop("invalid datatype")
  if (pmatch(datatype, DATATYPE) == -1) stop("ambiguous datatype")
  
  datatype <- match.arg(datatype, DATATYPE)
  
  if (sum(nT <= 0) != 0) stop("Number of sampling units should be a positive value.", call. = FALSE)
  
  if (datatype == "incidence_raw") {
    
    if (inherits(data, c("data.frame", "matrix")) & is.null(nT)) data = list('Assemblage_1' = as.matrix(data))
    
    if (inherits(data, c("data.frame", "matrix")) & !is.null(nT)) {
      
      if (ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
      
      datalist = list()
      
      ntmp <- 0
      
      for (i in 1:length(nT)) {
        
        datalist[[i]] <- data[, (ntmp + 1):(ntmp + nT[i])]
        
        ntmp <- ntmp + nT[i]
      }
      
      names(datalist) = names(nT)
      
      data = datalist
      
      if (is.null(names(data))) names(data) = paste0("Assemblage_", 1:length(data))
    }
    
    if (inherits(data, c("numeric", "integer", "double"))) data = list('Assemblage_1' = as.matrix(data))
    
    if (inherits(data, "list")) {
      
      data = lapply(data, function(i) as.matrix(i))
      
      if ( sum(sapply(data, function(y) sum(rowSums(y) > 0)) < 5) > 0 & !empirical) stop("To ensure reliable results, iNEXT.scaling requires sufficient data; the number of observed species should be at least five. 
", call. = FALSE)
      
      if ( sum(sapply(data, sum) == 0) > 0 & !empirical) stop("Data values are all zero in some assemblages. Please remove these assemblages.", call. = FALSE)
      
    }
    
    if (sum(sapply(data, function(x) sum(x > 1))) != 0) stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
    
    if (sum(sapply(data, ncol) <= 3) != 0 & !empirical) stop("The number of sampling units in some assemblages is too small. Please provide additional sampling unit data.", call. = FALSE)
    
  }
  
  if (datatype == "abundance") {
    
    if (inherits(data, "list")) {
      
      if (length(data) == 1) {
        
        dat = as.matrix(data[[1]])
        
        if (is.null(names(data))) colnames(dat) = "Assemblage_1" else colnames(dat) = names(data)
        
        data = dat
        
      } else {
        
        region_names = if (is.null(names(data))) paste0("Assemblage_", 1:length(data)) else names(data)
        
      }
      
    } else if (inherits(data, c("numeric", "integer", "double"))) {
      
      data = list("Assemblage_1" = data)
      
    } else if (inherits(data, c("data.frame", "matrix"))){
      
      datalist <- lapply(1:ncol(data), function(i) data[,i])
      
      if (is.null(colnames(data))) names(datalist) <-  paste0("Assemblage", 1:ncol(data)) else names(datalist) <- colnames(data)
      
      data <- datalist
    }
    
    if ( ((datatype == "abundance" & sum(sapply(data, sum) < 5) > 0 ) |
          (datatype == "incidence_raw" & sum(sapply(data, sum) < 5) > 0 )) & !empirical  ) stop("To ensure reliable results, iNEXT.scaling requires sufficient data; the number of observed species should be at least five.
", call. = FALSE)
    
    if ( (datatype == "abundance" & sum(sapply(data, sum) == 0) > 0) |
         (datatype == "incidence_raw" & sum(sapply(data, sum) == 0) > 0)) stop("Data values are all zero in some assemblages. Please remove these assemblages.", call. = FALSE)
    
  }
  
  return(list(datatype, data))
  }

check.rho <- function(data, rho) {
  
  if (is.null(rho)) stop('Please set the rho in the function.', call. = FALSE)
  
  if (sum((rho < 0) | (rho > 1) | (is.numeric(rho) == F)) > 0) stop('Please enter value between zero and one for sampling fraction (rho).', call. = FALSE)
  
  if (length(rho) == 1) rho = rep(rho, length(data))
  
  if (length(rho) != length(data)) stop('The length of rho vector does not equal to the number of assemblages. Please check the setting of rho.', call. = FALSE)
  
  return(rho)
}

check.q <- function(q) {
  
  if (!inherits(q, "numeric")) stop("invalid class of order q, q should be a postive value/vector of numeric object", call. = FALSE)
  
  if (min(q) < 0) {
    warning("ambigous of order q, we only compute postive q", call. = FALSE)
    q <- q[q >= 0]
  }
  
  return(q)
}

check.conf <- function(conf) {
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F)) stop('Please enter value between zero and one for confident interval.', call. = FALSE)
  
  return(conf)
  }

check.nboot <- function(nboot) {
  
  if ((nboot < 0) | (is.numeric(nboot) == F)) stop('Please enter non-negative integer for nboot.', call. = FALSE)
  
  return(nboot)
  }

check.base <- function(base) {
  
  BASE <- c("size", "coverage")
  
  if (is.na(pmatch(base, BASE))) stop("invalid datatype")
  
  if (pmatch(base, BASE) == -1) stop("ambiguous datatype")
  
  base <- match.arg(base, BASE)
  
  return(base) 
  }

check.size <- function(data, rho, datatype, size, endpoint, knots) {
  
  if (length(knots) != length(data)) knots <- rep(knots, length(data))
  
  if (is.null(size)) {
    
    if (is.null(endpoint)) {
      
      if (datatype == "abundance") {
        
        endpoint <- sapply(1:length(data), function(i) {
          
          index = ifelse(rho[i] < 0.2, 2, 3)
          min(sum(data[[i]]) / rho[i], index * sum(data[[i]]))
        })
        
      } else if (datatype == "incidence_freq") {
        
        endpoint <- sapply(1:length(data), function(i) {
          
          index = ifelse(rho[i] < 0.2, 2, 3)
          min(data[[i]][1] / rho[i], index * data[[i]][1])
        })
        
      } else if (datatype == "incidence_raw"){
        
        endpoint <- sapply(1:length(data), function(i) {
          
          index = ifelse(rho[i] < 0.2, 2, 3)
          min(ncol(data[[i]]) / rho[i], index * ncol(data[[i]]))
        })
        
      }
      
    } else {
      
      if(length(endpoint) != length(data)){
        endpoint <- rep(endpoint, length(data))
      }
      
    }
    
    size <- lapply(1:length(data), function(i){
      
      if (datatype == "abundance") {
        
        ni <- sum(data[[i]])
        
      } else if (datatype == "incidence_freq") {
        
        ni <- data[[i]][1]
        
      } else if (datatype == "incidence_raw") {
        
        ni <- ncol(data[[i]])
      }
      
      if (endpoint[i] <= ni){
        
        mi <- floor(seq(1,endpoint[i], length.out = knots[i]))
        
      } else {
        
        mi <- floor(c(seq(1, ni, length.out = floor(knots[i]/2)), seq(ni + 1, endpoint[i], length.out = knots[i] - floor(knots[i]/2))))
      }
      
      if(sum(mi < 0) > 0) stop("Sample size (or number of sampling units) cannot be a negative value.", call. = FALSE)
      
      unique(mi)
    })
    
  } else {
    
    if (inherits(size, c("numeric", "integer", "double"))) size <- list(size = size)
    
    if (length(size) != length(data)) size <- lapply(1:length(data), function(x) size[[1]])
    
    size <- lapply(1:length(data),function(i){
      
      if(datatype == "abundance") {
        
        ni <- sum(data[[i]])
        
      } else if (datatype == "incidence_freq") {
        
        ni <- data[[i]][1]
        
      } else if (datatype == "incidence_raw") {
        
        ni <- ncol(data[[i]])
      }
      
      if ( (sum(size[[i]] == ni) == 0) & (sum(size[[i]] > ni) != 0) & (sum(size[[i]] < ni) != 0) ) mi <- sort(c(ni,size[[i]])) else mi <- sort(size[[i]])
      
      if(sum(mi < 0) > 0) stop("Sample size (or number of sampling units) cannot be a negative value.", call. = FALSE)
      
      unique(mi)
    })
  }
  
  
  return(size)
}

check.level <- function(data, rho, datatype, base, level) {
  
  index = sapply(rho, function(i) ifelse(i < 0.2, 2, 3))
  
  if (is.null(level) & base == 'size') {
    
    if (datatype == "abundance") {
      
      level <- sapply(1:length(data), function(i) min(index[i] * sum(data[[i]]), sum(data[[i]]) / rho[i]))
      
    } else if(datatype == "incidence_freq") {
      
      level <- sapply(1:length(data), function(i) min(index[i] * data[[i]][1], data[[i]][1] / rho[i]))
      
    } else if(datatype == "incidence_raw") {
      
      level <- sapply(1:length(data), function(i) min(index[i] * ncol(data[[i]]), ncol(data[[i]]) / rho[i]))
      
    }
    
    level <- min(level)
    
  } else if (is.null(level) & base == 'coverage') {
    
    if (datatype=='abundance') {
      
      level <- sapply(1:length(data), function(i) Coverage(data[[i]], rho = rho[i], datatype, m = index[i] * sum(data[[i]])))
      
    } else if (datatype=='incidence_freq') {
      
      level <- sapply(1:length(data), function(i) Coverage(data[[i]], rho = rho[i], datatype, m = index[i] * data[[i]][1]))
      
    } else if (datatype=='incidence_raw') {
      
      level <- sapply(1:length(data), function(i) Coverage(data[[i]], rho = rho[i], datatype, m = index[i] * ncol(data[[i]])))
      
    }
    
    level <- min(level)
  }
  
  if(base == "size" & sum(level < 0) > 0) stop("Sample size (or number of sampling units) cannot be a negative value.", call. = FALSE)
  
  if(base == "coverage" & sum(level < 0 | level > 1) > 0) stop("The sample coverage values should be between zero and one.", call. = FALSE)  
  
  return(level)
}


