#' Data information for reference samples
#' 
#' \code{DataInfoscaling} provides basic data information for diversity based on a reference sample.
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified. 
#' @param rho the sampling fraction can be input as a vector for each assemblage, or specify a single numeric value to apply to all assemblages.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data in a single matrix/data.frame) a vector of positive integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc.
#' 
#' @return a data.frame including basic data information.\cr\cr 
#' For abundance data, basic information includes assemblage name (\code{Assemblage}), sample size in the reference sample (\code{n}), 
#' total abundance in the overall assemblage (\code{N}), sampling fraction of the reference sample (\code{rho}), 
#' observed species richness in the reference sample (\code{S.obs}), 
#' sample coverage estimates of the reference sample (\code{SC(n)}), sample coverage estimate for twice the reference sample size (\code{SC(2n)}),
#' the first five species abundance counts (\code{f1}--\code{f5}).\cr
#'  
#' For incidence data, the basic information includes assemblage name (\code{Assemblage}), number of sampling units (\code{T}), 
#' total number of sampling units in the overall assemblage (\code{T*}), sampling fraction of the reference sample (\code{rho}), 
#' total number of incidences (\code{U}), observed species richness (\code{S.obs}), 
#' sample coverage estimates of the reference sample (\code{SC(T)}), sample coverage estimate for twice the reference sample size
#' (\code{SC(2T)}), as well as the first five species incidence frequency counts (\code{Q1}--\code{Q5}).\cr
#'  
#' 
#' @examples
#' # Taxonomic diversity for abundance data
#' data(Brazil_rainforest_abun_data)
#' DataInfoscaling(Brazil_rainforest_abun_data, rho = 0.3, datatype = "abundance")
#' 
#' 
#' # Taxonomic diversity for incidence data
#' data(Fish_incidence_data)
#' DataInfoscaling(Fish_incidence_data, rho = 0.3, datatype = "incidence_raw")
#' 
#'
#' @export
DataInfoscaling <- function(data, rho, datatype = "abundance", nT = NULL){
  
  checkdatatype = check.datatype(data, datatype, nT = nT, empirical = TRUE)
  datatype = checkdatatype[[1]]
  data = checkdatatype[[2]]
  rho = check.rho(data, rho)
  
  out <- TDinfo(data, rho, datatype)
  
  return(out)
}


#' @useDynLib iNEXT.scaling, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



#' iNterpolation and EXTrapolation with three dimensions of biodiversity
#' 
#' \code{iNEXTscaling} mainly computes standardized estimates with a common sample size or sample coverage for orders q = 0, 1 and 2. It also computes relevant information/statistics.\cr\cr 
#' Relevant data information is summarized in the output \code{$TDInfo}. Diversity estimates for rarefied and extrapolated samples are provided in the output \code{$TDiNextEst}, which includes two data frames (\code{"$size_based"} and \code{"$coverage_based"}) based on two different standardizations; in the size-based standardization, all samples are standardized to a common target sample size, whereas the in the latter standardization, all samples are standardized to a common target level of sample coverage. The asymptotic diversity estimates for q = 0, 1 and 2 are provided in the list \code{$TDAsyEst}.\cr\cr 
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified.
#' @param rho the sampling fraction can be input as a vector for each assemblage, or specify a single numeric value to apply to all assemblages.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
#' If \code{NULL}, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}.
#' @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
#' If \code{NULL}, then \code{endpoint} \code{=} double reference sample size.
#' @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
#' each knot represents a particular sample size for which diversity estimate will be calculated.  
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXTscaling()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNEXTscaling()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data in a single matrix/data.frame) a vector of positive integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' 
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import tibble
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats optimize
#' @importFrom stats dhyper
#' @importFrom stats quantile
#' @importFrom grDevices hcl
#' 
#' @return a list of three objects: \cr\cr
#' (1) \code{$TDInfo} for summarizing data information for q = 0, 1 and 2. Refer to the output of \code{DataInfoscaling} for details. \cr\cr
#' (2) \code{$TDiNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics. There are two data frames: \code{"$size_based"} and \code{"$coverage_based"}. \cr\cr
#'    In \code{"$size_based"}, the output includes:
#'    \item{Assemblage}{the name of assemblage.} 
#'    \item{Order.q}{the diversity order of q.}
#'    \item{m, mT}{the target sample size (or number of sampling units for incidence data).}
#'    \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target sample size is less than, equal to, or greater than the size of the reference sample.}
#'    \item{qTD}{the estimated diversity estimate.}
#'    \item{qTD.LCL and qTD.UCL}{the bootstrap lower and upper confidence limits for the diversity of order q at the specified level (with a default value of 0.95).}
#'    \item{SC}{the standardized coverage value.}
#'    \item{SC.LCL, SC.UCL}{the bootstrap lower and upper confidence limits for coverage at the specified level (with a default value of 0.95).}
#'  Similar output is obtained for \code{"$coverage_based"}. \cr\cr
#' (3) \code{$TDAsyEst} for showing asymptotic diversity estimates along with related statistics: 
#'    \item{Assemblage}{the name of assemblage.} 
#'    \item{qTD}{the diversity order of q.}
#'    \item{TD_obs}{the observed diversity.}
#'    \item{TD_asy}{the asymptotic diversity estimate.}
#'    \item{s.e.}{standard error of asymptotic diversity.}
#'    \item{qTD.LCL and qTD.UCL}{the bootstrap lower and upper confidence limits for asymptotic diversity at the specified level (with a default value of 0.95).}
#' 
#' 
#' @examples
#' \donttest{
#' # Compute standardized estimates of taxonomic diversity for abundance data with order q = 0, 1, 2
#' data(Brazil_rainforest_abun_data)
#' output_TD_abun <- iNEXTscaling(Brazil_rainforest_abun_data, rho = 0.3, 
#'                                q = c(0, 1, 2), datatype = "abundance")
#' output_TD_abun
#' 
#' 
#' # Compute standardized estimates of taxonomic diversity for incidence data with order q = 0, 1, 2
#' data(Fish_incidence_data)
#' output_TD_inci <- iNEXTscaling(Fish_incidence_data, rho = 0.3, 
#'                                q = c(0, 1, 2), datatype = "incidence_raw")
#' output_TD_inci
#' 
#' }
#' 
#' 
#' @export
iNEXTscaling <- function(data, rho, q = c(0, 1, 2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, nboot = 50, conf = 0.95, nT = NULL) {
  
  data.original = data
  datatype.original = datatype
  
  checkdatatype = check.datatype(data, datatype, nT = nT) 
  datatype = checkdatatype[[1]]
  data = checkdatatype[[2]]
  
  rho = check.rho(data, rho)
  q = check.q(q)
  conf = check.conf(conf)
  nboot = check.nboot(nboot)
  size = check.size(data, rho, datatype, size, endpoint, knots)
  
  Fun <- function(x, rho, q, size, assem_name){
    
    if(datatype == "abundance"){
      
      out <- iNEXT.Ind(Spec=x, rho, q=q, m=size, knots=knots, nboot=nboot, conf=conf)
    }
    
    if(datatype == "incidence_raw"){
      
      if(sum(x)==0) stop("Zero incidence frequencies in one or more sample sites")
      
      out <- iNEXT.Sam(Spec=x, rho, q=q, t=size, knots=knots, nboot=nboot, conf=conf)
    }
    
    out <- lapply(out, function(out_) cbind(Assemblage = assem_name, out_))
    
    out
  }
  
  z <- qnorm(1-(1-conf)/2)
  
  if(is.null(names(data))){
    names(data) <- sapply(1:length(data), function(i) paste0('assemblage',i))
  }
  out <- lapply(1:length(data), function(i) {
    tmp <- Fun(data[[i]], rho[i], q, size[[i]], names(data)[i])
    tmp
  })
  
  out <- list(size_based = do.call(rbind, lapply(out, function(out_){out_[[1]]})),
              coverage_based = do.call(rbind, lapply(out, function(out_){out_[[2]]})))
  
  index <- rbind(asyTD(data, rho, datatype, c(0, 1, 2), nboot, conf),
                 obsTD(data, rho, datatype, c(0, 1, 2), nboot, conf))
  index = index[order(index$Assemblage),]
  LCL <- index$qTD.LCL[index$Method=='Asymptotic']
  UCL <- index$qTD.UCL[index$Method=='Asymptotic']
  index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qTD')
  index <- cbind(index,se = (UCL - index$Asymptotic)/z,LCL,UCL)
  index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
  index[,3:4] = index[,4:3]
  colnames(index) <- c("Assemblage", "qTD", "TD_obs", "TD_asy", "s.e.", "qTD.LCL", "qTD.UCL")
  
  
  out$size_based$Assemblage <- as.character(out$size_based$Assemblage)
  out$coverage_based$Assemblage <- as.character(out$coverage_based$Assemblage)
  
  info <- DataInfoscaling(data.original, rho, datatype.original, nT)
  
  out <- list("TDInfo"=info, "TDiNextEst"=out, "TDAsyEst"=index)
  
  
  if (datatype != 'abundance'){
    out[[2]]$size_based <- rename(out[[2]]$size_based, c("mT" = "nT"))
    out[[2]]$coverage_based <- rename(out[[2]]$coverage_based, c("mT" = "nT"))
  }
  
  class(out) <- c("iNEXTscaling")
  
  return(out)
}


#' ggplot2 extension for an iNEXTscaling object
#' 
#' \code{ggiNEXTscaling} is a \code{ggplot} extension for an \code{iNEXTscaling} object to plot sample-size- and coverage-based rarefaction/extrapolation sampling curves along with a bridging sample completeness curve.
#' @param output an \code{iNEXTscaling} object computed by \code{iNEXTscaling}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).            
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation (\code{facet.var = "None"}); 
#'  a separate plot for each diversity order (\code{facet.var = "Order.q"}); 
#'  a separate plot for each assemblage (\code{facet.var = "Assemblage"}); 
#'  a separate plot for each combination of diversity order and assemblage (\code{facet.var = "Both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var = "None"}); 
#'  use different colors for diversity orders (\code{color.var = "Order.q"}); 
#'  use different colors for assemblages/sites (\code{color.var = "Assemblage"}); 
#'  use different colors for combinations of diversity order and assemblage (\code{color.var = "Both"}).  
#' @return a \code{ggplot2} object for sample-size-based rarefaction/extrapolation curve (\code{type = 1}), sample completeness curve (\code{type = 2}), and coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' 
#' 
#' @examples
#' \donttest{
#' # Plot three types of curves of taxonomic diversity with facet.var = "Assemblage"
#' # for abundance data with order q = 0, 1, 2
#' data(Brazil_rainforest_abun_data)
#' output_TD_abun <- iNEXTscaling(Brazil_rainforest_abun_data, rho = 0.3, 
#'                                q = c(0, 1, 2), datatype = "abundance")
#' ggiNEXTscaling(output_TD_abun, facet.var = "Assemblage")
#' 
#' 
#' # Plot three types of curves of taxonomic diversity for incidence data with order q = 0, 1, 2
#' data(Fish_incidence_data)
#' output_TD_inci <- iNEXTscaling(Fish_incidence_data, rho = 0.3, 
#'                                q = c(0, 1, 2), datatype = "incidence_raw")
#' ggiNEXTscaling(output_TD_inci)
#' }
#' 
#' 
#' @export
ggiNEXTscaling = function(output, type = 1:3, facet.var = "Assemblage", color.var = "Order.q"){
  
  class = 'TD'
  plottable = output$TDiNextEst
  plottable$size_based = rename(plottable$size_based, c('qD' = 'qTD', 'qD.LCL' = 'qTD.LCL', 'qD.UCL' = 'qTD.UCL'))
  plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qTD', 'qD.LCL' = 'qTD.LCL', 'qD.UCL' = 'qTD.UCL'))
  
  SPLIT <- c("None", "Order.q", "Assemblage", "Both")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  TYPE <-  c(1, 2, 3)
  if(sum(!(type %in% TYPE)) >= 1)
    stop("invalid plot type")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  
  if(facet.var == "Order.q") color.var <- "Assemblage"
  if(facet.var == "Assemblage") color.var <- "Order.q"
  
  if ('m' %in% colnames(plottable$size_based) & 'm' %in% colnames(plottable$coverage_based)) datatype = 'abundance'
  if ('mT' %in% colnames(plottable$size_based) & 'mT' %in% colnames(plottable$coverage_based)) datatype = 'incidence'
  
  
  out = lapply(type, function(i) type_plot(x_list = plottable, i, class, datatype, facet.var, color.var))
  if (length(type) == 1) out = out[[1]]
  
  return(out)
}


type_plot = function(x_list, type, class, datatype, facet.var, color.var) {
  
  x_name <- colnames(x_list$size_based)[3]
  xlab_name <- ifelse(datatype == "incidence", "sampling units", "individuals")
  ylab_name = "Taxonomic diversity"
  
  if (type == 1) {
    
    output <- x_list$size_based
    output$y.lwr <- output$qD.LCL
    output$y.upr <- output$qD.UCL
    id <- match(c(x_name, "Method", "qD", "qD.LCL", "qD.UCL", "Assemblage", "Order.q"), names(output), nomatch = 0)
    output[,1:7] <- output[, id]
    
    xlab_name <- paste0("Number of ", xlab_name)
    
  } else if (type == 2) {
    
    output <- x_list$size_based
    if (length(unique(output$Order.q)) > 1) output <- subset(output, Order.q == unique(output$Order.q)[1])
    output$y.lwr <- output$SC.LCL
    output$y.upr <- output$SC.UCL
    id <- match(c(x_name, "Method", "SC", "SC.LCL", "SC.UCL", "Assemblage", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(output), nomatch = 0)
    output[,1:10] <- output[, id]
    
    xlab_name <- paste0("Number of ", xlab_name)
    ylab_name <- "Sample coverage"
    
  } else if (type == 3) {
    
    output <- x_list$coverage_based %>% tibble
    output$y.lwr <- output$qD.LCL
    output$y.upr <- output$qD.UCL
    id <- match(c("SC", "Method", "qD", "qD.LCL", "qD.UCL", "Assemblage", "Order.q", x_name), names(output), nomatch = 0)
    output[,1:8] <- output[, id]
    
    xlab_name <- "Sample coverage"
    
  }
  
  if (facet.var == "None" & color.var == "None" & length(unique(output$Order.q)) > 1 & length(unique(output$Assemblage)) > 1) {
    
    color.var <- "Order.q"
    facet.var <- "Assemblage"
    warning ("invalid color.var and facet.var setting, the iNEXTscaling object consists multiple orders and assemblage, change setting as Order.q and Assemblage")
    
  } else if (facet.var == "None" & color.var == "None" & length(unique(output$Order.q)) > 1) {
    
    color.var <- "Order.q"
    warning ("invalid color.var setting, the iNEXTscaling object consists multiple orders, change setting as Order.q")
    
  } else if (facet.var == "None" & color.var == "None" & length(unique(output$Assemblage)) > 1) { 
    
    color.var <- "Assemblage" 
    warning ("invalid color.var setting, the iNEXTscaling object consists multiple assemblage, change setting as Assemblage")
  }
  
  
  title <- c("Sample-size-based sampling curve", "Sample completeness curve", "Coverage-based sampling curve")[type]
  colnames(output)[1:7] <- c("x", "Method", "y", "LCL", "UCL", "Assemblage", "Order.q")
  
  if (color.var == "None") {
    
    if (levels(factor(output$Order.q)) > 1 & length(unique(output$Assemblage)) > 1) {
      warning ("invalid color.var setting, the iNEXTscaling object consists multiple assemblages and orders, change setting as Both")
      color.var <- "Both"
      output$col <- output$shape <- paste(output$Assemblage, output$Order.q, sep="-")
      
    } else if (length(unique(output$Assemblage)) > 1) {
      warning ("invalid color.var setting, the iNEXTscaling object consists multiple assemblages, change setting as Assemblage")
      color.var <- "Assemblage"
      output$col <- output$shape <- output$Assemblage
    } else if (levels(factor(output$Order.q)) > 1){
      warning ("invalid color.var setting, the iNEXTscaling object consists multiple orders, change setting as Order.q")
      color.var <- "Order.q"
      output$col <- output$shape <- factor(output$Order.q)
    } else {
      output$col <- output$shape <- rep(1, nrow(output))
    }
  } else if (color.var == "Order.q") {    
    
    output$col <- output$shape <- factor(output$Order.q)
  } else if (color.var == "Assemblage") {
    
    if (length(unique(output$Assemblage)) == 1) {
      warning ("invalid color.var setting, the iNEXTscaling object do not consist multiple assemblages, change setting as Order.q")
      output$col <- output$shape <- factor(output$Order.q)
    }
    output$col <- output$shape <- output$Assemblage
    
  } else if (color.var == "Both") {
    
    if (length(unique(output$Assemblage)) == 1) {
      warning ("invalid color.var setting, the iNEXTscaling object do not consist multiple assemblages, change setting as Order.q")
      output$col <- output$shape <- factor(output$Order.q)
    }
    output$col <- output$shape <- paste(output$Assemblage, output$Order.q, sep="-")
  }
  
  if (type == 2) output$col = output$shape = output$Assemblage
  
  data.sub = output
  tmp = output %>% filter(Method == "Observed") %>% mutate(Method = "Extrapolation")
  output$Method[output$Method == "Observed"] = "Rarefaction"
  output = rbind(output, tmp)
  output$lty <- factor(output$Method, levels = c("Rarefaction", "Extrapolation"))
  output$col <- factor(output$col)
  data.sub <- data.sub[which(data.sub$Method == "Observed"),]
  
  if (length(unique(output$Assemblage)) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(output$Assemblage))-8))
  }
  
  g <- ggplot(output, aes_string(x = "x", y = "y", colour = "col")) + 
    geom_line(aes_string(linetype = "lty"), lwd=1.5) +
    geom_point(aes_string(shape = "shape"), size=5, data = data.sub) +
    geom_ribbon(aes_string(ymin = "y.lwr", ymax = "y.upr", fill = "factor(col)", colour = "NULL"), alpha = 0.2) +
    scale_fill_manual(values = cbPalette) +
    scale_colour_manual(values = cbPalette) +
    guides(linetype = guide_legend(title = "Method"),
           colour = guide_legend(title = "Guides"), 
           fill = guide_legend(title = "Guides"), 
           shape = guide_legend(title = "Guides"))
  
  g = g + theme_bw() + 
    labs(x = xlab_name, y = ylab_name) + 
    ggtitle(title) + 
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          text = element_text(size = 16),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    guides(linetype = guide_legend(keywidth = 2.5))
  
  
  if (facet.var == "Order.q") {
    
    if(length(levels(factor(output$Order.q))) == 1 & type != 2){
      warning("invalid facet.var setting, the iNEXTscaling object do not consist multiple orders.")      
    } else {
      odr_grp <- labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      
      g <- g + facet_wrap( ~ Order.q, nrow = 1, labeller = odr_grp)
      
      if (color.var == "Both") {
        g <- g + guides(colour = guide_legend(title = "Guides", ncol = length(levels(factor(output$Order.q))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
      if(type == 2){
        g <- g + theme(strip.background = element_blank(), strip.text.x = element_blank())
        
      }
    }
  }
  
  if(facet.var == "Assemblage"){
    
    if(length(unique(output$Assemblage)) == 1) {
      warning("invalid facet.var setting, the iNEXTscaling object do not consist multiple assemblages")
    }else{
      g <- g + facet_wrap(. ~ Assemblage, nrow = 1)
      
      if(color.var == "Both"){
        g <- g + guides(colour = guide_legend(title = "Guides", nrow = length(levels(factor(output$Order.q)))),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  
  if(facet.var == "Both"){
    
    if(length(levels(factor(output$Order.q))) == 1 | length(unique(output$Assemblage)) == 1){
      warning("invalid facet.var setting, the iNEXTscaling object do not consist multiple assemblages or orders.")
    }else{
      odr_grp <- labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      
      g <- g + facet_wrap(Assemblage ~ Order.q, labeller = odr_grp)
      
      if(color.var == "both"){
        g <- g +  guides(colour = guide_legend(title = "Guides", nrow = length(levels(factor(output$Assemblage))), byrow = TRUE),
                         fill = guide_legend(title = "Guides"))
      }
    }
  }
  
  return(g)
}


#' Compute diversity estimates with a particular set of sample sizes/coverages
#' 
#' \code{estimatescaling} computes diversity (Hill-Chao number with q = 0, 1 and 2) with a particular set of user-specified levels of sample sizes or sample coverages. If no sample sizes or coverages are specified, this function by default computes diversity estimates for the minimum sample coverage or minimum sample size among all samples extrapolated to double reference sizes.
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified. 
#' @param rho the sampling fraction can be input as a vector for each assemblage, or specify a single numeric value to apply to all assemblages.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param base selection of sample-size-based (\code{base = "size"}) or coverage-based (\code{base = "coverage"}) rarefaction and extrapolation.
#' @param level A numerical vector specifying the particular sample sizes or sample coverages (between 0 and 1) for which diversity estimates (q =0, 1 and 2) will be computed. \cr
#' If \code{base = "coverage"} (default) and \code{level = NULL}, then this function computes the diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes. \cr
#' If \code{base = "size"} and \code{level = NULL}, then this function computes the diversity estimates for the minimum sample size among all samples extrapolated to double reference sizes. 
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data in a single matrix/data.frame) a vector of positive integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' 
#' @return a data.frame of diversity table including the following arguments: (when \code{base = "coverage"})
#' \item{Assemblage}{the name of assemblage.}
#' \item{Order.q}{the diversity order of q.}
#' \item{SC}{the target standardized coverage value.}
#' \item{m, mT}{the corresponding sample size (or number of sampling units) for the standardized coverage value.}
#' \item{qTD}{the estimated diversity of order q for the target coverage value. The estimate for complete coverage (when \code{base = "coverage"} and \code{level = 1}, or \code{base = "size"} and \code{level = Inf}) represents the estimated asymptotic diversity.}
#' \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample.}
#' \item{s.e.}{standard error of diversity estimate.}
#' \item{qTD.LCL and qTD.UCL}{the bootstrap lower and upper confidence limits for the diversity of order q at the specified level (with a default value of 0.95).}
#' Similar output is obtained for \code{base = "size"}. \cr\cr
#' 
#' 
#' @examples
#' # Taxonomic diversity for abundance data with two target coverages (93% and 97%)
#' data(Brazil_rainforest_abun_data)
#' output_est_cov_abun <- estimatescaling(Brazil_rainforest_abun_data, rho = 0.3, q = c(0, 1, 2), 
#'                                        datatype = "abundance", base = "coverage", 
#'                                        level = c(0.93, 0.97))
#' output_est_cov_abun
#' 
#' 
#' data(Brazil_rainforest_abun_data)
#' output_est_size_abun <- estimatescaling(Brazil_rainforest_abun_data, rho = 0.3, q = c(0, 1, 2), 
#'                                         datatype = "abundance", base = "size", 
#'                                         level = c(1000, 2000))
#' output_est_size_abun 
#' 
#' 
#' # Taxonomic diversity for incidence data with two target coverages (97.5% and 99%)
#' data(Fish_incidence_data)
#' output_est_cov_inci <- estimatescaling(Fish_incidence_data, rho = 0.3, q = c(0, 1, 2), 
#'                                        datatype = "incidence_raw", base = "coverage", 
#'                                        level = c(0.975, 0.99))
#' output_est_cov_inci
#' 
#' 
#' data(Fish_incidence_data)
#' output_est_size_inci <- estimatescaling(Fish_incidence_data, rho = 0.3, q = c(0, 1, 2), 
#'                                         datatype = "incidence_raw", base = "coverage", 
#'                                         level = c(0.975, 0.99))
#' output_est_size_inci
#' 
#' 
#' @export
estimatescaling <- function(data, rho, q = c(0,1,2), datatype = "abundance", base = "coverage", level = NULL, nboot = 50, conf = 0.95, nT = NULL) {
  
  checkdatatype = check.datatype(data, datatype, nT = nT)
  datatype = checkdatatype[[1]]
  data = checkdatatype[[2]]
  
  rho = check.rho(data, rho)
  q = check.q(q)
  conf = check.conf(conf)
  nboot = check.nboot(nboot)
  base = check.base(base)
  level = check.level(data, rho, datatype, base, level)
  
  if (base == "size") {
    
    out <- invSize(data, rho, q, datatype, size = level, nboot, conf = conf)
    
  } else if (base == "coverage") {
    
    out <- invChat(data, rho, q, datatype, C = level, nboot, conf = conf)
  }
  out$qTD.LCL[out$qTD.LCL<0] <- 0
  
  if (datatype == "incidence_raw") out <- rename(out, c("mT" = "nT"))
  
  return(out)
}


#' Asymptotic diversity and observed diversity of order q
#' 
#' \code{ObsAsyscaling} computes observed and asymptotic diversity of order q between 0 and 2 (in increments of 0.2) for diversity; these values with different order q can be used to depict a q-profile in the \code{ggObsAsyscaling} function.\cr\cr 
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified. 
#' @param rho the sampling fraction can be input as a vector for each assemblage, or specify a single numeric value to apply to all assemblages.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, by = 0.2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data in a single matrix/data.frame) a vector of positive integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param method Select \code{'Asymptotic'} or \code{'Observed'}.
#' 
#' @return a data frame including the following information/statistics: 
#' \item{Assemblage}{the name of assemblage.}
#' \item{Order.q}{the diversity order of q.}
#' \item{qTD}{the estimated asymptotic diversity or observed diversity of order q.} 
#' \item{s.e.}{standard error of diversity.}
#' \item{qTD.LCL and qTD.UCL}{the bootstrap lower and upper confidence limits for the diversity of order q at the specified level (with a default value of 0.95).}
#' \item{Method}{\code{"Asymptotic"} means asymptotic diversity and \code{"Observed"} means observed diversity.}
#' 
#' 
#' @examples
#' # Compute the observed and asymptotic taxonomic diversity for abundance data
#' # with order q between 0 and 2 (in increments of 0.2 by default)
#' data(Brazil_rainforest_abun_data)
#' output_ObsAsy_TD_abun <- ObsAsyscaling(Brazil_rainforest_abun_data, rho = 0.3, 
#'                                        datatype = "abundance")
#' output_ObsAsy_TD_abun
#' 
#' 
#' # Compute the observed and asymptotic taxonomic diversity for incidence data
#' # with order q between 0 and 2 (in increments of 0.2 by default).
#' data(Fish_incidence_data)
#' output_ObsAsy_TD_inci <- ObsAsyscaling(Fish_incidence_data, rho = 0.3, 
#'                                        datatype = "incidence_raw")
#' output_ObsAsy_TD_inci
#' 
#' 
#' @export
ObsAsyscaling <- function(data, rho, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, nT = NULL, method = c('Asymptotic', 'Observed')) {
  
  if ("Asymptotic" %in% method) checkdatatype = check.datatype(data, datatype, nT = nT) else
    checkdatatype = check.datatype(data, datatype, nT = nT, empirical = TRUE)
  
  datatype = checkdatatype[[1]]
  data = checkdatatype[[2]]
  
  rho = check.rho(data, rho)
  q = check.q(q)
  conf = check.conf(conf)
  nboot = check.nboot(nboot)
  
  
  if (sum(method == "Asymptotic") == length(method)) 
    
    out = asyTD(data, rho, datatype, q, nboot, conf) else if (sum(method == "Observed") == length(method)) 
      
      out = obsTD(data, rho, datatype, q, nboot, conf) else if (sum(method == c("Asymptotic", "Observed")) == length(method)) 
        
        out = rbind(asyTD(data, rho, datatype, q, nboot, conf), 
                    obsTD(data, rho, datatype, q, nboot, conf))
  
  return(out)
}


#' ggplot2 extension for plotting q-profile
#'
#' \code{ggObsAsyscaling} is a \code{ggplot2} extension for an \code{ObsAsyscaling} object to plot q-profile (which depicts the observed diversity and asymptotic diversity estimate with respect to order q) for q between 0 and 2 (in increments of 0.2).\cr\cr 
#' In the plot of profiles, only confidence intervals of the asymptotic diversity will be shown when both the observed and asymptotic diversity estimates are computed.
#' 
#' @param output the output of the function \code{ObsAsyscaling}.\cr
#' @return a q-profile based on the observed diversity and the asymptotic diversity estimate.\cr\cr
#' 
#' @examples
#' # Plot q-profile of taxonomic diversity for abundance data
#' # with order q between 0 and 2 (in increments of 0.2 by default).
#' data(Brazil_rainforest_abun_data)
#' output_ObsAsy_TD_abun <- ObsAsyscaling(Brazil_rainforest_abun_data, rho = 0.3, 
#'                                        datatype = "abundance")
#' ggObsAsyscaling(output_ObsAsy_TD_abun)
#' 
#' 
#' # Plot q-profile of taxonomic diversity for incidence data
#' # with order q between 0 and 2 (in increments of 0.2 by default)
#' data(Fish_incidence_data)
#' output_ObsAsy_TD_inci <- ObsAsyscaling(Fish_incidence_data, rho = 0.3, 
#'                                        datatype = "incidence_raw")
#' ggObsAsyscaling(output_ObsAsy_TD_inci)
#' 
#' 
#' @export
ggObsAsyscaling <- function(output){
  
  if (sum(unique(output$Method) %in% c("Asymptotic", "Observed")) == 0)
    stop("Please use the output from specified function 'ObsAsyscaling'")
  
  out = ggplot(output, aes(x = Order.q, y = qTD, colour = Assemblage, fill = Assemblage))
  
  if (length(unique(output$Method)) == 1) {
    out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qTD.LCL, ymax = qTD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
    
    if (unique(output$Method == 'Asymptotic')) out = out + labs(x = 'Order q', y = 'Asymptotic taxonomic diversity')
    if (unique(output$Method == 'Observed')) out = out + labs(x = 'Order q', y = 'Observed taxonomic diversity')
  } else {
    out = out + geom_line(aes(lty = Method), size = 1.5) + 
      geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qTD.LCL, ymax = qTD.UCL), linetype = 0, alpha = 0.2)
    
    out = out + labs(x = 'Order q', y = 'Taxonomic diversity')
  }
  
  if (length(unique(output$Assemblage)) <= 8){
    
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    
  }else{
    
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(output$Assemblage))-8))
  }
  
  out = out +
    scale_colour_manual(values = cbPalette) + theme_bw() + 
    scale_fill_manual(values = cbPalette) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          text = element_text(size = 16),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    guides(linetype = guide_legend(keywidth = 2.5))
  
  return(out)
}


ggplotColors <- function(g){
  d <- 360/g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d,g - 1))) # Create cumulative sums to define hue values
  hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}


## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("qTD", "Assemblage", "qTD.LCL", "qTD.UCL",
                         "Method", "Order.q"))




