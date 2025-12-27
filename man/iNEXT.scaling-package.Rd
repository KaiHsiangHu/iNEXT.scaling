\name{iNEXT.scaling-package}
\alias{iNEXT.scaling-package}
\alias{iNEXT.scaling}
\docType{package}
\title{
Interpolation and Extrapolation for upscaling biodiversity
}
\description{
\code{iNEXT.scaling} (iNterpolation and EXTrapolation for upscaling biodiversity) is a sequel to \code{iNEXT} (Hsieh et al., 2016).
\code{iNEXT.scaling} estimates Hill numbers (also called taxonomic diversity) under sampling without replacement model. An online version "iNEXT.scaling Online"
(https://kaihsiang-hu.shinyapps.io/iNEXT_scaling/) is also available for users without an R background.

\code{iNEXT.scaling} features two statistical analyses (non-asymptotic and asymptotic):

(1) A non-asymptotic approach based on interpolation and extrapolation for Hill numbers

\code{iNEXT.scaling} computes the estimated Hill numbers for standardized samples with a common sample size or sample completeness. This approach aims to compare diversity estimates for equally-large (with a common sample size) or equally-complete (with a common sample coverage) samples; it is based on the seamless rarefaction and extrapolation (R/E) sampling curves of Hill numbers for q = 0, 1 and 2. \code{iNEXT.scaling} offers three types of R/E sampling curves:

* Sample-size-based (or size-based) R/E sampling curves: This type of sampling curve plots the diversity estimates with respect to sample size.

* Coverage-based R/E sampling curves: This type of sampling curve plots the diversity estimates with respect to sample coverage.

* Sample completeness curve: This curve depicts how sample coverage varies with sample size. The sample completeness curve provides a bridge between the size- and coverage-based R/E sampling curves.

(2) An asymptotic approach to infer asymptotic diversity for Hill numbers

\code{iNEXT.scaling} computes the estimated asymptotic diversity and also plots diversity profiles (q-profiles) for q between 0 and 2, in comparison with the observed diversity. Typically, the asymptotic estimates for q \eqn{\geq} 1 are reliable, but for q < 1 (especially for q = 0, species richness), the asymptotic estimates represent only lower bounds. \cr

There are six main functions in this package: \cr

1. \code{iNEXTscaling} computes standardized diversity of order q = 0, 1 and 2 for rarefied and extrapolated samples at specified sample coverage values and sample sizes. \cr

2. \code{ggiNEXTscaling} visualizes the output from the function \code{iNEXTscaling}. \cr

3. \code{ObsAsyscaling} computes observed and asymptotic diversity of order q between 0 and 2 (in increments of 0.2) for Hill numbers. \cr

4. \code{ggObsAsyscaling} visualizes the output from the function \code{ObsAsyscaling}. \cr

5. \code{estimatescaling} computes diversity of order q = 0, 1 and 2 with a particular set of user-specified level of sample sizes or sample coverage values. \cr

6. \code{DataInfoscaling} provides basic data information based on the observed data. \cr

}
\author{
Kai-Hsiang Hu and Anne Chao \cr

Maintainer: Kai-Hsiang Hu <river19951224@gmail.com>
}
\keyword{package}
