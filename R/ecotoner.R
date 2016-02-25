#' Tool to extract geographic and geometric measures of ecotones for R.
#'
#' @section Methods to locate a transect:
#'
#' \enumerate{
#'		\item Select random points in edge cells of vegetation 1
#'		\item Locate a second point within a neighborhood of specified size by one of four methods:
#'		\itemize{
#'			\item at the highest elevation (\code{\linkS4class{EcotonerSettings@transect_type == 1}})
#'			\item following the steepest slope outside a minimum distance (\code{\linkS4class{EcotonerSettings@transect_type == 2}})
#'			\item by the flowpath (\code{\linkS4class{EcotonerSettings@transect_type == 3}})
#'			\item by sufficient area of homogeneous aspect (\code{\linkS4class{EcotonerSettings@transect_type == 4}})
#'		}
#'		\item Extend line connecting these two points in both directions
#'		\item Trim line to the section between the lowest and highest elevation
#'	}
#'
#' @section Determine if transect is a suitable ecotone / ecological boundary:
#'
#' \enumerate{
#'		\item Transect is empty (i.e., trivial condition)
#'		\item Definition of 'zone of ecotone / ecological boundary' is not met
#'	}
#'
#' @section Definition of 'zone of ecotone / ecological boundary':
#'
#' Operational definition of zone of ecological boundary (v5): 
#' Process:
#' \enumerate{
#'		\item Apply definition to vegetation 1 and 2
#'		\item Apply definition to majority type of vegetation 1 and majority type of vegetation 2
#' }
#' Definition:
#'	\itemize{
#'		\item Transect orientation (from low to high elevation, i.e., from left to right)
#'		\item Density = rolling mean [over 90 m = 3 cell widths] of the mean column density (over the 200 cell wide band) at a 0.01 resolution
#'		\item Overall start :=
#'		\itemize{
#'			\item Leftmost position within minimum zone (= first run of sufficiently low values, <= \code{\linkS4class{EcotonerSettings@vegDefiningDensityTransect_low}}) of veg2 such that
#'			\item There is at least one point of the maximum zone (= last run of sufficiently high values, >= \code{\linkS4class{EcotonerSettings@vegDefiningDensityTransect_high}}) of veg1 to the right and that
#'			\item All of veg1 between overall start and the rightmost point of the maximum zone of veg1 have sufficient density (i.e., > \code{\linkS4class{EcotonerSettings@vegDefiningDensityTransectExtended_min}})
#'		}
#'		\item Overall end :=
#'		\itemize{
#'			\item Rightmost position within minimum zone (=last low run after the last high run) of veg1 such that
#'			\item Overall start is on the left,
#'			\item There is at least one point of the maximum zone (=first high run after the last low run) of veg2 to the left and that
#'			\item All of veg2 between the lefmost point of the maximum zone of veg1 and overall end have sufficient density (i.e., > \code{\linkS4class{EcotonerSettings@vegDefiningDensityTransectExtended_min}})
#'		}
#'		\item No remaining BSE or TF vegetation should attain a maximal density (i.e., < \code{\linkS4class{EcotonerSettings@vegDefiningDensityTransectExtended_min}})
#'		\item Zone of ecological boundary := from overall start to overall end
#'	}
#'
#' @docType package
#' @name ecotoner
NULL

##------ Package level variables
etr_vars <- new.env()


##------ Package uses S4 classes - they are defined in package:methods
#' @importFrom methods setClass setGeneric setMethod setValidity slot setReplaceMethod new initialize validObject
NULL
