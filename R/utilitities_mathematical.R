#' Determine whether a number is odd.
#'
#' @param x An integer (vector).
#' @return A logical: TRUE, if \code{x} is odd; FALSE, otherwise
#' @examples
#' is_odd((-3):3)
is_odd <- function(x) x %% 2 != 0


#' Calculate Root Mean Square (Error) RMS(E).
#'
#' @param x A numeric vector.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return If \code{x} are residuals, then RSME else RMS of \code{x}
#' @aliases RMS RMSE
#' @seealso \code{\link[stats]{sd}}
#' @examples
#' x <- runif(30)
#' y <- 0.1 + 0.3 * x + rnorm(30, 0, 0.5)
#' r <- resid(lm(y ~ x))
#' stats::sd(r)
#' root_mean_square(r)
root_mean_square <- function(x, na.rm = FALSE) {
	# Root Mean Square (Error, if x = residuals)
	sqrt(mean(x^2, na.rm = na.rm))
}

quantile_95th_of_abs <- function(x, na.rm = FALSE) {
	quantile(abs(x), probs = 0.95, na.rm = na.rm)
}

max_of_abs <- function(x, na.rm = FALSE) {
	max(abs(x), na.rm = na.rm)
}


#Circular functions: int=number of units in circle, e.g., for days: int=365; for months: int=12
circ_mean <- function(x, int, na.rm = FALSE) {
	if (!all(is.na(x))) {
		circ <- 2 * pi / int
		x_circ <- circular::circular(x * circ, type = "angles", units = "radians", rotation = "clock", modulo = "2pi")
		x_int <- circular::mean.circular(x_circ, na.rm = na.rm) / circ
		res <- round(as.numeric(x_int) - 1, 13) %% int + 1	# map 0 -> int; rounding to 13 digits: 13 was empirically derived for int={12, 365} and x=c((-1):2, seq(x-5, x+5, by=1), seq(2*x-5, 2*x+5, by=1)) assuming that this function will never need to calculate for x > t*int with t>2
	} else {
		res <- NA
	}
	
	res
}

circ_sum <- function(x, y, int) {
	circ <- 2 * pi / int
	as.numeric(circular::circular((x + y) * circ, type = "angles", units = "radians", rotation = "clock", modulo = "2pi") / circ)
}

circ_range <- function(x, int, na.rm = FALSE) {
	if (!all(is.na(x))) {
		circ <- 2 * pi / int
		x_circ <- circular::circular(x * circ, type = "angles", units = "radians", rotation = "clock", modulo = "2pi")
		x_int <- range(x_circ, na.rm = na.rm) / circ
		res <- as.numeric(x_int)
	} else {
		res <- NA
	}
	
	res
}

circ_sd <- function(x, int, na.rm = FALSE) {
	if (length(x) - sum(is.na(x)) > 1) {
		if (sd(x, na.rm = TRUE) > 0) {
			circ <- 2 * pi / int
			x_circ <- circular::circular(x * circ, type = "angles", units = "radians", rotation = "clock", modulo = "2pi")
			x_int <- circular::sd.circular(x_circ, na.rm = na.rm) / circ
			res <- as.numeric(x_int)
		} else {
			res <- 0
		}
	} else {
		res <- NA
	}
	
	res
}


#' Calculate tangential curvature.
#'
#' @param x A numeric vector of length 9, i.e., the 8 cells plus the center of a 3 x 3 square grid.
#' @param ... Not used, but required so that \code{tang_curvature} can be used as function argument to \code{\link{raster::focal}}
#' @return A numeric value in the range [‑1,1] with a sign that ensures that positive curvature equates to convex forms and negative curvature equates to concave forms
#' @examples
#' if (require(raster, quietly = TRUE)) {
#'     x <- raster::raster(matrix(1:36, ncol = 6), xmn=0, xmx=1, ymn=0, ymx=1)
#'     tcurv <- raster::focal(x, w = 3, fun = tang_curvature, pad = TRUE, padValue = NA)
#' }
tang_curvature <- function(x, ..., res_m) {
	#Zeverbergen, L. W., and C. R. Thorne. 1987. Quantitative Analysis of Land Surface Topography. Earth Surface Processes and Landforms 12: 47–56.
	#Schmidt, J., Evans, I.S. & Brinkmann, J. (2003) Comparison of polynomial models for land surface curvature calculation. International Journal of Geographical Information Science, 17, 797-814.
	#http://www.spatialanalysisonline.com/HTML/index.html?profiles_and_curvature.htm
	#http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00q90000000t000000
	stopifnot(length(x) == 9)
	D <- ((x[4] + x[6]) / 2 - x[5]) / res_m^2 # == fxx according to notation used by Schmidt et al.
	E <- ((x[2] + x[8]) / 2 - x[5]) / res_m^2 # == fyy
	F <- (-x[1] + x[3] + x[7] - x[9]) / (4 * res_m^2) # ?= fxy
	G <- (-x[4] + x[6]) / (2 * res_m^2) # == fx
	H <- (-x[8] + x[2]) / (2 * res_m^2) # == fy
	PLANC <- -2 * (D * H^2 + E * G^2 - F * G * H) / (G^2 + H^2) #eq. 18 for planform curvature (but see below that this is really tangential curvature): Zeverbergen et al.
	#TANGC <- - (D * H^2 - 2 * F * G * H + E * G^2) / ((G^2 + H^2) * sqrt(G^2 + H^2 + 1)) #own calculations based on Table 1 for tangential curvature: Schmidt et al.
	#I don't understand where the difference in the '2 *' between PLANC and TANGC comes from ...
	#However, "For example, Zevenbergen and Thorne (1987) used the term ‘planform curvature’ for tangential curvature, which is clearly not in plan view." (Schmidt et al.)
	PLANC	#multiply curvature values by ‑100 for convenience, as recommended by Zeverbergen and Thorne (1987), to give values in the approximate range [‑1,1] with a sign that ensures that positive curvature equates to convex forms and negative curvature equates to concave forms
}


rotate_coords <- function(coords, angle_rad = 0, center = c(0, 0)) {
	if (is.null(dim(coords))) coords <- matrix(unlist(coords), ncol = 2)
	
	c0 <- sweep(coords, MARGIN = 2, STATS = center, FUN = "-")
	
	res <- coords
	res[, 1] <- c0[, 1] * cos(angle_rad) - c0[, 2] * sin(angle_rad) + center[1]
	res[, 2] <- c0[, 1] * sin(angle_rad) + c0[, 2] * cos(angle_rad) + center[2]
	
	res
}


#' Determines the most common element(s)
#'
#' @param x A numeric vector.
#' @param use_mode A logical indicating whether to use \code{\link[modeest]{modeest::mfv}} if suitable.
#' @param return_single A logical indicating whether a single value should return even if several elements of \code{x} are the most common. If so, then one will be picked at random.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param seed An integer value or \code{NULL} to set the random number generator. \code{NULL} follows the behavior of \code{\link{set.seed}}.
#' @examples
#' x1 <- c(rep(NA, 3), 1:4)
#' x2 <- c(rep(1/3, 3), 1:4, rep(NA, 3))
#' majority(x1)
#' majority(x1, na.rm = TRUE)
#' majority(x1, return_single = FALSE, na.rm = TRUE)
#' majority(x1, use_mode = FALSE, return_single = FALSE, na.rm = TRUE)
#' majority(x1, use_mode = FALSE, return_single = TRUE, na.rm = TRUE)
#' majority(x1, use_mode = FALSE, return_single = FALSE, na.rm = FALSE)
#' majority(x2, return_single = FALSE)
majority <- function(x, use_mode = TRUE, return_single = TRUE, na.rm = FALSE, seed = NULL){
	ina <- is.na(x)
	
	if (all(ina)) {
		res <- NA
	} else {
		with_mode <- FALSE
		na_val <- NULL
		
		if (na.rm) {
			x <- na.exclude(x)
			with_mode <- TRUE
		}
		if (!na.rm && sum(ina) > 0 && !inherits(x, "integer")) {
			na_val <- min(x, na.rm = TRUE) - 100
			x[ina] <- na_val
			with_mode <- TRUE
		}
		with_mode <- if (use_mode) with_mode else FALSE 
		
		if (with_mode && requireNamespace("modeest", quietly = TRUE)) {
			res <- modeest::mfv(x)
		} else {
			temp <- table(x, useNA = if(na.rm) "no" else "always")
			res <- as.numeric(names(temp)[which(c(temp) == max(c(temp)))])
		}
		
		res[!is.finite(res)] <- NA
		if (!is.null(na_val) && is.finite(na_val)) {
			res[abs(res - na_val) <= sqrt(.Machine$double.eps)] <- NA
		}
		
		if (return_single && length(res) > 1) {
			if (!is.na(seed)) set.seed(seed)
			res <- res[sample(x = length(res), size = 1)]
		}
	}

	res
}



#---Danz et al. 2012: 'Spatial change in vegetation across the boundary'
sigmoidal <- function(x, b, c) {
	1 / (exp(-b * (x - c)) + 1)
}

sigmoidal_inv <- function(y, b, c) {
	c - 1 / b * log(1 / y - 1, base = exp(1))
}


m_glm <- function(family, data.) {
	w <- if (length(unique(data.[["y"]])) == 2 || identical(family[["family"]], "gaussian")) NULL else data.[["w"]]
	mfit <- try(stats::glm(y ~ x, data = data.[c("x", "y")], family = family, weights = w), silent = TRUE)
			
	if (!inherits(mfit, "try-error")) {
		mpred <- predict(mfit, newdata = data.[["newdata"]], type = "response", se.fit = TRUE)
		quals <- c(isConv = mfit$converged, edf = (temp <- extractAIC(mfit))[1], AIC = temp[2],
					logLik = logLik(mfit), deviance = deviance(mfit), df.resid = df.residual(mfit))
	} else {
		mpred <- NULL
		quals <- c(isConv = FALSE, edf = NA_integer_, AIC = NA_integer_, logLik = NA_integer_, deviance = NA_integer_, df.resid = NA_integer_)
	}
	
	list(m = mfit, preds = mpred, quals = quals)
}

m_glmm <- function(family, data.) {
	mfit <- try(lme4::glmer(y ~ x + (x|r), data = data.[c("x", "y", "r")], family = family), silent = TRUE)
	
	if (!inherits(mfit, "try-error")) {
		# unconditional (level-0 random effect) prediction
		mpred <- predict(mfit, newdata = data.[["newdata"]], re.form = ~ 0, type = "response")
		quals <- c(isConv = TRUE, edf = attr(logLik(mfit), "df"),
					lme4::llikAIC(mfit)[["AICtab"]][c("AIC", "logLik", "deviance", "df.resid")])
	} else {
		mpred <- NULL
		quals <- c(isConv = FALSE, edf = NA_integer_, AIC = NA_integer_, logLik = NA_integer_, deviance = NA_integer_, df.resid = NA_integer_)
	}
	
	list(m = mfit, preds = mpred, quals = quals)
}

m_sig <- function(data.) {
	mfit <- {i <- 1
				repeat {
					fit <- try(nls(y ~ sigmoidal(x, b, c), data = data., start = list(b = runif(1, -1, 1), c = runif(1, -1, 1))), silent = TRUE)
					if (i > 50 || !inherits(fit, "try-error")) break
					i <- i + 1
				}
			fit}
	
	if (!inherits(mfit, "try-error")) {
		mpred <- predict(mfit, newdata = data.[["newdata"]], type = "response")
		temp <- logLik(mfit)
		quals <- c(isConv = mfit[["convInfo"]][["isConv"]], edf = attr(temp, "df"), AIC = AIC(mfit),
					logLik = temp, deviance = deviance(mfit), df.resid = df.residual(mfit))
	} else {
		mpred <- NULL
		quals <- c(isConv = FALSE, edf = NA_integer_, AIC = NA_integer_, logLik = NA_integer_, deviance = NA_integer_, df.resid = NA_integer_)
	}

	list(m = mfit, preds = mpred, quals = quals)
}



#' Heaviside function.
#'
#' The piecewise constant Heaviside/unit step function.
#'
#' @param x A numeric (vector).
#' @param H0 Selects one of three versions of the Heaviside function. A value out of \code{{0, 0.5, 1}} corresponding to \code{Heaviside(0)}
#' @return A numeric vector of the same length as \code{x} with elements out of \code{{0, 1}}, respectively \code{{0, 0.5, 1}} if \code{H0 = 0.5}.
#' @examples
#' x <- (-2):2
#' cbind(x, Heaviside(x, H0 = 0))
#' cbind(x, Heaviside(x, H0 = 0.5))
#' cbind(x, Heaviside(x, H0 = 1))
#' \dontrun{cbind(x, Heaviside(x, H0 = 0.3))}
Heaviside <- function(x, H0 = 0L) {
	if (isTRUE(all.equal(H0, 0))) {
		# This appears to be the version used by Eppinga et al. 2013
		ifelse(x > 0, 1L, 0L)
	} else if (isTRUE(all.equal(H0, 1))) {
		ifelse(x >= 0, 1L, 0L)
	} else if (isTRUE(all.equal(H0, 0.5))) {
		# http://mathworld.wolfram.com/HeavisideStepFunction.html
		0.5 * (1 + sign(x))
	} else {
		stop(paste("The requested version of the 'Heaviside' function with H0 =", H0, "is not implemented."))
	}
}


#---Eppinga et al. 2013: eq. 17
transformation17 <- function(x) (2 * Heaviside(x) - 1) * log(1 + abs(x))

backtransformation17 <- function(y) (2 * Heaviside(y) - 1) * (exp(abs(y)) - 1)


erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


defuzz_2rad <- function(x) {
	x <- ifelse(x > -10*sqrt(.Machine$double.neg.eps) & x < 10*sqrt(.Machine$double.eps), 0, x)
	ifelse(x > 2*pi - 10*sqrt(.Machine$double.neg.eps) & x < 2*pi + 10*sqrt(.Machine$double.eps), 2*pi, x)
}

circle_2rad <- function(x){
	x <- ifelse(x <= 2*pi, x, x %% (2*pi + sqrt(.Machine$double.eps))) #I add sqrt(.Machine$double.eps) so that (x*2*pi) %% (2*pi) is equal to 2*pi
	ifelse(x >= 0, x, 2*pi + x)
}

rotate_azimuth_matrix <- function(mat, angle) {
	mat <- mat - angle	# positive angles means clockwise
	mat <- circle_2rad(mat)
	mat <- defuzz_2rad(mat)
	mat[mat == 0] <- 2*pi #previous step has direction 2*pi and not 0
	mat
}

transform_01_to_R <- function(x, ymin, ymax) ymin + x * (ymax - ymin)

transform_R_to_01 <- function(y, ymin, ymax) (y - ymin) / (ymax - ymin)


