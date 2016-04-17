#' Simplify a list to a data.frame
#'
#' Convert a list of elements with the same number of columns to a data.frame (discarding NULL elements) and return it unchanged otherwise.
#'
#' @param x A R object. See \code{Value}.
#' @param showWarnings A logical value. If \code{TRUE} then warnings are collected and printed if there are any.
#'
#' @return A data.frame if x is a list of elements with the same number of columns. If the list elements do not all have the same number of columns, then a list without the NULL elements is returned. If x is not a list, then x is returned unchanged.
#' @seealso \code{\link{simplify2array}}
#' @export
simplify2result <- function(x, showWarnings = TRUE) {
	if (inherits(x, "list")) {
		len0 <- sapply(x, length) == 0
		if (any(len0)) x <- x[!len0]
	
		dims <- sapply(x, function(d) if (inherits(d, "data.frame") || inherits(d, "matrix")) dim(d) else c(1, length(d)))
		if (all(diff(dims[2, ]) == 0)) { # all must have same number of columns
			cumrows <- c(0, cumsum(dims[1, ]))
			res <- as.data.frame(matrix(NA, 
						nrow = cumrows[1 + length(x)],
						ncol = dims[2, 1],
						dimnames = list(NULL, colnames(x[[1]]))))
			for (i in seq_along(x)) res[cumrows[i] + seq_len(dims[1, i]), ] <- x[[i]]
			x <- res
		}
	}
	
	if (showWarnings) {
		temp <- warnings()
		if (length(temp) > 0) {
			print(temp)
			# TODO: according to ?warnings it is undocumented and subject to change where last.warning is stored; yet, without flushing, we carry forward old warnings - until R has a flushWarnings function, we meddle with 'baseenv'
			try(assign("last.warning", NULL, envir = baseenv()))
		}
	}
		
	x
}


#' A load balancing version of clusterApply with increased level of balancing
#'
#' This function is a modified copy of \code{\link{parallel::parLapplyLB}}  (v3.2.4)
#' 
#' Bug report is submitted under \href{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=16792}{# 16792}
#'
#' @param cl A cluster object, created by this package or by package \href{https://CRAN.R-project.org/package=snow}{\pkg{snow}}. If \code{NULL}, use the registered default cluster.
#' @param X A vector (atomic or list).
#' @param fun A function or character string naming a function.
#' @param ... Additional arguments to pass to \code{fun}; beware of partial matching to earlier arguments.
#' @param chunksize An integer value. A suggestion on how large the number of tasks should be that are sent per job to a worker.
#'
#' @seealso \code{\link{parallel::parLapplyLB}} and \code{\link{parallel::parLapply}}
#' @return A list of the same length as \code{X} each with the returned value of \code{fun}.
#' @export
parLapplyLBc <- function (cl = NULL, X, fun, ..., chunksize = 10L) {
	# modified from parallel::parLapplyLB  (v3.2.4)
    cl <- parallel:::defaultCluster(cl)
    
    # x <- splitList(X, length(cl)) # this is the original code, but if length(x) == length(cl), then 'dynamicClusterApply' does not load balance the tasks
    # instead the modified version creates a list of indices much longer than there are workers
    chunksize <- as.integer(chunksize)
    x <- if (chunksize == 1L) X else parallel:::splitList(X, max(length(cl), round(length(X) / chunksize)))
    
    do.call(c, parallel::clusterApplyLB(cl, x = x, fun = lapply, fun, ...), quote = TRUE)
}

#' Create a function to apply a function to each transect either in parallel with load balancing or in a serial process depending on the core settings
#'
#' @param ecotoner_settings An object of class 'EcotonerSettings' of which the cores_N(.) will be used.
#' @return A function
#' @export
fun_to_apply_foreach_transect <- function(ecotoner_settings) {
	if (cores_N(ecotoner_settings) > 1)
		# load balancing is reproducible if each transect sets its own unique random seed
		function(X, FUN, ...) simplify2result(parLapplyLBc(cl, X, FUN, ...))
	else
		function(X, FUN, ...) simplify2result(lapply(X, FUN, ...))

}

#' Creates a path
#' 
#' This function first checks if the path exists. If the path does not exist, it attempts to create the path using
#' \code{\link{dir.create}} with different default values.
#'
#' @inheritParams base::dir.create
#'
#' @return The function returns invisible \code{TRUE} if the path exists
#' or if the path could successfully (and recursively by default) be created. If the path does not exist and the path cannot
#' be created it throws an error with a hopefully useful message.
#'
#' @seealso \code{\link{dir.create}}
#' @export
dir_create <- function(path, showWarnings = FALSE, recursive = TRUE, mode = "0777") {
	if (!dir.exists(path) && !dir.create(path, recursive = recursive, showWarnings = showWarnings, mode = mode)) {
		temp <- sys.calls()
		nframe <- sys.nframe() - 1
		tcall <- if (length(temp) < nframe || nframe < 1) temp[[1]] else temp[[nframe]]
		stop(deparse(tcall), ": failed to create directory ", sQuote(path), call. = FALSE)
	}
	invisible(TRUE)
}

#' @export
flag_itransect <- function(ecotoner_settings, iflag) {
	paste0("Transect", formatC(iflag, format = "d", flag = "0", width = floor(log(transect_N(ecotoner_settings), base = 10) + 1)))
}

#' @export
flag_basename <- function(ecotoner_settings, iflag, b) {
	paste0("Fig_", iflag, "_Neighbors", neighborhoods(ecotoner_settings)[b], "_")
}

#' @export
fname_etfailed <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_failed.RData"))
}

#' @export
fname_etlocated <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_etransect.RData"))
}

#' @export
fname_etnone <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_notransect.RData"))
}

#' @export
fname_etsearching <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_etsearching.RData"))
}

#' @export
fname_etmeasured <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_et_measured.RData"))
}

#' @export
write_ecotoner_row <- function(data_row, filename, tag_fun = "", tag_id = "") {
	res <- if (all(apply(data_row, 1:2, class) %in% c("numeric", "logical", "integer", "character"))) {
				appT <- file.exists(filename)
				write.table(data_row, file = filename, append = appT, sep = ",", dec = ".", qmethod = "double", row.names = FALSE, col.names = !appT)
				"written"
			} else {
				warning(tag_fun, ": ", tag_id, " contains unsuitable elements, e.g., a 'list'", immediate. = TRUE)
				"error"
			}

	invisible(res)
}

#' @export
N_of_location_calls <- function(ecotoner_settings) transect_N(ecotoner_settings) * neighbors_N(ecotoner_settings)

#' @export
N_of_measure_calls <- function(ecotoner_settings) N_of_location_calls(ecotoner_settings) * length(get("migtypes", envir = etr_vars))

#' @export
add_new_timing <- function(..., time_h, filename) {
	x_cols <- matrix(c(as.numeric(unlist(list(...))), round(time_h, 2)), nrow = 1)
	appT <- file.exists(filename)
	write.table(x_cols, file = filename, append = appT, sep = ",", dec = ".", qmethod = "double", row.names = FALSE, col.names = !appT)
}

#' @export
get_max_timing <- function(filename, add_hours = 2) {
	add_hours + if (file.exists(filename)) {
					x <- read.csv(filename, header = TRUE, row.names = NULL)
					x <- x[, ncol(x)]
					max(x, na.rm = TRUE)
				} else 0
}
