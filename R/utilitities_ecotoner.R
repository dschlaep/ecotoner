#' @export
simplify2result <- function(x) {
	len0 <- sapply(x, length) == 0
	if (sum(len0) > 0) x <- x[!len0]
	
    dim1 <- if (inherits(x[[1]], "data.frame") || inherits(x[[1]], "matrix")) {
				dim(x[[1]])
			} else {
				c(1, length(x[[1]]))
			}
	rows1 <- seq_len(dim1[1])
	res <- as.data.frame(matrix(NA, nrow = length(x) * dim1[1], ncol = dim1[2], dimnames = list(NULL, colnames(x[[1]]))))
	for (i in seq_along(x)) res[(i - 1) * dim1[1] + rows1, ] <- x[[i]]
	
	temp <- warnings()
	if (length(temp) > 0) print(temp)
	
	res
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
N_of_measure_calls <- function(ecotoner_settings) transect_N(ecotoner_settings) * neighbors_N(ecotoner_settings) * length(get("migtypes", envir = etr_vars))
