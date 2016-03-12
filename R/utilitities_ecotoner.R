#' @export
simplify2result <- function(x) {
	len0 <- sapply(x, length) == 0
	if (sum(len0) > 0) x <- x[!len0]
	
	dim1 <- dim(x[[1]])
	rows1 <- seq_len(dim1[1])
	res <- as.data.frame(matrix(NA, nrow = length(x) * dim1[1], ncol = dim1[2], dimnames = list(NULL, colnames(x[[1]]))))
	for (i in seq_along(x)) res[(i - 1) * dim1[1] + rows1, ] <- x[[i]]
	
	temp <- warnings()
	if (length(temp) > 0) print(temp)
	
	res
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
