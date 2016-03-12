
flag_itransect <- function(ecotoner_settings, iflag) {
	paste0("Transect", formatC(iflag, format = "d", flag = "0", width = floor(log(transect_N(ecotoner_settings), base = 10) + 1)))
}

flag_basename <- function(ecotoner_settings, iflag, b) {
	paste0("Fig_", iflag, "_Neighbors", neighborhoods(ecotoner_settings)[b], "_")
}

fname_etfailed <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_failed.RData"))
}

fname_etlocated <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_etransect.RData"))
}


fname_etnone <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_notransect.RData"))
}


fname_etsearching <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_etsearching.RData"))
}

fname_etmeasured <- function(ecotoner_settings, iflag) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_et_measured.RData"))
}
