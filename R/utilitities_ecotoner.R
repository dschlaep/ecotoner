
flag_itransect <- function(i, ecotoner_settings) {
	paste0("Transect", formatC(i, format = "d", flag = "0", width = floor(log(transect_N(ecotoner_settings), base = 10) + 1)))
}

flag_basename <- function(iflag, b, ecotoner_settings) {
	paste0("Fig_", iflag, "_Neighbors", neighborhoods(ecotoner_settings)[b], "_")
}

fname_etfailed <- function(iflag, ecotoner_settings) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_notransect.RData"))
}

fname_etlocated <- function(iflag, ecotoner_settings) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_etransect.RData"))
}

fname_etmeasured <- function(iflag, ecotoner_settings) {
	file.path(dir_out_dat(ecotoner_settings), paste0(iflag, "_et_measured.RData"))
}
