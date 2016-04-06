# demo
i	# 1
b	# 1
ic	# 3


dir_temp <- "~/Dropbox/Work_Stuff/2_Research/Software/GitHub_Projects"
dir_tools <- file.path(dir_temp, "ecotoner", "tools", "Location_EcotoneCriteria")

# Load/save data
ftemp <- file.path(dir_tools, "data_i1b1ic3_plot_ecotone_criteria.RData")
if (file.exists(ftemp)) {
	load(ftemp)
} else {
	save(esets, i, b, ic, temp_stline_cands, res_m,
		dens_Veg1ab, end1_toLeft, cv1, dens_Veg2ab, end2_toLeft, cv2, veg1main_leftmost, veg1max_rightmost, veg_density_extended_min, veg2max_leftmost, veg2main_rigthmost, veg_density_extended_min, croptemp, 
		distX1, distX2, elev, elev2, dens_Veg1, dens_BSE, dens_Veg2, dens_TF, start1_m, end1_m, start2_m, end2_m,
		file = ftemp)
}



# Panels a and b from objects of function 'determine_ecotone_linear_extent' for 'limits1_ZoneEcolBoundary'
# Panel c same plot as by function 'plot_transect_density_profiles'

add_panel <- function(distX, dens_Veg, end_toLeft, cv, col) {
	ldens <- zoo::rollmean(dens_Veg, k = 3, fill = "extend", align = "center")
	plot(distX, ldens, type = "l", xlab = "", ylab = "Vegetation density (rel.)", ylim = c(0, 1), lwd = 2, col = col, axes = FALSE)
	axis(side = 1, pos = 0, labels = FALSE)
	axis(side = 2, pos = 0)
	abline(h = c(veg_density_low, veg_density_extended_min, veg_density_high), col = c("darkgray", "green", "blue"))

	if (end_toLeft) {
		with(cv, rect(res_m * ezone[1], 0, res_m * ezone[length(ezone)], veg_density_high, col = adjustcolor("blue", alpha = 0.3), lty = 2))
		with(cv, rect(res_m * szone[1], 0, res_m * szone[length(szone)], veg_density_low, col = adjustcolor("darkgray", alpha = 0.3), lty = 2))
	} else {
		with(cv, rect(res_m * ezone[1], 0, res_m * ezone[length(ezone)], veg_density_low, col = adjustcolor("darkgray", alpha = 0.3), lty = 2))
		with(cv, rect(res_m * szone[1], 0, res_m * szone[length(szone)], veg_density_high, col = adjustcolor("blue", alpha = 0.3), lty = 2))
	}
}		

pdf(height = 7.1, width = 4, file = file.path(dir_tools, "Fig_EcotoneCriteria.pdf"))
op_old <- par(mfrow = c(3, 1), mar = c(0.75, 2.5, 1, 2.1), mgp = c(1.25, 0.5, 0), cex = cex <- 1)

# Panel a (Veg1)
add_panel(distX1, dens_Veg1ab, end1_toLeft, cv1, col = "orange")
with(cv1, rect(res_m * veg1main_leftmost, 0, res_m * veg1max_rightmost, veg_density_extended_min, col = adjustcolor("green", alpha = 0.3), lty = 2))
abline(v = res_m * croptemp, lty = 3, lwd = 3, col = "black")
mtext(side = 3, text = "(a)", font = 2, adj = 0.045, line = -0.25)

# Panel b (Veg2)
par(mar = c(1.25, 2.5, 0.5, 2.1))
add_panel(distX1, dens_Veg2ab, end2_toLeft, cv2, col = "darkgreen")
with(cv2, rect(res_m * veg2max_leftmost, 0, res_m * veg2main_rigthmost, veg_density_extended_min, col = adjustcolor("green", alpha = 0.3), lty = 2))
abline(v = res_m * croptemp, lty = 3, lwd = 3, col = "black")
mtext(side = 3, text = "(b)", font = 2, adj = 0.045, line = -0.25)

# Panel c (all)
par(mar = c(2.75, 2.5, 0.1, 2.1))
panel_transect_density_profiles(distX1, distX2, elev, dens_Veg1, dens_BSE, dens_Veg2, dens_TF, start1_m, end1_m, start2_m, end2_m, cex)
matlines(distX1, t(raster::as.matrix(elev2)), col = adjustcolor("black", alpha.f = 0.05))
lines(distX1, elev, lwd = 2, col = "black")
mtext(side = 3, text = "(c)", font = 2, adj = 0.045, line = -0.25)

options(op_old)
dev.off()
