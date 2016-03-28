# conversion of script to package
libraries  <- c("devtools")
l <- lapply(libraries, FUN=function(lib) stopifnot(require(lib, character.only=TRUE, quietly=FALSE)))

dir_dev <- "~/Dropbox/Work_Stuff/2_Research/Software/GitHub_Projects"
setwd(file.path(dir_dev, "ecotoner"))
devtools::document()
devtools::load_all()

if (FALSE) { # Vignettes
	#devtools::use_vignette("example")
	devtools::build_vignettes()
}

if (FALSE) { # Install package from github
	devtools::install_github("dschlaep/ecotoner")
}
