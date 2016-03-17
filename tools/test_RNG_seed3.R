rm(list=ls(all=TRUE))

f1 <- function(x) {
	kind10 <- RNGkind()[1]
	if (exists(".Random.seed")) {
		rs10 <- .Random.seed
		rm(.Random.seed)
	} else rs10 <- NULL
	set.seed(x)
	kind11 <- RNGkind()[1]
	sample(1)
	rs11 <- .Random.seed
	if (exists(".Random.seed")) rm(.Random.seed)
	
	
	f2 <- function(x) {
		kind0 <- RNGkind()[1]
		if (exists(".Random.seed")) {
			rs0 <- .Random.seed
			rm(.Random.seed)
		} else rs0 <- NULL
		RNGkind(kind = "L'Ecuyer-CMRG")
		set.seed(x)
		kind1 <- RNGkind()[1]
		sample(1)
		rs1 <- .Random.seed
		if (exists(".Random.seed")) rm(.Random.seed)
	
		list(rseed0 = rs0, kind0 = kind0, rseed1 = rs1, kind1 = kind1)
	}
	f2res <- f2(x)

	kind20 <- RNGkind()[1]
	if (exists(".Random.seed")) {
		rs20 <- .Random.seed
		rm(.Random.seed)
	} else rs20 <- NULL
	RNGkind(kind = "Wichmann-Hill")
	set.seed(x)
	kind21 <- RNGkind()[1]
	sample(1)
	rs21 <- .Random.seed
	if (exists(".Random.seed")) rm(.Random.seed)
	
	list(f1a = list(rseed0 = rs10, kind0 = kind10, rseed1 = rs11, kind1 = kind11),
		 f2 = f2res,
		 f1b = list(rseed0 = rs20, kind0 = kind20, rseed1 = rs21, kind1 = kind21))
}

kinds <- list()
rs <- list()

RNGkind(kind = "Marsaglia-Multicarry")
set.seed(5)
kinds[[1]] <- RNGkind()[1]
sample(1)
rs[[1]] <- .Random.seed

res1 <- f1(5)
kinds[[2]] <- RNGkind()[1]
rs[[2]] <- .Random.seed

res2 <- f1(5)
kinds[[3]] <- RNGkind()[1]
rs[[3]] <- .Random.seed
RNGkind("default")
kinds[[4]] <- RNGkind()
rs[[4]] <- .Random.seed
RNGkind(kind = kinds[[3]][1])
kinds[[5]] <- RNGkind()[1]
rs[[5]] <- .Random.seed
.Random.seed <- rs[[1]]
kinds[[6]] <- RNGkind()[1]
rs[[6]] <- .Random.seed

temp <- sapply(kinds, print)
	#[1] "Marsaglia-Multicarry"
	#[1] "Wichmann-Hill"
	#[1] "Wichmann-Hill"
	#[1] "Mersenne-Twister" "Inversion"       
	#[1] "Wichmann-Hill"
	#[1] "Marsaglia-Multicarry"

temp <- sapply(rs, function(x) str(x))
	# int [1:3] 401 1686395561 149529984
	# int [1:4] 400 16752 21308 2998
	# int [1:4] 400 16752 21308 2998
	# int [1:626] 403 624 1767525564 1072758413 1870004202 1063702627 -833821112 26087337 -2059985066 -1726908961 ...
	# int [1:4] 400 13587 11201 18402
	# int [1:3] 401 1686395561 149529984

unlist(res1)[grepl("kind", names(unlist(res1)))]
	#             f1a.kind0              f1a.kind1               f2.kind0 
	#"Marsaglia-Multicarry" "Marsaglia-Multicarry" "Marsaglia-Multicarry" 
	#              f2.kind1              f1b.kind0              f1b.kind1 
	#       "L'Ecuyer-CMRG"        "L'Ecuyer-CMRG"        "Wichmann-Hill" 
unlist(res2)[grepl("kind", names(unlist(res2)))]
	#      f1a.kind0       f1a.kind1        f2.kind0        f2.kind1       f1b.kind0 
	#"Wichmann-Hill" "Wichmann-Hill" "Wichmann-Hill" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" 
	#      f1b.kind1 
	#"Wichmann-Hill" 


# ==> setting RNGkind(kind = "XXX") updates .Random.seed to match kind, but this update is not reversible, i.e., switching from kind = X1 to X2 and back to X1 does not get back to the same .Random.seed

# ==> assigning an appropriate value to .Random.seed updates RNGkind() accordingly

# ==> setting RNGkind in any function carries forward to child-, parent-, and the global environment




#########################################
#########################################
library(parallel)
setwd("~/Dropbox/Work_Stuff/2_Research/Software/GitHub_Projects/ecotoner/tools")
cl <- parallel::makeCluster(2, type = "FORK", outfile = "RNG_test3.txt")
#cl <- parallel::makeCluster(2, type = "SOCK", outfile = "RNG_test3.txt")

RNGkind(kind = "default")
kinds <- list()
rs <- list()
set.seed(5)
kinds[[1]] <- RNGkind()[1]
sample(1)
rs[[1]] <- .Random.seed

res3 <- parallel::parLapply(cl, 1:4, f1)
names(res3) <- paste0("run", 1:4)
parallel::stopCluster(cl)

kinds[[2]] <- RNGkind()[1]
sample(1)
rs[[2]] <- .Random.seed


temp <- sapply(kinds, print)
	#[1] "Mersenne-Twister"    
	#[1] "Mersenne-Twister"   

temp <- sapply(rs, function(x) str(x))
	# int [1:626] 403 1 -1802397187 26112759 1167250462 -480616744 -1319230917 2098946009 -1500478104 630330518 ...
	# int [1:626] 403 2 -1802397187 26112759 1167250462 -480616744 -1319230917 2098946009 -1500478104 630330518 ...
identical(rs[[1]][-2], rs[[2]][-2])
	# TRUE

unlist(res3)[grepl("kind", names(unlist(res3)))]
	#        run1.f1a.kind0         run1.f1a.kind1          run1.f2.kind0 
	#"Marsaglia-Multicarry" "Marsaglia-Multicarry" "Marsaglia-Multicarry" 
	#         run1.f2.kind1         run1.f1b.kind0         run1.f1b.kind1 
	#       "L'Ecuyer-CMRG"        "L'Ecuyer-CMRG"        "Wichmann-Hill" 

	#        run2.f1a.kind0         run2.f1a.kind1          run2.f2.kind0 
	#       "Wichmann-Hill"        "Wichmann-Hill"        "Wichmann-Hill" 
	#         run2.f2.kind1         run2.f1b.kind0         run2.f1b.kind1 
	#       "L'Ecuyer-CMRG"        "L'Ecuyer-CMRG"        "Wichmann-Hill" 

	#        run3.f1a.kind0         run3.f1a.kind1          run3.f2.kind0 
	#"Marsaglia-Multicarry" "Marsaglia-Multicarry" "Marsaglia-Multicarry" 
	#         run3.f2.kind1         run3.f1b.kind0         run3.f1b.kind1 
	#       "L'Ecuyer-CMRG"        "L'Ecuyer-CMRG"        "Wichmann-Hill" 

	#        run4.f1a.kind0         run4.f1a.kind1          run4.f2.kind0 
	#       "Wichmann-Hill"        "Wichmann-Hill"        "Wichmann-Hill" 
	#         run4.f2.kind1         run4.f1b.kind0         run4.f1b.kind1 
	#       "L'Ecuyer-CMRG"        "L'Ecuyer-CMRG"        "Wichmann-Hill" 


# ==> RNGkind() on the cluster at execution time does not and is not affected by the global environment's RNGkind()

# ==> RNGkind() on each cluster remains to the next function call

# ==> the initial state of RNGkind() on the cluster
#		- cluster type = "FORK: then it depends on what RNGkind() was in the environment when the cluster was created
#		- cluster type = "SOCK": then it will be the default value, i.e., "Mersenne-Twister"






#########################################
#########################################
library(parallel)
setwd("~/Dropbox/Work_Stuff/2_Research/Software/GitHub_Projects/ecotoner/tools")
kinds <- list()
rs <- list()

f2 <- function(i, pseeds) {
	set_RNG_stream(pseeds[[1 + (1 + i) %% 2]])
	f1(i)
}

RNGkind(kind = "Marsaglia-Multicarry")
set.seed(5)
kinds[[1]] <- RNGkind()[1]
sample(1)
rs[[1]] <- .Random.seed

#cl <- parallel::makeCluster(2, type = "FORK", outfile = "RNG_test3.txt")
cl <- parallel::makeCluster(2, type = "SOCK", outfile = "RNG_test3.txt")
clusterExport(cl, varlist = c("f2", "f1", "set_RNG_stream")) # only necessary for SOCKET-clusters

RNGkind(kind = "default")
set.seed(5)
kinds[[2]] <- RNGkind()[1]
sample(1)
rs[[2]] <- .Random.seed

pseeds <- prepare_RNG_streams(N = 2, iseed = 5)
res5 <- parallel::parLapply(cl, 1:4, f2, pseeds)
names(res5) <- paste0("run", 1:4)
parallel::stopCluster(cl)

kinds[[3]] <- RNGkind()[1]
sample(1)
rs[[3]] <- .Random.seed


temp <- sapply(kinds, print)
	#[1] "Marsaglia-Multicarry"
	#[1] "Mersenne-Twister"    
	#[1] "Mersenne-Twister"   

temp <- sapply(rs, function(x) str(x))
	# int [1:3] 401 1686395561 149529984
	# int [1:626] 403 1 -1802397187 26112759 1167250462 -480616744 -1319230917 2098946009 -1500478104 630330518 ...
	# int [1:626] 403 2 -1802397187 26112759 1167250462 -480616744 -1319230917 2098946009 -1500478104 630330518 ...
identical(rs[[2]][-2], rs[[3]][-2])
	# TRUE

unlist(res5)[grepl("kind", names(unlist(res5)))]
	# run1.f1a.kind0  run1.f1a.kind1   run1.f2.kind0   run1.f2.kind1  run1.f1b.kind0 
	#"L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" 
	# run1.f1b.kind1  run2.f1a.kind0  run2.f1a.kind1   run2.f2.kind0   run2.f2.kind1 
	#"Wichmann-Hill" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" 
	# run2.f1b.kind0  run2.f1b.kind1  run3.f1a.kind0  run3.f1a.kind1   run3.f2.kind0 
	#"L'Ecuyer-CMRG" "Wichmann-Hill" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" 
	#  run3.f2.kind1  run3.f1b.kind0  run3.f1b.kind1  run4.f1a.kind0  run4.f1a.kind1 
	#"L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "Wichmann-Hill" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" 
	#  run4.f2.kind0   run4.f2.kind1  run4.f1b.kind0  run4.f1b.kind1 
	#"L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "L'Ecuyer-CMRG" "Wichmann-Hill" 

sapply(res5, function(x) sapply(x, function(y) y$rseed0))
	#             run1        run2        run3        run4
	# [1,]         407         407         407         407
	# [2,]  1157214768 -1828978842  1157214768 -1828978842
	# [3,] -1674567567   -44835162 -1674567567   -44835162
	# [4,] -1532971138  -788600610 -1532971138  -788600610
	# [5,] -1249749529   520285260 -1249749529   520285260
	# [6,]  1302496508   327419947  1302496508   327419947
	# [7,]  -253670963  2055199760  -253670963  2055199760
	# [8,]         407         407         407         407
	# [9,]  -169270483 -1619336578  1225564623  -224501472
	#[10,]  -442010614  -714750745  -987490876 -1260231007
	#[11,]   471000765  2099518112  -397489685  1524376091
	#[12,]  -222347416   158863565   540074546   921285527
	#[13,]  1489374793 -1093870294   617851915 -1965393172
	#[14,]  1855986652  1639909500 -2043818655  2035071489
	#[15,]         407         407         407         407
	#[16,]  -169270483 -1619336578  1225564623  -224501472
	#[17,]  -442010614  -714750745  -987490876 -1260231007
	#[18,]   471000765  2099518112  -397489685  1524376091
	#[19,]  -222347416   158863565   540074546   921285527
	#[20,]  1489374793 -1093870294   617851915 -1965393172
	#[21,]  1855986652  1639909500 -2043818655  2035071489


# ==> function 'prepare_RNG_streams' in combination with 'set_RNG_stream' works as designed both with FORK- and SOCKET-clusters



#########################################
#########################################
RNGkind(kind = "Marsaglia-Multicarry")
str(.Random.seed)
set.seed(5, kind = "L'Ecuyer-CMRG")
RNGkind()
str(.Random.seed)
