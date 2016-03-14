library(parallel)

seed <- 123

runs <- list(sample(10, 10), sample(10, 10), sample(10, 10))
print(runs)

foo <- function(x, seed = NULL) {
	if (is.null(seed)) Sys.sleep(runif(1, min=0.1, max=30))

	print("foo: 1"); print(str(.Random.seed))
	if (!is.na(seed)) set.seed(seed)
	res1 <- sample(10, 3)
	
	print("foo: 2"); print(str(.Random.seed))
	if (!is.na(seed)) set.seed(seed)
	print(res <- c(res1, sample(10, 3)))
	res
}


cl <- makeCluster(10, type = "FORK", outfile = "log_seed.txt")
RNGkind("L'Ecuyer-CMRG")


print("foo1, set seed")
clusterSetRNGStream(cl, seed)
print(RNGkind()); print(str(rseed <- .Random.seed))
res <- parSapplyLB(cl, runs[[1]], foo, seed = rseed)
print(RNGkind()); print(str(rseed <- .Random.seed))
print(t(res))

print("foo1, no seed")
clusterSetRNGStream(cl, seed)
print(RNGkind()); print(str(rseed <- .Random.seed))
res1 <- parSapplyLB(cl, runs[[2]], foo, seed = NULL)
print(RNGkind()); print(str(rseed <- .Random.seed))
print(t(res1))

print("foo1, no seed: repeat")
clusterSetRNGStream(cl, seed)
res2 <- parSapplyLB(cl, runs[[3]], foo, seed = NULL)
identical(res1, res2)
print(t(res2))

