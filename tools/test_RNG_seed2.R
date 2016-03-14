N1 <- 1e4
N2 <- 1e4
Nb <- 1e1
bi <- seq_len(Nb)

cl <- parallel::makeCluster(20, type = "FORK")

resT <- parallel::parSapply(cl, seq_len(N1), function(i) {set.seed(i); sample(N2)})

blocks <- lapply(seq_len(floor(N2 / Nb)), function(i) (i - 1) * Nb + bi)

dups <- parallel::parSapply(cl, seq_along(blocks), function(i, blocks, resT) {temp <- t(resT[blocks[[i]], ]); anyDuplicated(temp)}, blocks, resT)

print(sum(dups))

parallel::stopCluster(cl)
