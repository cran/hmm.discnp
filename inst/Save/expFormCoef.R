expFormCoef <- function(pea) {
   pea <- pea/sum(pea)
   iz  <- pea == 0
   pnz <- pea[!iz]
   m   <- length(pnz)
   M   <- matrix(pnz[-1],nrow=m-1,ncol=m-1) - diag(m-1)
   v   <- c(0,log(solve(M,-pnz[-1])))
   u   <- rep(-Inf,length(pea))
   u[!iz] <- v
   u
}
