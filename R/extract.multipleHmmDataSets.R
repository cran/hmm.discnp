"[.multipleHmmDataSets" <- function(x,i) {
   cls <- class(x)
   class(x) <- "list"
   y <- x[i]
   class(y) <- cls
   y
}
