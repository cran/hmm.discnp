ordinal <- function(k) paste(k, ordinalsuffix(k), sep="")

ordinalsuffix <- function(k) {
  last <- abs(k) %% 10
  lasttwo <- abs(k) %% 100
  isteen <- (lasttwo > 10 & lasttwo < 20)
  ending <- ifelse(isteen, "th",
                   ifelse(last == 1, "st",
                          ifelse(last == 2, "nd",
                                 ifelse(last == 3, "rd",
                                        "th"))))
  return(ending)
}
