# -----------------------------------
# Pre-defined function for statistics
# -----------------------------------
twosidePVttestEX <- function(data1,data2){
#' Two-sided directional t-test returning a signed p-value
#'
#' @param data1 Numeric vector for group 1
#' @param data2 Numeric vector for group 2
#' @return Signed p-value: positive if mean(data1) > mean(data2), negative otherwise
#' @examples
#' twosidePVttestDIR(rnorm(10, mean=5), rnorm(10, mean=3))
  if(length(data1) == 0 | length(data2) == 0){
    return(NA)
  }else{
    m1 <- mean(data1,na.rm=T)
    m2 <- mean(data2,na.rm=T)
    if(m1 > m2){
      TEST <- t.test(data1,data2,alternative="greater")
      PV <-  TEST$p.value
      score <-  TEST$statistic
    }else{
      TEST <- t.test(data1,data2,alternative="less")
      PV <-  TEST$p.value
      score <-  TEST$statistic
    }
    return(c(PV, score, m1, m2))
  }
}

twosidePVttestDIR <- function(data1,data2){
#' Two-sided directional t-test returning a signed p-value
#'
#' @param data1 Numeric vector for group 1
#' @param data2 Numeric vector for group 2
#' @return Signed p-value: positive if mean(data1) > mean(data2), negative otherwise
#' @examples
#' twosidePVttestDIR(rnorm(10, mean=5), rnorm(10, mean=3))
  if(length(data1) == 0 | length(data2) == 0){
    return(NA)
  }else{
    m1 <- mean(data1,na.rm=T)
    m2 <- mean(data2,na.rm=T)
    if(m1 > m2){
      TEST <- t.test(data1,data2,alternative="greater")
      PV <-  TEST$p.value
    }else{
      TEST <- t.test(data1,data2,alternative="less")
      PV <-  -TEST$p.value
    }
    return(PV)
  }
}
