#Useful fucntions to decompose life expectancy using Horiuchi's method


#labesls for causes of death
cause_names<-c("1"="Circulatory diseases", "2"="Neoplasm","3"="Diabetes",
               "4"="Homicide and other violent causes","5"="Other external causes",
               "6"= "Conditions of the perinatal period","7"="Respiratory diseases",
               "8"="Infectious diseases","9"="Digestive diseases","10"="Others")

#labels for age-groups 
age_names<-c("0"="0","1"="1-4","5"="5-9","10"="10-14","15"="15-19","20"="20-24",
             "25"="25-29","30"="30-34","35"="35-39","40"="40-44","45"="45-49",
             "50"="50-54","55"="55-59","60"="60-64","65"="65-69","70"="70-74",
             "75"="75-79","80"="80-84","85"="85+")     


#we need a function for the standard deviation of the age at death distribution and life expectancy
e0.frommx <- function(nmx =  mx, sex=1, age = c(0, 1, seq(5, 85, 5)), nax = NULL){
  n   <- c(diff(age), 999)
  
  if (is.null(nax)) {
    nax <- 0.5 * n
    if (n[2] == 4) {
      if (sex == 1) {
        if (nmx[1] >= 0.107) {
          nax[1] <- 0.33
          nax[2] <- 1.352
        }
        else {
          nax[1] <- 0.045 + 2.684 * nmx[1]
          nax[2] <- 1.651 - 2.816 * nmx[1]
        }
      }
      if (sex == 2) {
        if (nmx[1] >= 0.107) {
          nax[1] <- 0.35
          nax[2] <- 1.361
        }
        else {
          nax[1] <- 0.053 + 2.8 * nmx[1]
          nax[2] <- 1.522 - 1.518 * nmx[1]
        }
      }
    }
  }
  nqx          <- (n * nmx)/(1 + (n - nax) * nmx)
  nqx          <- c(nqx[-(length(nqx))], 1)
  nqx[nqx > 1] <- 1
  
  npx <- 1 - nqx
  lx <- cumprod(c(1, npx))
  ndx <- -diff(lx)
  lxpn <- lx[-1]
  nLxpn <- n * lxpn + ndx * nax
  nLx <- c(nLxpn[-length(nLxpn)], lxpn[length(lxpn)-1]/nmx[length(nmx)])
  Tx <- rev(cumsum(rev(nLx)))
  lx <- lx[1:length(age)]
  ex <- Tx/lx
  e0 <- ex[1]
  
  return(e0)
}


Decomp <-function (func, rates1, rates2, N, ...) {
  y1 <- func(rates1, ...)
  y2 <- func(rates2, ...)
  d <- rates2 - rates1
  n <- length(rates1)
  delta <- d/N
  x <- rates1 + d * matrix(rep(0.5:(N - 0.5)/N, length(rates1)), 
                           byrow = TRUE, ncol = N)
  cc <- matrix(0, nrow = n, ncol = N)
  for (j in 1:N) {
    for (i in 1:n) {
      z <- rep(0, n)
      z[i] <- delta[i]/2
      cc[i, j] <- func((x[, j] + z), ...) - func((x[, j] - 
                                                    z), ...)
    }
  }
  return(rowSums(cc))
}


sd.frommx <- function(nmx =  mx, sex=1, age = c(0, 1, seq(5, 85, 5)), nax = NULL){
  n   <- c(diff(age), 999)
  
  if (is.null(nax)) {
    nax <- 0.5 * n
    if (n[2] == 4) {
      if (sex == 1) {
        if (nmx[1] >= 0.107) {
          nax[1] <- 0.33
          nax[2] <- 1.352
        }
        else {
          nax[1] <- 0.045 + 2.684 * nmx[1]
          nax[2] <- 1.651 - 2.816 * nmx[1]
        }
      }
      if (sex == 2) {
        if (nmx[1] >= 0.107) {
          nax[1] <- 0.35
          nax[2] <- 1.361
        }
        else {
          nax[1] <- 0.053 + 2.8 * nmx[1]
          nax[2] <- 1.522 - 1.518 * nmx[1]
        }
      }
    }
  }
  nqx          <- (n * nmx)/(1 + (n - nax) * nmx)
  nqx          <- c(nqx[-(length(nqx))], 1)
  nqx[nqx > 1] <- 1
  
  npx <- 1 - nqx
  lx <- cumprod(c(1, npx))
  ndx <- -diff(lx)
  lxpn <- lx[-1]
  nLxpn <- n * lxpn + ndx * nax
  nLx <- c(nLxpn[-length(nLxpn)], lxpn[length(lxpn)-1]/nmx[length(nmx)])
  Tx <- rev(cumsum(rev(nLx)))
  lx <- lx[1:length(age)]
  ex <- Tx/lx
  nax[length(nax)] <- ex[length(ex)]
  vx <- sum(ndx*(age+nax-ex[1L])^2)
  sd <- sqrt(vx)
  return(sd)
}

sdfrommxc <- function(mxcvec,sex=1){
  dim(mxcvec) <- c(19,length(mxcvec)/19)
  mx          <- rowSums(mxcvec)
  sd.frommx(mx,sex)
}

