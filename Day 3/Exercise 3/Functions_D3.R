# Decomposition class - EDSD 2022
# FUNCTIONS AND INFORMATION FOR DAY 3

# Information for graphs
# Labesls for causes of death
cause_names<-c("1"="Infectious and respiratory", "2"="Neoplasm","3"="Circulatory and heart",
               "4"="Birth","5"="Diabetes",
               "6"= "Homicide","7"="Other external",
               "8"="Other","9"="Undefined")


# Labels for age-groups 
age_names<-c("0"="0","1"="1-4","5"="5-9","10"="10-14","15"="15-19","20"="20-24",
             "25"="25-29","30"="30-34","35"="35-39","40"="40-44","45"="45-49",
             "50"="50-54","55"="55-59","60"="60-64","65"="65-69","70"="70-74",
             "75"="75-79","80"="80-84","85"="85+")     


# Functions to use with DemoDecomp functions
# Life expectancy at birth
e0.frommx <- function(nmx =  mx, sex=1, nax = NULL){
  age <- 0:(length(nmx)-1)
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

# Lifespan disparity
edagger.frommx <- function(mx,sex=1){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  if (mx[i.openage] < 0.5 | is.na(mx[i.openage])) mx[i.openage] = mx[i.openage - 1]*1.1 
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  ax[i.openage] <- if (mx[i.openage] == 0) 0.5 else 1/mx[i.openage]
  
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  
  
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  #dx[i.openage] <-0
  Lx 				    <- lx - (1 - ax) * dx
  Lx[i.openage ]	    <- dx[i.openage ] * ax[i.openage ]
  Lx[is.na(Lx)] <- 0
  Tx 				    <- c(rev(cumsum(rev(Lx[1:OPENAGE]))),0) + Lx[i.openage]
  ex 				    <- Tx / lx
  ex[is.na(ex)] <- 0
  ex[i.openage] <- if (ex[OPENAGE] == 0) 0 else ax[i.openage]
  
  v        <- (ax*c(ex[-1L],0) + (1-ax)*ex)
  v[length(ex)] <- ex[length(ex)]
  v <- dx*v
  e.dagger <- rev(cumsum(rev(v)))/lx
  e.dagger[16]
}

AKm02a0 <- function(m0, sex = "m"){
  sex <- rep(sex, length(m0))
  ifelse(sex == "m", 
         ifelse(m0 < .0230, {0.14929 - 1.99545 * m0},
                ifelse(m0 < 0.08307, {0.02832 + 3.26201 * m0},.29915)),
         # f
         ifelse(m0 < 0.01724, {0.14903 - 2.05527 * m0},
                ifelse(m0 < 0.06891, {0.04667 + 3.88089 * m0}, 0.31411))
  )
}


edaggerfrommxc <- function(mxcvec,sex=1){
  dim(mxcvec) <- c(110,length(mxcvec)/110)
  mx          <- rowSums(mxcvec)
  edagger.frommx(mx,sex)
}


function (func, pars1, pars2, N, ...) 
{
  y1 <- e0.frommx(pars1)
  y2 <- func(pars2, ...)
  d <- pars2 - pars1
  n <- length(pars1)
  delta <- d/N
  x <- pars1 + d * matrix(rep(0.5:(N - 0.5)/N, n), byrow = TRUE, 
                          ncol = N)
  cc <- matrix(0, nrow = n, ncol = N)
  zeros <- rep(0, n)
  for (j in 1:N) {
    DD <- diag(delta/2)
    for (i in 1:n) {
      cc[i, j] <- func((x[, j] + DD[, i]), ...) - func((x[, 
                                                          j] - DD[, i]), ...)
    }
  }
  return(rowSums(cc))
}
