# H: Hurst index
H <- 0.25
TH <- 2*H

# n: Level of truncation
n <- 11
dt <- 2^{-(n+1)}

# Setseed
setSeed <- 112233

# Given (m, k), t^1_{m, k} = 2k/2^{m+1}, t^2_{m, k} = (2k+1)/2^{m+1}, t^3_{m, k} = (2k+2)/2^{m+1}
# Until n-th level, there will be 2^(n+1)-1 many of pairs (m, k) for m = 0, ..., n, k = 0, ..., 2^m-1
# The pair (m, k) corresponds to l = 2^m+k entry of a vector (or matrix)

# Make the covariance matrix sigma
sigma <- matrix(0.0, nrow = 2^(n+1)-1, ncol = 2^(n+1)-1)

for(m in 0:n){
  for(k in 0:(2^m-1)) {
    for(mp in 0:n){
      for(kp in 0:(2^mp-1)) {
        if(m == mp & k == kp) {
          # l is a column number corresponding to (m, k)
          l <- 2^m+k
          sigma[l, l] <- (abs(2-2^(2*H-1))*2^((m+1)*(1-2*H)))
        }
        else {
          l <- 2^m+k
          lp <- 2^mp+kp
          
          t1mk <- (2*k)/2^{m+1}
          t2mk <- (2*k+1)/2^{m+1}
          t3mk <- (2*k+2)/2^{m+1}
          t1mpkp <- (2*kp)/2^{mp+1}
          t2mpkp <- (2*kp+1)/2^{mp+1}
          t3mpkp <- (2*kp+2)/2^{mp+1}
          
          xi11 <- 0.5*(abs(t1mk-t2mpkp)^{TH}+abs(t2mk-t1mpkp)^{TH}-abs(t1mk-t1mpkp)^{TH}-abs(t2mk-t2mpkp)^{TH})
          xi22 <- 0.5*(abs(t2mk-t3mpkp)^{TH}+abs(t3mk-t2mpkp)^{TH}-abs(t2mk-t2mpkp)^{TH}-abs(t3mk-t3mpkp)^{TH})
          xi21 <- 0.5*(abs(t2mk-t2mpkp)^{TH}+abs(t3mk-t1mpkp)^{TH}-abs(t2mk-t1mpkp)^{TH}-abs(t3mk-t2mpkp)^{TH})
          xi12 <- 0.5*(abs(t1mk-t3mpkp)^{TH}+abs(t2mk-t2mpkp)^{TH}-abs(t1mk-t2mpkp)^{TH}-abs(t2mk-t3mpkp)^{TH})
          
          sigma[l, lp] <- ((xi11+xi22-xi21-xi12)*(2^(0.5*(m+mp))))
        }
      }
    }
  }
}

# correlation matrix corresponding to sigma
cor_sigma <- cov2cor(sigma)


### Generate uniform random variables using Copula with the given correlation matrix cor_sigma
library(copula)

#d: dimension
d <- (2^(n+1)-1)
mcop <- normalCopula(P2p(cor_sigma), dim = d, dispstr="un")
cop1 <- mvdc(copula = mcop, margins = c("unif"), paramMargins = list(list(min = -sqrt(3), max = sqrt(3))), marginsIdentical = TRUE)
uniforms <- rMvdc(1, cop1)

#Scale uniforms(-sqrt(3), sqrt(3)) with relevant variances
theta <- rep(0, d)
for(i in 1:d) {
  theta[i] <- uniforms[i]*sqrt(sigma[i, i])
}


# Compute the function values over 2^(n+1)+1 discrete time points of [0, 1]
data <- data.frame(time = seq(0,1,by=dt))
fct = rep(0,2^(n+1)+1)

for(m in 0:n){
  slope <- 2^{0.5*(m)}
  for(k in 0:(2^m-1)) {
    l <- 2^m+k
    
    t1 <- (2*k)/2^{m+1}
    t2 <- (2*k+1)/2^{m+1}
    t3 <- (2*k+2)/2^{m+1}
    
    t1_index <- (t1 * (2^{n+1}) + 1)
    t2_index <- (t2 * (2^{n+1}) + 1)  # Add 1 because time includes the origin as the first element
    t3_index <- (t3 * (2^{n+1}) + 1)
    
    for(i in t1_index:t2_index) {
      fct[i] <- fct[i] + theta[l] * (i-t1_index)*dt*slope
    }
    for(i in (t2_index+1):t3_index) {
      fct[i] <- fct[i] + theta[l] * (t3_index-i)*dt*slope
    }
  }
}

# Store function values in a dataframe
data$fct <- fct
colnames(data) <- c("time", "fct")

# 1/H-variation of data$fct
p <- (1/H)
data$pvar <- c(0, cumsum(abs(diff(data$fct))^p))

# Scaled quadratic variation of data$fct
data$scaled_qv <- c(0, (1:(2^(n+1)))^(2*H-1)*cumsum(abs(diff(data$fct))^2))