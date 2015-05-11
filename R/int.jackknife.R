# Jackknife estimation

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

int.jackknife <- function(x, fun) {  #for multiple, object to "leave one out" in row (populations in columns)
  
  clase <- class(x)
  if(clase == "matrix" | clase == "data.frame") {
    method <- "multiple" 
    N <- nrow(x)
  } else if(class(x) == "vector" | class(x) == "integer" | class(x) == "numeric") {
    method <- "single"
    N <- length(x)
  } else {
    stop("x is not of class matrix, data.frame or vector")
  }
  
  crit <- abs(qt(0.975, N  -1))
  
  if(method == "sigle") {
    obs <- fun(x)
    jack <- vector()
    for(i in 1:N) {
      jack[i] <- fun(x[-i])
    }
    #mean of jackknife values
    pseudo <- (N * obs) - ((N - 1) * jack)                   #pseudo-value
    theta <- mean(pseudo)
    bias <- (N - 1) * (obs - theta)
    sd.jack <- sqrt(var(pseudo) / N)
    
    interval <- sd.jack * crit
    CI <- c(theta - interval, theta + interval)
    names(CI) <- c("lwr", "uppr")
    
  } else if(method == "multiple") {
    obs <- apply(x, 2, fun)
    obs2 <- matrix(obs, nrow= N,ncol=length(obs), byrow=T)
    jack <- x - x
    for(i in 1:N) {
      jack[i, ] <- apply(x[-i, ], 2, fun)
    }
    
    pseudo <- (N * obs2) - ((N - 1) * jack)
    theta <- apply(pseudo, 2, mean)
    bias <- (N - 1) * (obs - theta)
    pseudo.variance <- apply(pseudo, 2, var)		
    sd.jack <- sqrt(pseudo.variance / N)
    interval <- sd.jack * crit
    CI<- rbind(theta - interval, theta + interval)
    colnames(CI) <- colnames(x)
    rownames(CI) <- c("lwr", "uppr")
  }
  
  result <- list(obs = obs, 
                 theta = theta,
                 sd = sd.jack,
                 bias = bias,
                 CI = CI)
  
  result
}
