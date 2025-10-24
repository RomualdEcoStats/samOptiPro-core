

# nlchoose - Binomial coefficient calculator ---------------
nlchoose <- nimbleFunction(
  run = function(n = double(0), k = double(0)){
    returnType(double(0))
    k <- round(k)
    n <- round(n)
    if(n==k | k==0) return(0)
    else if(n < k) return(-Inf)
    else return(sum(log((n-k+1):n)) - sum(log(1:k)))
  }
) 

# dhyperg - hypergeometric density function ---------------
dhyperg <- nimbleFunction(
  run = function(x = double(0), nx = double(0), N1 = double(0), N0 = double(0),
                 log = integer(0, default = 1)){
    returnType(double(0))
    logProb <- nlchoose(N1,x) + nlchoose(N0,nx-x) - nlchoose(N0+N1,nx)
    if(log) return(logProb)
    else return(exp(logProb))
  })

# rhyperg - hypergeometric RNG -------------
rhyperg <- nimbleFunction(
  run = function(n = integer(0, default = 1), nx = double(0), N1 = double(0), N0 = double(0)){
    returnType(double(0))
    if(n!=1) print('rhyperg only allows n=1; using n = 1')
    dev <- runif(1,0,1)
    minx <- max(nx - N0, 0)
    maxx <- min(nx, N1)
    xs <- numeric(length = maxx-minx+1)
    ps <- xs
    xs[1] <- minx
    ps[1] <- dhyperg(minx, nx, N1, N0, log = 0)
    
    for(i in (minx+1):maxx){
      xs[i-minx+1] <- i
      ps[i-minx+1] <- dhyperg(i, nx, N1, N0, log = 0) + ps[i-minx]
    }
    
    # return(max(xs[dev > ps]))
    return(xs[which(! dev > ps)[1]])
  }
)

# register Nimble functions to R's global environment ----------------
assign('nlchoose', nlchoose, .GlobalEnv)
assign('dhyperg', dhyperg, .GlobalEnv)
assign('rhyperg', rhyperg, .GlobalEnv)
