#Bonus (3)
#Pricing European Up-and-Out Put Option
Tri.Am.Put <- function(K,T,S,r,H,X_rebate,N){
  dt = T/N
  nu = r-0.5*sig^2
  dx = sqrt(sig^2*dt + (nu*dt)^2)
  pu = 0.5 * ( (sig^2*dt + nu^2 *dt^2)/dx^2 + nu*dt/dx )
  pm = 1.0 -   (sig^2*dt + nu^2 *dt^2)/dx^2 
  pd = 0.5 * ( (sig^2*dt + nu^2 *dt^2)/dx^2 - nu*dt/dx )
  #precompute constants
  disc = exp(-r*dt)
  
  #initialise asset prices at maturity N
  St = array(0,dim = c(N+1))
  St[1] <- S * exp(-dx*N)
  for(j in 2:N+1){
    St[j] = St[j-1]*exp(2*(dx))
  }
  
  #initialise option values at maturity
  C = array(0,dim = c(N+1,2*N+1))
  for(j in 1:(2*N+1)){
    
    if(St[j]<H)
      C[N+1,j] <- max(0.0,K-St[j])
    else
      C[N+1,j] <- 0.0
  }
  
  #step back through the tree applying the barrier and early exercise condition
  for(i in N:1){
    for(j in 1:(2*N+1)){
      St[j] = St[j]*exp(dx)
      if(St[j]<H){
        C[i,j] = disc*(pu*C[i+1,j+1] + pm*C[i+1,j]+pd*C[i+1,j-1])
        
        #apply the early exercise condition
        C[i,j]=max(C[i,j],K-St[j])
      }
      else
        C[i,j] = X_rebate
    }
  }
  C[1,1]  
}