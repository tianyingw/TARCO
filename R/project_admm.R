# ADMM solver for the max-norm PSD projection.
# Primal: min ||Q-S||_max, Q >= 0.
# Dual: max -<S,Y>, Y >= 0, sum(abs(Y)) <= 1.
# Acceptance uses feasible primal/dual bounds, not ADMM residuals alone.
project_admm <- function(S,eps=1e-6,max_iter=100000L,verbose=TRUE) {
  S<-(S+t(S))/2;p<-nrow(S);magnitude<-max(1,max(abs(S)));T<-S/magnitude
  psd<-function(M) {
    e<-eigen((M+t(M))/2,symmetric=TRUE)
    e$vectors%*%(pmax(e$values,0)*t(e$vectors))
  }
  ball<-function(V,radius) {
    a<-abs(as.numeric(V))
    if(sum(a)<=radius)return(V)
    u<-sort(a,decreasing=TRUE);cu<-cumsum(u)
    j<-max(which(u>(cu-radius)/seq_along(u)))
    sign(V)*pmax(abs(V)-(cu[j]-radius)/j,0)
  }
  X<-matrix(0,p,p);U<-X;rho<-1/p;alpha<-1.5
  bestQ<-diag(p)*1e-12;upper<-max(abs(bestQ-T));lower<-0;bestY<-X
  started<-Sys.time();trace<-list()
  for(k in seq_len(max_iter)) {
    Q<-psd(T+X-U)
    relaxed<-alpha*Q+(1-alpha)*(T+X)
    oldX<-X;V<-relaxed-T+U
    X<-V-ball(V,1/rho)
    U<-U+relaxed-T-X
    if(k%%25L==0L || k==1L) {
      feasibleQ<-Q+diag(p)*1e-12
      candidate<-max(abs(feasibleQ-T))
      if(candidate<upper){upper<-candidate;bestQ<-feasibleQ}
      Y<-psd(rho*U)+diag(p)*1e-12
      Y<-Y/max(1,sum(abs(Y)))
      bound<--sum(T*Y)
      if(bound>lower){lower<-bound;bestY<-Y}
      gap<-upper-lower
      if(k%%500L==0L || gap<=eps*(1+max(abs(upper),abs(lower)))) {
        if(verbose){cat("ADMM",k,"upper",upper,"lower",lower,"gap",gap,"rho",rho,"\n");flush.console()}
        trace[[length(trace)+1L]]<-c(iter=k,upper=upper,lower=lower,gap=gap,rho=rho)
      }
      if(gap>=-1e-10 && gap<=eps*(1+max(abs(upper),abs(lower)))) {
        # Independently check the stored best primal/dual pair before accepting.
        stopifnot(min(eigen(bestQ,symmetric=TRUE,only.values=TRUE)$values)>-1e-10,
          min(eigen(bestY,symmetric=TRUE,only.values=TRUE)$values)>-1e-10,
          sum(abs(bestY))<=1+1e-10,
          abs(max(abs(bestQ-T))-upper)<1e-10,
          abs(-sum(T*bestY)-lower)<1e-10)
        return(list(Q=magnitude*bestQ,status="optimal",adjustment=1e-12*magnitude,
          distance=magnitude*upper,solver="ADMM_primal_dual_certificate",
          certificate=list(upper=upper,lower=lower,gap=gap,scale=magnitude,
            dual=bestY,iterations=k,elapsed=Sys.time()-started,trace=trace)))
      }
      r<-sqrt(sum((Q-T-X)^2));d<-rho*sqrt(sum((X-oldX)^2))
      newrho<-if(r>10*d)min(rho*2,1e4) else if(d>10*r)max(rho/2,1e-8) else rho
      U<-U*(rho/newrho);rho<-newrho
    }
  }
  stop("ADMM did not attain the required primal-dual gap: ",upper-lower)
}
