######################################################################
##### Pfeffermann and Sverchkov's method
######################################################################

#-----------------------------------------------------------------
# Maximising with respect to beta only
#-----------------------------------------------------------------

## Loglikelihood

logLikFunPS_FP <- function(beta.start,y,Xs,sigmae,Yvec,cvec) {
  beta <- beta.start
  p<-length(beta)
  n <- length(y)
  sigmae2<-sigmae^2
  Y1<-Yvec[1]
  Y2<-Yvec[2]
  c1<-cvec[1]
  c2<-cvec[2]
  c3<-cvec[3]
  n1<-sum(y<Y1)
  n2<-sum((y>=Y1)&(y<Y2))
  n3<-sum(y>=Y2)
  Xst<-t(Xs)

  M<-1
  iterMax<-1000
  iter<-0
  
  while((M>0.0001)&(iter<iterMax)){
    
  iter<-iter+1
  #cat("iter=",iter,"\n")

  sumbetai<-numeric(p)
  for (i in 1:n){
    Xsibeta<-Xs[i,]%*%beta
    p1<-pnorm((Y1-Xsibeta)/sigmae)
    p2<-pnorm((Y2-Xsibeta)/sigmae)
    di<-c1*p1+c2*(p2-p1)+c3*(1-p2)

    d1<-dnorm(as.numeric((Y1-Xsibeta)/sigmae))
    d2<-dnorm(as.numeric((Y2-Xsibeta)/sigmae))
    Ddibeta<- ( (-c1)*d1-c2*(d2-d1)+c3*d2)*Xs[i,]/sigmae
    
    sumbetai<-sumbetai+Ddibeta/di[1,1]
  }
  
  betak<-solve(Xst%*%Xs)%*%(Xst%*%matrix(y,nr=n,nc=1)-sigmae2*matrix(sumbetai,nr=p,nc=1))
  M<-max(abs((betak-beta)/beta))
  beta<-betak
  #print(betak)
  }
  
  return(betak)
}



######################################################################
##### Molina and Ghosh's method
######################################################################

# Optimizing only beta, theta known

library(numDeriv)

logLikFunMG_FP <- function(beta.start,theta,y,Xs,sigmae,Yvec,cvec) {
#beta.start=beta.in;theta=theta;y=y;Xs=Xs;sigmae=sigmae;Yvec=Yvec;cvec=cvec
  beta <- beta.start
  p<-length(beta)
  n <- length(y)
  sigmae2<-sigmae^2
  Y1<-Yvec[1]
  Y2<-Yvec[2]
  c1<-cvec[1]
  c2<-cvec[2]
  c3<-cvec[3]
  n1<-sum(y<Y1)
  n2<-sum((y>=Y1)&(y<Y2))
  n3<-sum(y>=Y2)
  Xst<-t(Xs)
 
  #cat("MG method","\n")
  
  M<-1
  iterMax<-1000
  iter<-0
  
  while((M>0.0001)&(iter<iterMax)){
    
    iter<-iter+1
    #cat("iter=",iter,"\n")
    
    sumbetai<-numeric(p)
    #proddfi<-1
    for (i in 1:n){
      
      Xsibeta<-0
      for (k in 1:p){
        Xsibeta<-Xsibeta+Xs[i,k]*beta[k]
      }
      
      p1<-pnorm((Y1-Xsibeta)/sigmae)
      p2<-pnorm((Y2-Xsibeta)/sigmae)
      di<-c1*p1+c2*(p2-p1)+c3*(1-p2)
      
      d1<-dnorm((Y1-Xsibeta)/sigmae)
      d2<-dnorm((Y2-Xsibeta)/sigmae)
      Ddibeta<-( (-c1)*d1-c2*(d2-d1)+c3*d2)*Xs[i,]/sigmae

      sumbetai<-sumbetai+Ddibeta/di
      
      #dfi<-1-( (c1^2)*p1+(c2^2)*(p2-p1)+(c3^2)*(1-p2) )/di     
      #proddfi<-proddfi*(dfi)
    }
    
    # Now we do the numerical derivatives
    A<-function(x){
        proddfi.s<-1
        for (i in 1:n){
      
        Xsibeta.s<-0
        for (k in 1:p){
            Xsibeta.s<-Xsibeta.s+Xs[i,k]*x[k]
        }
      
        p1.s<-pnorm((Y1-Xsibeta.s)/sigmae)
        p2.s<-pnorm((Y2-Xsibeta.s)/sigmae)
        di.s<-c1*p1.s+c2*(p2.s-p1.s)+c3*(1-p2.s)

        dfi.s<-1-( (c1^2)*p1.s+(c2^2)*(p2.s-p1.s)+(c3^2)*(1-p2.s) )/di.s
      
        proddfi.s<-proddfi.s*dfi.s
        #cat("dfi=",dfi.s,"proddfi=",proddfi.s,"\n")
        }
        return(proddfi.s)
    }
    
    DAbeta<-grad(A,x=beta)
    Abeta<-A(beta)

    #betak<-solve(Xst%*%Xs)%*%(Xst%*%matrix(y,nr=n,nc=1)-sigmae2*matrix(sumbetai,nr=p,nc=1)+sigmae2*theta*matrix(DAbeta,nr=p,nc=1)/(1+theta*Abeta))
    betak<-solve(Xst%*%Xs)%*%(Xst%*%matrix(y,nr=n,nc=1)-sigmae2*matrix(sumbetai,nr=p,nc=1)-sigmae2*theta*matrix(DAbeta,nr=p,nc=1)/(1+theta*Abeta))
    M<-max(abs((betak-beta)/beta))
    beta<-betak
    #print(betak)
  }# End of while
  
  return(betak)
}

