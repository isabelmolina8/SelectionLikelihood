######################################################################
# Simulation experiment: Comparison of Molina-Ghosh method with 
# Pfeffermann-Sverchkov method (P-S 1999, Shankya) for estimation 
# of regression coefficients in birthweight (grams) in terms of 
# gestational age (weeks) with data from the U.S. NMIHS used by 
# Korn and Graubard (1995).
######################################################################

# Clean workspace and set directory

rm(list=ls())
setwd("C:/Isa/TrabajoMalay/Simulations_NMIHS1988_2020-02-13")

# Load required libraries

#library(maxLik)
source("lmInformative_FixedPoint_2020-02-13.R")

# Read new data set without missing values

datosf<-read.table("Data_Small.txt",header=T)
head(datosf)

dim(datosf)[1]
#[1] 9448

### Analysing the survey weights
### to see if the design is informative

# Original survey weight

sweight<-datosf$baseweight/100000

# Inclusion probabilities

pii<-1/sweight

# Boxplot of inclusion probabilities by race

boxplot(pii~datosf$race)
boxplot(1/sweight~datosf$race)

# Blacks (race=3) have the largest inclusion probabilities
# and whites (race=1) have the smallest ones

pii3<-as.vector(tapply(pii,datosf$race,mean))
pii3
#[1] 0.005288654 0.008450741 0.013232471
table(datosf$race)
#1    2    3 
#4250  295 4903 

# We select only whites

datosw<-datosf[datosf$race==1,]
dim(datosw)

# Take a SRS from the data to act as population
# We would like to have a small population

set.seed(2234)
N<-200
pop<-sample(1:dim(datosw)[1],size=N,replace=FALSE)
datosws<-datosw[pop,]

attach(datosws)
N<-dim(datosws)[1]
N

# Population mean
mean.t<-mean(bwg);mean.t
#[1] 3057.445

# Constants defining the expected inclusion probabilities

Y1<-1500
Y2<-2500
Yvec<-c(Y1,Y2)

# Group of birth weight

group<-numeric(N)
group[bwg<Y1]<-1 # Very low birth weight
group[(Y1<=bwg)&(bwg<Y2)]<-2 # Low birth weight
group[Y2<=bwg]<-3 # Normal birth weight
table(group)
datosw.n<-cbind(datosws,group)

sweight<-baseweight/100000
pii<-1/sweight
boxplot(pii~group)
boxplot(1/sweight~group)
# Very low and low birth weight infants have larger inclusion probabilities
# The sampling is very informative
piibyg<-tapply(pii,group,mean);piibyg
piibyg[1]/piibyg[3]
piibyg[2]/piibyg[3]

### Fit the regression model to the data (acting as population)
### to find the population values of the model parameters

# We plot the birthweight against gestage 

plot(gestage,bwg)
x11()
plot(gestage,log(bwg))

# We can see a quadratic or cubic relationship

# Regression fit by OLS to the original data (acting as population)

fit<-lm(bwg~gestage)

Xpop<-model.matrix(fit)
dim(Xpop)
p<-dim(Xpop)[2]

# Values of model parameters for the population

beta.t<-coef(fit)
beta.t

sigmae2.t<-sum(((bwg-fitted(fit))^2)/(N-p))
sigmae2.t
sigma(fit)^2

mean(Xpop%*%beta.t)
#[1] 3057.445
# Identical to true mean of bwg

D<-length(unique(group));D
Nd<-as.vector(table(group));Nd

# Definition of constants c1, c2 and c3

c1<-1/5
c2<-1/10
c3<-1/50
cvec<-c(c1,c2,c3)

nd<-round(cvec*Nd)
n<-sum(nd);n;nd

# Sort the data by group
datos.o<-rbind(datosw.n[group==1,],datosw.n[group==2,],datosw.n[group==3,])
attach(datos.o)

# Code with R functions

source("lmInformative_FixedPoint_2020-02-13.R")

# Set the seed for random number generation

#set.seed(3234)
set.seed(1234)

### Simulations start

M<-1000
beta.ols<-matrix(0,nr=M,nc=p)
beta.w<-matrix(0,nr=M,nc=p)
beta.ps<-matrix(0,nr=M,nc=p)
beta.mg<-matrix(0,nr=M,nc=p)
sm<-ht<-reg.w<-reg.ols<-reg.ps<-reg.mg<-rep(0,M)

for (sim in 1:M){ # Simulations start

#sim<-1
cat("sim=",sim,"\n\n")
  
# Draw a sample by systematic sampling within each stratum (group)

kvec<-1/cvec

samp<-NULL
sumNd<-0

# Draw a sample by SRS within each stratum

for (d in 1:D){
    samp<-c(samp,sample((sumNd+1):(sumNd+Nd[d]),size=nd[d],replace=FALSE))
    sumNd<-sumNd+Nd[d]
}

# Draw a sample by systematic sampling within each stratum

#for (d in 1:D){

#    kd<-kvec[d]
#    rd<-sample(1:kvec[d],size=1)+sumNd
    
#    samp1<-NULL
#    samp1[1]<-rd
    
#    ld<-1+floor((Nd[d]-rd)/kd)
#    for (j in 1:(ld-1)){
#        samp1[j+1]<-rd+j*kd
#    }
#    samp<-c(samp,samp1)
#    sumNd<-sumNd+Nd[d]
#}

datos.s<-datos.o[samp,]
bwg.s<-datos.s$bwg
gestage.s<-datos.s$gestage
group.s<-datos.s$group

# Sample mean

sm[sim]<-mean(bwg.s)

# HT estimator

ht.sum<-0
for (d in 1:D){
    ht.sum<-ht.sum+kvec[d]*sum(bwg.s[group.s==d])
}
ht[sim]<-ht.sum/N

# OLS fit 

fit.s<-lm(bwg.s~gestage.s)
Xs<-model.matrix(fit.s)
beta0<-coef(fit.s)
beta.ols[sim,]<-beta0
cat("OLS fit",beta.ols[sim,],"\n")

# Regression estimator based on OLS fit

reg.ols[sim]<-mean(Xpop%*%beta0)

# We set the value of sigmae for P-S method to the true value

#sigmae2<-sum(((bwg.s-fitted(fit.s))^2)/(n.s[sim]-p))
#sigmae<-sqrt(sigmae2)

sigmae<-sqrt(sigmae2.t)
sigmae2<-sigmae2.t

# WLS: pseudo-likelihood

Xst<-t(Xs)
ss<-as.vector(table(group.s))
W<-diag(rep(kvec,ss))
beta.w[sim,]<-solve(Xst%*%W%*%Xs)%*%Xst%*%W%*%matrix(bwg.s,nr=sum(ss),nc=1)

# Regression estimator based on OLS fit

reg.w[sim]<-mean(Xpop%*%beta.w[sim,])

# P-S method (parametric)

y<-bwg.s

beta.in<-beta0
opt.ps<-logLikFunPS_FP(beta.start=beta.in,y,Xs,sigmae,Yvec,cvec)
beta.ps[sim,]<-opt.ps[,1]
cat("P-S fit",beta.ps[sim,],"\n")

# Regression estimator based on PS fit

reg.ps[sim]<-mean(Xpop%*%beta.ps[sim,])

# M-G method

#theta<-(-0.2)
theta<-(-1)
#theta<-0.8

opt.mg<-logLikFunMG_FP(beta.start=beta.in,theta,y,Xs,sigmae,Yvec,cvec)
beta.mg[sim,]<-opt.mg[,1]
cat("M-G fit",beta.mg[sim,],"\n")

# Regression estimator based on MG fit

reg.mg[sim]<-mean(Xpop%*%beta.mg[sim,])

} # End of cicle for simulations


Abeta.ols<-colMeans(beta.ols)
Abeta.w<-colMeans(beta.w)
Abeta.ps<-colMeans(beta.ps)
Abeta.mg<-colMeans(beta.mg)
Abeta.ols
Abeta.w
Abeta.ps
Abeta.mg
beta.t

SDbeta.ols<-numeric(p)
SDbeta.w<-numeric(p)
SDbeta.ps<-numeric(p)
SDbeta.mg<-numeric(p)
for (k in 1:p){
  SDbeta.ols[k]<-sd(beta.ols[,k])
  SDbeta.w[k]<-sd(beta.w[,k])
  SDbeta.ps[k]<-sd(beta.ps[,k])
  SDbeta.mg[k]<-sd(beta.mg[,k])
}

sumres<-data.frame(beta.t,Abeta.ols,Abeta.w,Abeta.ps,Abeta.mg,SDbeta.ols,SDbeta.w,SDbeta.ps,SDbeta.mg)
round(sumres,3)


mean.t
sm.m<-mean(sm)
ht.m<-mean(ht)
reg.ols.m<-mean(reg.ols)
reg.w.m<-mean(reg.w)
reg.ps.m<-mean(reg.ps)
reg.mg.m<-mean(reg.mg)
sm.m;ht.m;reg.ols.m;reg.w.m;reg.ps.m;reg.mg.m

reg.ols.m-mean.t
reg.w.m-mean.t
reg.ps.m-mean.t
reg.mg.m-mean.t

sm.sd<-sd(sm)
ht.sd<-sd(ht)
reg.ols.sd<-sd(reg.ols)
reg.w.sd<-sd(reg.w)
reg.ps.sd<-sd(reg.ps)
reg.mg.sd<-sd(reg.mg)
sm.sd;ht.sd;reg.ols.sd;reg.w.sd;reg.ps.sd;reg.mg.sd


resmean<-data.frame(mean.t,sm.m,reg.ols.m,reg.w.m,reg.ps.m,reg.mg.m,sm.sd,reg.ols.sd,reg.w.sd,reg.ps.sd,reg.mg.sd)
resmean

100*abs((sumres[,2:5]-sumres[,1])/sumres[,1])

#write.table(sumres,"SummaryResults_M1000_SYS_n10_tm1.txt")
#write.table(resmean,"SummaryResultsMean_M1000_SYS_n10_tm1.txt")

write.table(sumres,"SummaryResults_M1000_SRS_n10_tm1.txt")
write.table(resmean,"SummaryResultsMean_M1000_SRS_n10_tm1.txt")

