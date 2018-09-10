#############################################################################                                                                       ##
## Accurate error control in high dimensional association testing using  ##
##  conditional false discovery rates                                    ##
##                                                                       ##
## Analyse results of simulations, draw plots, and analyse GWAS data     ##
##                                                                       ##
## James Liley and Chris Wallace                                         ##
##                                                                       ##
###########################################################################
#
# This R script will deterministically generate all plots and outputs used
#  in the paper above. Outputs are generated in the same order that they 
#  appear in the paper, with supplementary material generated according to
#  the first referenced in the main  text. 
#
# This script requires all requisite packages to be installed and the 
#  working directory to be set to the file one level up from where this 
#  script is saved (it should contain subdirectories ./code, ./data, 
#  ./simulations and ./plots). All outputs are written to ./outputs.
#
# Simulations are generated using the R script ./code/run_simulation.R. 
#  Each run of the simulation is dependent on a random seed. A matrix of 
#  completed simulation outputs with associated random seeds is stored in 
#  ./data/cfdrsimmatrix.RData.
#
# Input data from the GWAS analysis is included in the R object 
#  ./data/gwas_data.RData. This file contains summary statistics for JIA
#  and T1D along with chromosome number and position (build 37) for all 
#  variants. Data were originally sourced from the ImmunoBase repository,
#  (www.immunobase.org) and have already undergone quality control and 
#  genomic control as per Liley et al 2015. The MHC region has already 
#  removed.
#
# The analysis of GWAS data runs quite slowly. It will only run if the 
#  file ./data/gwas_cfdr_data.RData does not already exist. A copy of
#  this file is already included in the ./data directory. To re-run the 
#  GWAS analysis, delete this file and run this script. The file 
#  generated will be identical to that included in the directory.
#  

###########################################################################
## Packages, directories and switches #####################################
###########################################################################

# Packages

library(cfdr)
library(mnormt)
library(mgcv)
library(pbivnorm)
library(MASS)
library(fields)
library(matrixStats)
library(plotrix)


# Directories

# This variable gives the location of the folder cfdr_control. It should 
#  contain a document README, and four folders  named ./code, ./data, 
#  ./simulations, and ./outputs. Here it is assumed to be the current working 
#  directory
cfdr_dir="./"

# This variable gives the location to save outputs (plots, tables) to.
output_dir=paste0(cfdr_dir,"outputs/")



# Switches

# set to F to draw plots rather than saving them as PDFs
save_pdf=T 




###########################################################################
## Read simulation data ###################################################
###########################################################################

# Load file from directory ./data/. See file ./data/cfdrsimmatrix_README 
#  for details.
load(paste0(cfdr_dir,"data/cfdrsimmatrix.RData"))



###########################################################################
## Figure showing failure of cFDR<alpha control at alpha ##################
###########################################################################

# This figure does not require any data

if (save_pdf) pdf(paste0(output_dir,"ml.pdf"),width=5,height=5)

set.seed(1)

plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="P",ylab="Q",xaxs="i",yaxs="i")

q=seq(0,1,length.out=500); k=1/10; cf=0.5; p=cf*k/ (q - cf*q + k)
lines(p,q,col="red")

polygon(c(0,0,p,0,0),c(0,0,q,1,1),col="lightgray",border="red")
polygon(c(0,p[200],p[200],0),c(q[200],q[200],0,0),col="darkgray",border="black")

points(abs(rnorm(50))/50,abs(rnorm(50))/50,pch=16,cex=0.5)
points(p[200],q[200])

text(p[200]/2,q[200]/2,"M")
text(0.05,0.7,"L")


text(p[200]+0.03,q[200]+0.03,"p,q")

if (save_pdf) dev.off()




###########################################################################
## Behaviour of L-curves using FDR-control methods 1,4 ####################
###########################################################################

# L examples

# Simulation
set.seed(1) # reproducibility

nsnp=10000; n1p=100; n1q=100; n1pq=100; sp=2; sq=2

zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(nsnp-n1p-n1pq,sd=1))
zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(nsnp-n1p-n1pq-n1q,sd=1))

# P-values
p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))


# Subset of curves at which to plot curves/points
xmin=5e-4
sub=which(p< xmin)
highlight=c(3,10,16) # highlight these curves




###########################################################################
## Plot of points using FDR control method 2 ##############################
###########################################################################
 
if (save_pdf) pdf(paste0(output_dir,"Lplot_pcut.pdf"),width=6,height=5)

# Compute L-curves using method 2 (mode=0). This step is slow due to the high resolution
vx2=vl(p,q,indices=sub,adj=F,mode=0,nv=25000); 

# Plot title
ex=expression(atop(L[S](alpha[i]),alpha[i]==paste("min"[p >= p[i]],"{",{widehat(cFDR)} 
                       [paste(S+"(p,q"[i],")")](p[i],q[i]),"}")))

plot(0,0,type="n",xlab="P",ylab="Q",xlim=c(0,xmin),ylim=c(0,1),main=ex,
     xaxs="i",yaxs="i")

for (i in 1:length(sub)) lines(vx2$x[i,],vx2$y,col="gray")

for (i in highlight) {
  lines(vx2$x[i,],vx2$y,col="blue"); points(p[sub[i]],q[sub[i]],cex=2,col="blue")
}

points(p,q,col="red",cex=0.5)

if (save_pdf) dev.off()



###########################################################################
## Examination of 'points on stalks' effect (Supplementary figure) ########
###########################################################################

if (save_pdf) pdf(paste0(output_dir,"local_l.pdf"),width=6,height=5)

i=1
plot(vx2$x[i,],vx2$y,type="l",xlim=c(0,2*p[sub[i]]),ylim=c(0,0.05),xlab="P",ylab="Q"); 
points(p[sub[i]],q[sub[i]],col="blue",pch=16,cex=2)

k=200
points(vx2$x[i,5945-k],vx2$y[5945-k],col="black",pch=1,cex=2)
text(vx2$x[i,5945-k]*1.2,vx2$y[5945-k]*1.2,"A",cex=2)

points(vx2$x[i,5945+k],vx2$y[5945+k],col="black",pch=1,cex=2)
text(vx2$x[i,5945+k]*0.8,vx2$y[5945+k]*0.8,"B",cex=2)

if (save_pdf) dev.off()



###########################################################################
## Plot of points using FDR control method 4 ##############################
###########################################################################

if (save_pdf) pdf(paste0(output_dir,"Lplot_pcut_leaveout.pdf"),width=6,height=5)

vx4=vl(p,q,indices=sub,adj=F,mode=1); 

ex=expression(atop({L}[paste("S-(p"[i],",q"[i],")")](alpha[i])
                          ,alpha[i]==paste("min"[p >= p[i]],"{",{widehat(cFDR)} 
                          [paste(S+"(p,q"[i],")-(p"[i],",q"[i],")")] 
                          (p[i],q[i]),"}"))) 


plot(0,0,type="n",xlab="P",ylab="Q",xlim=c(0,xmin),ylim=c(0,1), main=ex)

for (i in 1:length(sub)) lines(vx4$x[i,],vx4$y,col="gray")

for (i in highlight) {
  lines(vx4$x[i,],vx4$y,col="blue"); points(p[sub[i]],q[sub[i]],cex=2,col="blue")
}


points(p,q,col="red",cex=0.5)

if (save_pdf) dev.off()





###########################################################################
## Numerical FDR comparison between methods ###############################
###########################################################################

# The output of this section of code is saved in ./outputs as a text file
#  ./outputs/FDR_comparison.txt.

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
attach(rx)

# Shorthands
fp=fdp_p; # FDP for p-values (BH procedure)
f1=fdp_cf1_fdr1_adj0_dist1; # FDP using method 1
f2=fdp_cf1_fdr2_adj0_dist1; # FDP using method 2
f3=fdp_cf1_fdr3_adj0_dist1; # FDP using method 3
f3b=fdp_cf1_fdr3b_adj0_dist1; # FDP using method 3b
f4=fdp_cf1_fdr4_adj0_dist1; # FDP using method 4

f3_1=fdp_cf1_fdr3_adj0_dist1_fold1; # FDP using method 3, fold 1
f3_2=fdp_cf1_fdr3_adj0_dist1_fold2; # FDP using method 3, fold 2
f3_3=fdp_cf1_fdr3_adj0_dist1_fold3; # FDP using method 3, fold 3

# We separately consider the cases n1p+ n1pq=0 and n1p+n1pq>0, since 
#  failure of FDR control tends to occur most severely at n1p+n1pq=0
xx=n1p+n1pq; 
x0=which(xx==0)
x1=which(xx>0)

alpha=0.1 # universal


# FDR control at n1p+n1pq=0. In this case, any rejection of the 
#  null constitutes a false discovery, so the FDP is either 0 or 1.
#  We thus use a one-tailed binomial test of proportion to see if 
#  the FDR exceeds alpha (0.1).
n=length(x0)
prop.test(sum(f1[x0]),n,alpha,alternative="greater") # controlled using method 1
prop.test(sum(f2[x0]),n,alpha,alternative="greater") # NOT controlled using method 2
prop.test(sum(f3[x0]),n,alpha,alternative="greater") # NOT controlled using method 3
prop.test(sum(f3b[x0]),n,alpha,alternative="greater") # controlled using method 3b
prop.test(sum(f4[x0]),n,alpha,alternative="greater") # controlled using method 4
prop.test(sum(c(f3_1[x0],f3_2[x0],f3_3[x0])),
          3*n,alpha,alternative="greater") # controlled within-fold using method 3a

# When used on f1,f2,f3,f3b,f4, this test has power 90% to detect a deviation of 9.5% at 0.05
#  significance:
POWER=0.9; P=0.05; n=length(x0); alpha=0.1
uniroot(function(x) qbinom(1-POWER,n,x)-qbinom(1-P,n,alpha),interval=c(0.05,0.5))$root / alpha # proportional change in alpha detectable



# FDR control at n1p+n1pq>0. FDP is roughly normal for large n1p+n1pq; we
#  use a one-sided t-test here. In general, the FDR is controlled at a 
#  level lower than alpha when n1p+n1pq>0 (see figures). 
t.test(f1[x1]-alpha,alternative="greater") # controlled using method 1
t.test(f2[x1]-alpha,alternative="greater") # NOT controlled using method 2
t.test(f3[x1]-alpha,alternative="greater") # NOT controlled using method 3
t.test(f3b[x1]-alpha,alternative="greater") # controlled using method 3b
t.test(f4[x1]-alpha,alternative="greater") # controlled using method 4
t.test(c(f3_1[x0],f3_2[x0],f3_3[x0])-alpha,
         alternative="greater") # controlled within-fold using method 3a
# When used on f1,f2,f3,f3b,f4, this test has power to detect a deviation of 11% at 0.05 significance:
POWER=0.9; P=0.05; n=length(x1); SD=sd(c(f1,f2,f3,f3b,f4)); alpha=0.1
(power.t.test(n=n,delta=NULL,sd=SD,sig=P,power=POWER,alt="one")$delta / alpha) + 1 # proportional change in alpha detectable



# FDR comparison to p-values with n1p+n1pq=0: the FDR of the B-H method at n1p+n1pq=0 is exactly alpha,
#  so we can use a two-sided binomial test of proportion. Method 3a cannot easily be compared in this
#  way within-fold as it uses a different sample size to the B-H method.
n=length(x0)
prop.test(sum(f1[x0]),n,alpha,alternative="two") # no difference using method 1
prop.test(sum(f2[x0]),n,alpha,alternative="two") # FDR higher using method 2
prop.test(sum(f3[x0]),n,alpha,alternative="two") # FDR higher using method 3
prop.test(sum(f3b[x0]),n,alpha,alternative="two") # no difference using method 3b
prop.test(sum(f4[x0]),n,alpha,alternative="two") # no difference using method 4
# When used on f1,f2,f3,f3b,f4, this test has power 90% to detect a deviation of 10.5% at 0.05
#  significance:
POWER=0.9; P=0.05; n=length(x0); alpha=0.1
uniroot(function(x) qbinom(1-POWER,n,x)-qbinom(1- P/2,n,alpha),interval=c(0.05,0.5))$root / alpha # proportional change in alpha detectable





# FDR comparison to p-values with n1p+n1pq>0. The FDR of the B-H method in
#  this case is generally <0.1, so we compare pairwise. We use a paired 
#  Wilcoxon rank sum test. The sign of the difference is given by the sign
#  of the pseudomedian; that is, if the pseudomedian of wilcox.test(x,y,...) 
#  is >0, then in general x>y.
wilcox.test(f1[x1],fp[x1],paired=T,conf.int=T) # FDP SMALLER using method 1 (ie overconservative)
wilcox.test(f2[x1],fp[x1],paired=T,conf.int=T) # FDP larger using method 2 
wilcox.test(f3[x1],fp[x1],paired=T,conf.int=T) # FDP larger using method 3
wilcox.test(f3b[x1],fp[x1],paired=T,conf.int=T) # FDP equivalent using method 3b
wilcox.test(f4[x1],fp[x1],paired=T,conf.int=T) # FDP equivalent using method 4
# It is difficult to sensibly calculate the power of a Wilcoxon test; however,
#  the test is powerful. Approximately 20% of FDPs are equal in each comparison. We show below
#  that we have approximately 90% power to detect a 3% difference in 
#  Pr(FDP(method) > FDP(P)) and Pr(FDP(method) < FDP(P))
delta=0.03 # P(FDP(method) > FDP(P)) - P(FDP(method) < FDP(P))
equal=0.2 # P(FDP(method)=FDP(P))
n=length(x1); np=0; ntrial=1000; P=0.05; set.seed(1)
for (i in 1:ntrial) {
 s=sample(c(-1,0,1),n,prob=c((1-equal-delta)/2,equal,(1-equal+delta)/2),rep=T)
 pw=wilcox.test(s,rep(0,n),paired=T)$p.value; if (pw< P) np=np+1
}
np/ntrial # power to detect difference





###########################################################################
## Numerical TDR comparison between methods ###############################
###########################################################################

# The output of this section of code is saved in ./outputs as a text file
#  ./outputs/TDR_comparison.txt.

# Shorthands. We do not consider FDR controlling methods 2 or 3, since 
#  these failed to control FDR. We also do not consider the power of method
#  3 within-fold, as this will effectively be captured by considering the 
#  power of method 3b.
tp=tdr_p; # TDR for p-values (BH procedure)
t1=tdr_cf1_fdr1_adj0_dist1; # TDR using method 1
t3b=tdr_cf1_fdr3b_adj0_dist1; # TDR using method 3b
t4=tdr_cf1_fdr4_adj0_dist1; # TDR using method 4

# We only consider values with n1p+n1pq>0, since TDR is indeterminate if 
#  n1p+n1pq==0
x1=which(xx>0)


# TDR comparison between methods, using paired Wilcoxon rank sum tests.
#  The sign of the difference is given by the sign of the pseudomedian;
#  that is, if the pseudomedian of wilcox.test(x,y,...) is >0, then 
#  in general x>y.
wilcox.test(t1[x1],tp[x1],paired=T,conf.int=T) # Method 1 is more powerful than p-value
wilcox.test(t3b[x1],tp[x1],paired=T,conf.int=T) # Method 3b is more powerful than p-value
wilcox.test(t4[x1],tp[x1],paired=T,conf.int=T) # Method 3b is more powerful than p-value

wilcox.test(t4[x1],t1[x1],paired=T,conf.int=T) # Method 4 is more powerful than method 1
wilcox.test(t3b[x1],t1[x1],paired=T,conf.int=T) # Method 3b is more powerful than method 1

wilcox.test(t4[x1],t3b[x1],paired=T,conf.int=T) # Methods 3b and 4 are similarly powerful.



# It is difficult to sensibly calculate the power of a Wilcoxon test; however,
#  the test is powerful. Approximately 63% of FDPs are equal in each comparison. We show below
#  that we have approximately 90% power to detect a 2% difference in 
#  P(FDP(method A) > FDP(method B)) and P(FDP(method A) < FDP(method B))
delta=0.02 # P(FDP(method A) > FDP(method B)) - P(FDP(method A) < FDP(method B))
equal=0.63 # P(FDP(method A)=FDP(method B))
n=length(x1); np=0; ntrial=1000; P=0.05; set.seed(1)
for (i in 1:ntrial) {
 s=sample(c(-1,0,1),n,prob=c((1-equal-delta)/2,equal,(1-equal+delta)/2),rep=T)
 pw=wilcox.test(s,rep(0,n),paired=T)$p.value; if (pw< P) np=np+1
}
np/ntrial # power to detect difference



detach(rx)


###########################################################################
## Plot of FDR-controlling methods 1-4 FDR comparison #####################
###########################################################################

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
attach(rx)

# shorthand for relevant variables
f1=fdp_cf1_fdr1_adj0_dist1
f2=fdp_cf1_fdr2_adj0_dist1
f3=fdp_cf1_fdr3_adj0_dist1
f3b=fdp_cf1_fdr3b_adj0_dist1
f4=fdp_cf1_fdr4_adj0_dist1
fp=fdp_p
xx=n1p+n1pq

sdx=0.15; # standard deviation for smoothing
sdx=sdx*(max(xx)-min(xx)); # scale standard deviation

w=which(pmin(f1,f2,f3,f3b,f4,fp,xx) > -0.5) # include only complete observations
f1=f1[w]; f2=f2[w]; f3=f3[w]; f3b=f3b[w]; f4=f4[w]; fp=fp[w]; xx=xx[w]

u=sort(unique(xx)); 
f1x=0*u; f2x=0*u; f3x=0*u; f3bx=0*u; f4x=0*u; fpx=0*u; 
for (i in 1:length(u)) {
    wts=dnorm(xx-u[i],sd=sdx); wts=wts/sum(wts)
    f1x[i]=sum(f1*wts); f2x[i]=sum(f2*wts); f3x[i]=sum(f3*wts);
    f3bx[i]=sum(f3b*wts); f4x[i]=sum(f4*wts); fpx[i]=sum(fp*wts)
}

if (save_pdf) pdf(paste0(output_dir,"fdp_m1234_noci.pdf"),width=6,height=6)

plot(0,type="n",xlim=range(xx),ylim=c(0.05,0.25),xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab="FDR",bty="n")

su=sort(u); ou=order(u)
lines(su,fpx[ou],col="black",lwd=2,lty=1); 
lines(su,f1x[ou],col="red",lwd=2,lty=2); 
lines(su,f2x[ou],col="blue",lwd=2,lty=3); 
lines(su,f3x[ou],col="darkgray",lwd=2,lty=4); 
lines(su,f3bx[ou],col="darkgray",lwd=2,lty=5); 
lines(su,f4x[ou],col="purple",lwd=2,lty=6); 

legend(min(xx) + 0.7*(max(xx)-min(xx)),0.23,c("BH","1","2","3a","3b","4"),
       lty=c(1:6),col=c("black","red","blue","darkgray","darkgray", "purple"),lwd=2,bty="n")

abline(h=0.1,lty=1,lwd=1)


if (save_pdf) dev.off()

detach(rx)








###########################################################################
## Method 3 - FDR controlled in each fold, not overall ####################
###########################################################################

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
attach(rx)

sdx=0.1 # Gaussian smoothing width
f1=fdp_cf1_fdr3_adj0_dist1_fold1
f2=fdp_cf1_fdr3_adj0_dist1_fold2
f3=fdp_cf1_fdr3_adj0_dist1_fold3
f3b=fdp_cf1_fdr3b_adj0_dist1
fr=fdp_cf1_fdr3_adj0_dist1

xx=(n1p+n1pq)
sdx=sdx*(max(xx)-min(xx))# rescale kernel width

w=which(pmin(fr,f1,f2,f3,f3b) > -0.5)
f1=f1[w]; f2=f2[w]; f3=f3[w]; f3b=f3b[w]; fr=fr[w]; xx=xx[w];


u=sort(unique(xx)); 
f1x=0*u; f2x=0*u; f3x=0*u; f3bx=0*u; frx=0*u; 
for (i in 1:length(u)) {
    wts=dnorm(xx-u[i],sd=sdx); wts=wts/sum(wts)
    f1x[i]=sum(f1*wts); f2x[i]=sum(f2*wts); f3x[i]=sum(f3*wts);
    f3bx[i]=sum(f3b*wts); frx[i]=sum(fr*wts);
}

if (save_pdf) pdf(paste0(output_dir,"fdp_cf1_noci.pdf"),width=6,height=6)

plot(0,xlim=range(xx),ylim=c(0.05,0.25),ylab="FDR",xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),bty="n")


su=sort(u); ou=order(u)
lines(su,f1x[ou],col="blue",lwd=2,lty=3); 
lines(su,f2x[ou],col="blue",lwd=2,lty=3); 
lines(su,f3x[ou],col="blue",lwd=2,lty=3); 
lines(su,f3bx[ou],col="black",lwd=2,lty=1); 
lines(su,frx[ou],col="red",lwd=2,lty=2); 



legend("topright",c("Fold 1-3","Overall","3b", expression(alpha)),lty=c(3,2,1,4),col=c("blue","red","black","black"),bty="n")

abline(h=0.1,col="black",lty=4)

if (save_pdf) dev.off()

detach(rx)




###########################################################################
## Comparison of power by FDR control method ##############################
###########################################################################

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
attach(rx)

if (save_pdf) pdf(paste0(output_dir,"fdp_control.pdf"),width=6,height=6)

normalise=T # set to F to show absolute power rather than relative to fp
txp=tdr_p
tx1=tdr_cf1_fdr1_adj0_dist1
tx3b=tdr_cf1_fdr3b_adj0_dist1
tx4=tdr_cf1_fdr4_adj0_dist1
xx=n1p+n1pq

w=which(pmin(txp,tx1,tx3b,tx4)> -0.5 & # remove simulations which failed to finish in time
          xx>0 & # tdr is indeterminate if xx==0
          !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions
txp=txp[w]; tx1=tx1[w]; tx3b=tx3b[w]; tx4=tx4[w]; xx=xx[w]

ci=0.9 # confidence envelope width
sdx=0.15 # gaussian smoothing
mwd=2 # line width

xf=unique(xx)

fp=rep(0,length(xf));
f1=fp; f3b=fp; f4=fp

sdx=sdx*(max(xx)-min(xx))  #xx=xx/max(xx)

for (i in 1:length(xf)) {
   wts=dnorm(xx-xf[i],sd=sdx); wts=wts/sum(wts)
   fp[i]=sum(txp*wts)
   f1[i]=sum(tx1*wts)
   f3b[i]=sum(tx3b*wts)
   f4[i]=sum(tx4*wts)
}

if (normalise) {
  f1=f1-fp; f3b=f3b-fp; f4=f4-fp; fp=fp-fp
}

prange=range(c(fp,f1,f3b,f4))
if (normalise) yx=expression(paste(Delta,"(power)")) else yx="Power"
plot(0,type="n",xlim=range(xx),ylim=prange,xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=yx,bty="n")

lines(sort(xf),fp[order(xf)],col="black",lty=1,lwd=mwd)
lines(sort(xf),f1[order(xf)],col="black",lty=2,lwd=mwd)
lines(sort(xf),f3b[order(xf)],col="red",lty=3,lwd=mwd)
lines(sort(xf),f4[order(xf)],col="blue",lty=4,lwd=mwd)

legend(300,0.004,
       c("4", "3b","1", "BH"),lty=4:1,lwd=mwd,col=c("blue","red","black","black"),bg="white",bty="n")

if (save_pdf) dev.off()

detach(rx)






###########################################################################
## Assessment of approximation of null distribution #######################
###########################################################################

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
attach(rx)

if (save_pdf) pdf(paste0(output_dir,"fdp_approx_effect.pdf"),width=5,height=5)

fp=fdp_p

xx=n1p+n1pq; rr=range(xx)

plot(0,type="n",xlim=rr,ylim=c(0.07,0.35),ylab="FDR",yaxt="n",
     xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),bty="n")

axat=c(0.07,0.1,0.13)
axis(2,at=axat,label=axat,las=1)
axis(2,at=axat+0.1,label=axat,las=1)
axis(2,at=axat+0.2,label=axat,las=1)



for (ii in 1:3) {

dset=ii
fx=fdp_cf1_fdr4_adj0_dist1
gx=fdp_cf1_fdr4_adj0_dist2
xx=n1p+n1pq
  
w=which(pmin(fx,gx)> -0.5 & (dist1 %in% dset) ) # & (pi0_null+sigma_null!= 1.5) & (pi0_null + sigma_null)>0)

xx=xx[w]; fx=fx[w]; gx=gx[w]; fpx=fp[w]
sdx=0.15

ux=sort(unique(xx)); fy=rep(0,length(ux)); gy=fy; py=fy
sdx=sdx*(max(xx)-min(xx))
for (i in 1:length(ux)) {
    wts=dnorm(xx-ux[i],sd=sdx); wts=wts/sum(wts)
    fy[i]=sum(fx*wts);  gy[i]=sum(gx*wts); py[i]=sum(fpx*wts)
}

lines(sort(ux),fy[order(ux)] + 0.1*(3-ii),lwd=2,lty=2,col="red"); 
lines(sort(ux),gy[order(ux)] + 0.1*(3-ii) ,lwd=2,lty=3,col="blue"); 

}


abline(h=0.1,col="black")
abline(h=0.2,col="black")
abline(h=0.3,col="black")

legend("topright",c("True","Est."),lty=c(2,3),col=c("red","blue"),bty="n")

abline(h=0.15,col="white",lwd=5)
abline(h=0.25,col="white",lwd=5)

text(50,0.33,"Normal")
text(50,0.23,"t (3df)")
text(50,0.13,"Cauchy")


if (save_pdf) dev.off()







###########################################################################
## GWAS analysis ##########################################################
###########################################################################

###########################################################################
## Load data ##############################################################
###########################################################################

load(paste0(cfdr_dir,"data/gwas_data.RData"))

# P-values
p=gwas_results$JIA; q=gwas_results$T1D; 
names(p)=rownames(gwas_results); names(q)=names(p)
zp=-qnorm(p/2); zq=-qnorm(q/2); pq=cbind(p,q)

# Chromosome and position (build 37)
chrx=gwas_results$CHR; posx=gwas_results$POS; 
names(chrx)=names(p); names(posx)=names(p)

# Folds for cross-validation
fold=chrx # that was easy



###########################################################################
## cFDR and FDR control ###################################################
###########################################################################

# Parameters of null distribution
nullq_pars=fit.2g(pq[which(p>0.5),],rho=corr_null)$pars

# File to save in
vfile=paste0(cfdr_dir,"data/gwas_cfdr_data.RData")

# get co-ordinates of l-regions and v-values

if (!file.exists(vfile)) {
  xsub=which(p<5e-4) # only calculate v values for these SNPs, to save time
  
  nv=2000
  vx3x=matrix(-1,length(xsub),4+nv); rownames(vx3x)=names(p)[xsub]
  for (i in 1:max(fold)) {
    foldx=which(fold==i)
    isub=intersect(xsub,foldx)
    if (length(isub) > 0) {
      vxi=vl(p,q,indices=isub,fold=foldx,mode=2,adj=F,nv=nv)
      vx3x[match(isub,xsub),]=vxi$x; 
    }
    print(i)
  }
  
  Lcoords_jia_t1d=list(x=vx3x,y=vxi$y)
  v_jia_t1d=il(Lcoords_jia_t1d,pi0_null=nullq_pars[1],sigma_null=nullq_pars[2],rho_null=corr_null)

  save(xsub,nullq_pars,v_jia_t1d,Lcoords_jia_t1d,file=vfile)
  
} else load(vfile)




###########################################################################
## Hypothesis tests #######################################################
###########################################################################

vsub=v_jia_t1d; vxy=Lcoords_jia_t1d


# FDR control level
alpha= 5e-8 * length(p)/length(which(p< 5e-8))

# Benjamini-Hochberg procedure on p-values
rp=rank(p)
hit_p=names(p)[which(rp<max(rp[which(p/rp < alpha/length(p))]))]

# FDR-control method 3b on cFDR
v=rep(1,length(p)); names(v)=names(p);  v[xsub]=vsub; # v values; set to 1 for values for which L-regions were not computed
rv=rank(v)
hit_c=names(p)[which(rv<max(rv[which(v/rv < alpha/length(v))]))]

# cFDR hits in new regions (no p-value hit within 2MB)
hit_c_unique=c()
for (i in 1:length(hit_c)) {
  mp=which(chrx[hit_p]==chrx[hit_c[i]]) # p-value hits on the same chromosome
  if (length(mp)==0)  hit_c_unique=c(hit_c_unique,hit_c[i])
   else if (min(abs(posx[hit_p[mp]] - posx[hit_c[i]]))> 2e6) 
     hit_c_unique=c(hit_c_unique,hit_c[i]) 
}

# p-value hits lost in cFDR analysis (no cfdr hit within 2MB)
hit_p_unique=c()
for (i in 1:length(hit_p)) {
  mc=which(chrx[hit_c]==chrx[hit_p[i]]) # p-value hits on the same chromosome
  if (length(mc)==0)  hit_p_unique=c(hit_p_unique,hit_p[i])
   else if (min(abs(posx[hit_c[mc]] - posx[hit_p[i]]))> 2e6) 
     hit_p_unique=c(hit_p_unique,hit_p[i]) 
}

## Write to file
hit_c_tab=data.frame("RSID"=hit_c,"CHR"=chrx[hit_c],"POS"=posx[hit_c],
                     "P_JIA"=p[hit_c],"P_T1D"=q[hit_c],"V"=v[hit_c],
                     "UNIQUE_TO_CFDR"=hit_c %in% hit_c_unique)
write.table(hit_c_tab,file=paste0(output_dir,"cfdr_hits_tab.txt"),row.names=F,quote=F)

hit_p_tab=data.frame("RSID"=hit_p,"CHR"=chrx[hit_p],"POS"=posx[hit_p],
                     "P_JIA"=p[hit_p],"P_T1D"=q[hit_p],"V"=v[hit_p],
                     "UNIQUE_TO_PVAL"=hit_p %in% hit_p_unique)
write.table(hit_p_tab,file=paste0(output_dir,"pval_hits_tab.txt"),row.names=F,quote=F)


###########################################################################
## Plot of hits ###########################################################
###########################################################################

if (save_pdf) pdf(paste0(output_dir,"cfdr_jia_t1d_z.pdf"),width=6,height=6)

# don't draw every point near the origin, to cut down on plot size
rlim=3; osub=which(zp^2 + zq^2 > rlim^2)
plot(zp[osub],zq[osub],cex=0.3,xlim=c(0,8),ylim=c(0,10),xaxs="i",yaxs="i",
     main="JIA|T1D",xlab="Z (JIA)",ylab="Z (T1D)",bty="n")
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")

points(zp[hit_c],zq[hit_c],col="red",pch=3,cex=0.5)
points(zp[hit_p],zq[hit_p],col="blue",pch=5,cex=0.5)

rv=rank(v)
mv=which(xsub==which(rv==length(hit_c))); 
psub=3:2003
lines(-qnorm(vxy$x[mv,psub]/2),-qnorm(vxy$y[psub]/2),col="red",lwd=2,lty=1)
abline(v=-qnorm(2.5e-8),col="blue",lwd=2,lty=2)

legend("bottomleft",c("Datapoints", expression(paste("H"[0]," rej. by p-val")),
       expression(paste("H"[0]," rej. by cFDR"))),bg="white",
       col=c("black","blue","red"),pch=c(1,5,3),pt.cex=c(0.3,0.75,0.75))

if (save_pdf) dev.off()



###########################################################################
## Back-to-back Manhattan plot (Supplementary figure) #####################
###########################################################################


if (save_pdf) pdf(paste0(output_dir,"manhattan_jia_t1d.pdf"),width=12,height=8)

# basic function
manhattan=function(chr,pos,p,minp=0,col2=c("blue","red"),add=F,flip=F,...) {
 cpos=pos
 mpos=0
 for (i in 1:max(chr)) {
 cpos[which(chr==i)]=cpos[which(chr==i)]+mpos
 mpos=max(cpos[which(chr==i)])
 }
 if (flip) sgn=-1 else sgn=1
 w=which(-log10(p) >= minp)
 if (!add) plot(cpos[w],-sgn*log10(p)[w],col=col2[1+(chr[w] %% 2)],...) else points(cpos[w],-sgn*log10(p)[w],col=col2[1+(chr[w] %% 2)],...)
 cpos
}

cpost=manhattan(chrx,posx,q,xaxt="n",yaxt="n",minp=3,
               xlab="Chr",ylab=expression(paste("-log"[10],"P/v(L)")),col2=c("gray","darkgray"),
               pch=16,cex=0.5,ylim=c(-15,15),xaxs="i",bty="n"); names(cpost)=names(p)
cpos=manhattan(chrx,posx,p,minp=3,add=T,
               pch=16,cex=0.5); names(cpos)=names(p)
cpos2=manhattan(chrx,posx,v,xaxt="n",add=T,flip=T,
               pch=16,cex=0.5)
atx=c(); for (i in 1:max(chrx)) atx=c(atx,mean(cpos[which(chrx==i)]))
axis(1,1:max(chrx),at=atx)
axis(2,abs(5*(-3:3)),at=5*(-3:3))

lm=-log10(5e-4)
#polygon(c(0,0,max(cpos),max(cpos)),c(-lm,lm,lm,-lm),col="gray")
abline(h=-log10(5e-8))
abline(h=log10(max(v[hit_c])))

s=setdiff(hit_c,hit_p)
ds=c(); for (i in 1:length(s)) ds[i]=min(abs(cpos[hit_p]-cpos[s[i]]))
s=s[which(ds> 1e6)]
points(cpos[s],log10(v[s]),cex=4)

s=setdiff(hit_p,hit_c)
ds=c(); for (i in 1:length(s)) ds[i]=min(abs(cpos[hit_c]-cpos[s[i]]))
s=s[which(ds> 1e6)]
points(cpos[s],-log10(p[s]),cex=4)


if (save_pdf) dev.off()
