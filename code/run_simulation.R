###########################################################################
##                                                                       ##
## Assess improvements to false-discovery rate control in cFDR method    ##
##                                                                       ##
## James Liley, 15/2/18                                                  ##
##                                                                       ##
###########################################################################

###########################################################################
## Directory to write results to ##########################################
###########################################################################

# Assumes current working directory is folder containing file ./data, 
#  ./code, ./output and ./simulations. An example output with seed=1 is 
#  shown in ./simulations.
out_dir="./simulations/"

###########################################################################
## Packages and scripts ###################################################
###########################################################################

library(cfdr)
library(mnormt)
library(MASS)
library(fields)
library(matrixStats)


###########################################################################
## Simulation parameters ##################################################
###########################################################################

if (!exists("seed")) { #choose seed based on system time if it is not already set
  options(digits.secs=8)
  seed=as.numeric(substr(Sys.time(),21,27))
}
set.seed(seed) # random seed

distx=sample(3,1)  # 1 for normal, 2 for t (3df), 3 for Cauchy

par_cont=runif(1);
if (par_cont>0.5) { # select pars from a discrete distribution

 alpha= 0.1 # universal
 nsnp=sample(c(1000,10000),1) # number of variables

 n1p=sample(c(0,10,200),1) # number of variables associated ONLY with P (principal)
 n1q=sample(c(0,10,200),1) # number of variables associated ONLY with Q (conditional)
 n1pq=sample(c(0,10,200),1) # number of variables associated with both

 sp=sample(c(1.5,3),1) # scale of effect sizpe distribution for associations with P
 sq=sample(c(1.5,3),1) # scale of effect sizpe distribution for associations with Q

} else {

 alpha= 0.1 
 nsnp=round(10^runif(1,3,4)) 
 
 n1p=sample(200,1) # number of variables associated ONLY with P (principal)
 n1q=sample(200,1) # number of variables associated ONLY with Q (conditional)
 n1pq=sample(200,1) # number of variables associated with both

 sp=runif(1,1.5,3) # scale of effect sizpe distribution for associations with A
 sq=runif(1,1.5,3) # scale of effect sizpe distribution for associations with B

}



###########################################################################
## Simulation data ########################################################
###########################################################################

if (distx==1) {
  zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(nsnp-n1p-n1pq,sd=1))
  zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(nsnp-n1p-n1pq-n1q,sd=1))
}
if (distx==2) {
  zp=c(rt(n1pq,df=3)*sp,rt(n1p,df=3)*sp,rnorm(nsnp-n1p-n1pq,sd=1))
  zq=c(rt(n1pq,df=3)*sq,rnorm(n1p,sd=1),rt(n1q,df=3)*sq,rnorm(nsnp-n1p-n1pq-n1q,sd=1))
}
if (distx==3) {
  zp=c(rcauchy(n1pq,sc=sp),rcauchy(n1p,sc=sp),rnorm(nsnp-n1p-n1pq,sd=1))
  zq=c(rcauchy(n1pq,sc=sq),rnorm(n1p,sd=1),rcauchy(n1q,sc=sq),rnorm(nsnp-n1p-n1pq-n1q,sd=1))
}

# P-values
p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

mp=min(p[which(p>0)]); p[which(p==0)]=mp
mq=min(q[which(q>0)]); q[which(q==0)]=mq

zp=-qnorm(p/2); zq=-qnorm(q/2)

# H (hypothesis indicators)
if (n1pq+n1p > 0) h1a=1:(n1pq+n1p) else h1a=c()
h0a=setdiff(1:nsnp, h1a)

# Fit four-groups model
pi1_a=n1p/nsnp; pi1_b=n1q/nsnp; pi1_pl=n1pq/nsnp
p1x=max(min(1-pi1_a-pi1_b-pi1_pl,0.99),0.005)
mxit=1000
fitx=fit.4g(cbind(zp,zq),maxit=mxit)
pars=fitx$pars
conv=dim(fitx$hist)[1] < mxit

# Folds for method 3
nfold=3
fold=rep(1:nfold,1+floor(nsnp/nfold))[1:nsnp][order(runif(nsnp))]

# Parmameters for underlying distribution of zq|HP0
ff=tryCatch(fit.2g(q[which(p>0.5)]),
            error=function(e) list(pars=c(0.5,1)), warning=function(w) list(pars=c(0.5,1)))

pi0_null = c(1-(n1q/(nsnp-n1p-n1pq)),ff$pars[1])
sigma_null = c(sq,ff$pars[2])
dist_null = c(distx,1)

# Potential rejections
subx=which(log10(p)+log10(q) < -2)
cf=cfdr(p,q,sub=subx)
sub=which(cf<6*alpha)
if (length(sub)<30) sub=order(cf)[1:30]


###########################################################################
## Save procedure #########################################################
###########################################################################


vars=c("hit_p",as.vector(
  outer(as.vector(
    outer(as.vector(
     outer(as.vector(
         outer(0:1,1:5,function(x,y) paste0("_fdr",y,"_adj",x))),
         1:2,function(x,y) paste0(x,"_dist",y))),
    1:3,function(x,y) paste0("hit_cf",y,x))),
    c("",paste0("_fold",1:nfold)),function(x,y) paste0(x,y))))
vars=vars[which(
  (grepl("fdr3",vars) & grepl("fold",vars)) |
    (!grepl("fold",vars)))]
vars=gsub("fdr5","fdr3b",vars)

save_vec=function() {

fdp_vec=c()
t2r_vec=c()

for (i in 1:length(vars)) {
  if (exists(vars[i])) {
   hitx=get(vars[i])

   # false discovery proportions
   if (length(hitx)>0) {
     if (hitx[1]>0) {
       fdp=length(intersect(hitx,h0a))/length(hitx) 
     } else fdp=-1
   } else fdp=0 # false-discovery proportion by p-value
   fdp_vec=c(fdp_vec,fdp)
   
   # type-2 error rate
   t2r=length(intersect(hitx,h1a))/length(h1a)
   t2r_vec=c(t2r_vec,t2r)
  } else {
    fdp_vec=c(fdp_vec,-1)
    t2r_vec=c(t2r_vec,-1)
  }
}

state_vec=c(seed,alpha,nsnp,dist_null,n1p,n1q,n1pq,sp,sq,pars,conv,ff$pars)
names(state_vec)=c("seed","alpha","N","dist1","dist2","n1p","n1q","n1pq","sp","sq",
                   paste0("fit_",c("pi0","pi1","pi2","tau1","tau2","s1","s2","conv")),"pi0_null","sigma_null")
names(fdp_vec)=gsub("hit","fdp",vars)
names(t2r_vec)=gsub("hit","t2r",vars)

out=c(state_vec,fdp_vec,t2r_vec)

write(out,file=paste0(out_dir,"cfdrsim",seed,".txt"),ncol=length(out))
}


###########################################################################
## Rejections with p-value: Benjamini-Hochberg procedure ##################
###########################################################################

# Rejections with p-value: Benjamini-Hochberg procedure
rp=rank(p)
hit_p0=which(p< alpha*rank(p)/length(p))
if (length(hit_p0)> 0) hit_p=which(rp <= max(rp[hit_p0])) else hit_p=c()


###########################################################################
## Standard non-parametric FDR ############################################
###########################################################################

# Needed several times
vLf=vl(p,q,indices=sub,mode=0,adj=F); 
vLt=vl(p,q,indices=sub,mode=0,adj=T); 

iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

### FDR control method 1
# No adjustment for Pr(H0|Q), using true null distribution
hit_cf1_fdr1_adj0_dist1=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

# using estimated null distribution
hit_cf1_fdr1_adj0_dist2=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf1_fdr1_adj1_dist1=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf1_fdr1_adj1_dist2=({
  cf=cfdr(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})


save_vec()


### FDR control method 2
# No adjustment for Pr(H0|Q)
hit_cf1_fdr2_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr2_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()



# Adjustment for Pr(H0|Q)
hit_cf1_fdr2_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr2_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()




### FDR control method 3/3b
# No adjustment for Pr(H0|Q)
hit_cf1_fdr3_adj0_dist1 = ({
  hit=c()
  vv_all=rep(1,length(p)) 
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=F)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    vv=rep(1,length(p)); vv[sub2]=vvx[,1]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj0_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf1_fdr3b_adj0_dist1",hit_all) 
  hit
})

hit_cf1_fdr3_adj0_dist2 = ({
  hit=c()
  vv_all=rep(1,length(p)) 
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=F)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    vv=rep(1,length(p)); vv[sub2]=vvx[,2]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj0_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf1_fdr3b_adj0_dist2",hit_all) 
  hit
})


save_vec()

# Adjustment for Pr(H0|Q)
hit_cf1_fdr3_adj1_dist1 = ({
  vv_all=rep(1,length(p)) 
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=T)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)    
    vv=rep(1,length(p)); vv[sub2]=vvx[,1]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj1_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf1_fdr3b_adj1_dist1",hit_all) 
  hit
})
save_vec()

hit_cf1_fdr3_adj1_dist2 = ({
  vv_all=rep(1,length(p)) 
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    sub2=intersect(inf,sub)
    vl0=vl(p,q,indices=sub2,mode=2,fold=inf,adj=T)
    vvx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)    
    vv=rep(1,length(p)); vv[sub2]=vvx[,2]; vv=vv[inf]
    rv=rank(vv)
    out=vv*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf1_fdr3_adj1_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf1_fdr3b_adj1_dist2",hit_all) 
  hit
})
save_vec()





### FDR control method 4
# No adjustment for Pr(H0|Q)
vLf=vl(p,q,indices=sub,adj=F,mode=1)
iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

hit_cf1_fdr4_adj0_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr4_adj0_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()



# Adjustment for Pr(H0|Q)
vLt=vl(p,q,indices=sub,adj=T,mode=1)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)

hit_cf1_fdr4_adj1_dist1=({
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()

hit_cf1_fdr4_adj1_dist2=({
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(p)/rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})
save_vec()




###########################################################################
## Parametric cFDR, KDE ###################################################
###########################################################################

## general
vLf=vly(p,q,indices=sub,adj=F); 
vLt=vly(p,q,indices=sub,adj=T); 

iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)



### FDR control method 1
# No adjustment for Pr(H0|Q)
hit_cf2_fdr1_adj0_dist1=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf2_fdr1_adj0_dist2=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf2_fdr1_adj1_dist1=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf2_fdr1_adj1_dist2=({
  cf=cfdry(p,q,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

save_vec()


### FDR control method 2
# No adjustment for Pr(H0|Q)
hit_cf2_fdr2_adj0_dist1=({
 vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf2_fdr2_adj0_dist2=({
 vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf2_fdr2_adj1_dist1=({
 vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf2_fdr2_adj1_dist2=({
 vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

save_vec()


### FDR control method 3
# No adjustment for Pr(H0|Q)
hit_cf2_fdr3_adj0_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=F)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,1])
    out=vv[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj0_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  assign("hit_cf2_fdr3b_adj0_dist1",hit_all)
  hit
})

hit_cf2_fdr3_adj0_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=F)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,2])
    out=vv[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj0_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf2_fdr3b_adj0_dist2",hit_all) 
  hit
})


# Adjustment for Pr(H0|Q)
hit_cf2_fdr3_adj1_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=T)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,1])
    out=vv[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj1_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf2_fdr3b_adj1_dist1",hit_all) 
  hit
})

hit_cf2_fdr3_adj1_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    vl0=vly(p,q,mode=2,indices=inf,fold=inf,adj=T)
    vv=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vv[,2])
    out=vv[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()    
    assign(paste0("hit_cf2_fdr3_adj1_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vv[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf2_fdr3b_adj1_dist2",hit_all) 
  hit
})

save_vec()


### FDR control method 4
# Deprecated

hit_cf2_fdr4_adj0_dist1=-1
hit_cf2_fdr4_adj0_dist2=-1
hit_cf2_fdr4_adj1_dist1=-1
hit_cf2_fdr4_adj1_dist2=-1

save_vec()






###########################################################################
## Parametric cFDR, four-groups ###########################################
###########################################################################

## general
vLf=vlx(p,q,pars,indices=sub,adj=F); 
vLt=vlx(p,q,pars,indices=sub,adj=T); 

iLf=il(vLf,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
iLt=il(vLt,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)


### FDR control method 1
# No adjustment for Pr(H0|Q)
hit_cf3_fdr1_adj0_dist1=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,1])
  out=ccut*iLf[,1]/vM; out[which(iLf[,1]==0)]=0; out[which(iLf[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf3_fdr1_adj0_dist2=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLf$x[,which(vLf$y>0.1 & vLf$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLf[,2])
  out=ccut*iLf[,2]/vM; out[which(iLf[,2]==0)]=0; out[which(iLf[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})
save_vec()


# Adjustment for Pr(H0|Q)
hit_cf3_fdr1_adj1_dist1=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)])  # largest rectangle contained in L
  vM=pmin(vM,iLt[,1])
  out=ccut*iLt[,1]/vM; out[which(iLt[,1]==0)]=0; out[which(iLt[,1]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

hit_cf3_fdr1_adj1_dist2=({
  cf=cfdrx(p,q,pars,adj=F); ccut=cf[sub]
  vM=rowMins(vLt$x[,which(vLt$y>0.1 & vLt$y<1)]) # largest rectangle contained in L
  vM=pmin(vM,iLt[,2])
  out=ccut*iLt[,2]/vM; out[which(iLt[,2]==0)]=0; out[which(iLt[,2]==1)]=1
  
  vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  if (any(out<alpha)) hit=which(rv <= max(rv[sub[which(out<alpha)]])) else hit=c()
  hit
})

save_vec()



### FDR control method 2
# No adjustment for Pr(H0|Q)
hit_cf3_fdr2_adj0_dist1=({
 vr=rep(1,length(p)); vr[sub]=iLf[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf3_fdr2_adj0_dist2=({
 vr=rep(1,length(p)); vr[sub]=iLf[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

save_vec()

# Adjustment for Pr(H0|Q)
hit_cf3_fdr2_adj1_dist1=({
 vr=rep(1,length(p)); vr[sub]=iLt[,1]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

hit_cf3_fdr2_adj1_dist2=({
 vr=rep(1,length(p)); vr[sub]=iLt[,2]; rv=rank(vr)
  out=vr*length(vr)/rv
  if (any(out<alpha)) hit=which(rv <= max(rv[which(out<alpha)])) else hit=c()
  hit
})

save_vec()




### FDR control method 3. 
# No adjustment for Pr(H0|Q)
parsx=matrix(0,nfold,7)
hit_cf3_fdr3_adj0_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=F)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj0_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf3_fdr3b_adj0_dist1",hit_all) 
  hit
})

hit_cf3_fdr3_adj0_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    #parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=F)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj0_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  assign("hit_cf3_fdr3b_adj0_dist2",hit_all) 
  hit
})

# Adjustment for Pr(H0|Q)
hit_cf3_fdr3_adj1_dist1 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    #parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=T)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj1_dist1_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  assign("hit_cf3_fdr3b_adj1_dist1",hit_all) 
  hit
})

hit_cf3_fdr3_adj1_dist2 = ({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    #parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars)$pars #only need to do this once per fold
    vl0=vlx(p,q,indices=inf,pars=parsx[i,],adj=T)
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    assign(paste0("hit_cf3_fdr3_adj1_dist2_fold",i),hitx)
    hit=c(hit,hitx)
    vv_all[inf]=vx[,2] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  assign("hit_cf3_fdr3b_adj1_dist2",hit_all)
  hit
})


### FDR control method 4. For method 2, this corresponds to method 3b applied to local FDR
# Deprecated. Corresponds to LOCAL cfdr.

hit_cf3_fdr4_adj0_dist1=({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars,sgm=c(1,1000))$pars #only need to do this once per fold
    vl0=vlxl(p,q,indices=inf,pars=parsx[i,])
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,1])
    out=vx[,1]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c()
  hit_all
})

hit_cf3_fdr4_adj0_dist2=({
  vv_all=c()
  hit=c()
  for (i in 1:nfold) {
    inf=which(fold==i); outf=which(fold!=i)
    lf=length(inf)
    parsx[i,]=fit.4g(cbind(zp[outf],zq[outf]),pars=pars,sgm=c(1,1000))$pars #only need to do this once per fold
    vl0=vlxl(p,q,indices=inf,pars=parsx[i,])
    vx=il(vl0,pi0_null=pi0_null,sigma_null=sigma_null,dist=dist_null)
    rv=rank(vx[,2])
    out=vx[,2]*lf/rv
    if (any(out<alpha)) hitx=inf[which(rv <= max(rv[which(out<alpha)]))] else hitx=c()
    hit=c(hit,hitx)
    vv_all[inf]=vx[,1] 
  } 
  rva=rank(vv_all); out_all=vv_all*length(p)/rva 
  if (any(out_all<alpha)) hit_all=which(rva <= max(rva[which(out_all<alpha)])) else hit_all=c() 
  hit_all
})


# following is deprecated
hit_cf3_fdr4_adj1_dist1=-1
hit_cf3_fdr4_adj1_dist2=-1

save_vec()



