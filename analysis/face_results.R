library(scales)
setwd("~/set/your/directory")

# LOAD RESULTS ---------------------------------------------------------------

# create array to load results into
face_res  = array(dim=c(3,5,3,4,100),
                  dimnames = list(
                    paste0("R",c(5,15,30)),
                    paste0("DR",c("1",".8",".6",".4",".2")),
                    c("err","time","tau2"),
                    paste0("C",1:4),
                    paste0("D",1:100)
                  ))
dd = head(dim(face_res),-1)
dn = dimnames(face_res)

# load raw results
face_res_raw = read.csv("P2_faces_chains.csv", header=F)

# arrange raw results into array
for(ir in 1:dd[1])
  for(im in 1:dd[2])
    for(iq in 1:dd[3])
      for(ic in 1:dd[4]){
        iv = c(ir,im,iq,ic)
        cn = 1
        for(i in 1:length(dd))
          cn = cn + (iv[i]-1)*(prod(dd[-c(1:i)]))
        face_res[ir,im,iq,ic,] = face_res_raw[,cn]
      }


# PLOT CHAIN ERROR --------------------------------------------------------

par(mfrow=c(2,2), mar=c(2.5,4,2,0.3), las=1)

tmp = face_res
# setting some graphical parameters
ymm = 120 #ymax
ym = 500 #ymin
yl = "Error" #ylab
wd=0.25 #line weight
bl = c(186.4, 158.1) #HOOI comparison

plot(tmp["R15","DR1","err","C1",], type='l', col=alpha(1,1), ylim=c(ymm,ym), main="R = 15, DR = 1", ylab=yl, cex.main=0.95, xlab="Saved draw", lwd=wd)
lines(tmp["R15","DR1","err","C2",], type='l', col=alpha(2,1), lwd=wd)
lines(tmp["R15","DR1","err","C3",], type='l', col=alpha(3,1), lwd=wd)
lines(tmp["R15","DR1","err","C4",], type='l', col=alpha(4,1), lwd=wd)
legend('topleft', legend = c("Chain draws","Median","HOOI"), lwd=c(1,1.5,1.5), lty=c(1:2,4), seg.len=3, cex=0.8)
abline(h=bl[1], lty='dotdash', lwd=2)
abline(h=apply(tmp["R15","DR1","err",,],1,median), lty=2, lwd=2, col=1:4)
text(x=30, y=bl[1], labels=paste("HOOI error:",bl[1]), pos=1)

plot(tmp["R15","DR.4","err","C1",], type='l', col=alpha(1,1), ylim=c(ymm,ym), main="R = 15, DR = 0.4", ylab="", cex.main=0.95, xlab="Saved draw", lwd=wd)
lines(tmp["R15","DR.4","err","C2",], type='l', col=alpha(2,1), lwd=wd)
lines(tmp["R15","DR.4","err","C3",], type='l', col=alpha(3,1), lwd=wd)
lines(tmp["R15","DR.4","err","C4",], type='l', col=alpha(4,1), lwd=wd)
abline(h=bl[1], lty=2, lwd=2)
abline(h=apply(tmp["R15","DR.4","err",,],1,median), lty=2, lwd=2, col=1:4)
text(x=30, y=bl[1], labels=paste("HOOI error:",bl[1]), pos=1)

plot(tmp["R30","DR1","err","C1",], type='l', col=alpha(1,1), ylim=c(ymm,ym), main="R = 30, DR = 1", ylab="", cex.main=0.95, xlab="Saved draw", lwd=wd)
lines(tmp["R30","DR1","err","C2",], type='l', col=alpha(2,1), lwd=wd)
lines(tmp["R30","DR1","err","C3",], type='l', col=alpha(3,1), lwd=wd)
lines(tmp["R30","DR1","err","C4",], type='l', col=alpha(4,1), lwd=wd)
abline(h=apply(tmp["R30","DR1","err",,],1,median), lty=2, lwd=2, col=1:4)
abline(h=bl[2], lty=2, lwd=2)
text(x=30, y=bl[2], labels=paste("HOOI error:",bl[2]), pos=1)

plot(tmp["R30","DR.4","err","C1",], type='l', col=alpha(1,1), ylim=c(ymm,ym), main="R = 30, DR = 0.4", ylab="", cex.main=0.95, xlab="Saved draw", lwd=wd)
lines(tmp["R30","DR.4","err","C2",], type='l', col=alpha(2,1), lwd=wd)
lines(tmp["R30","DR.4","err","C3",], type='l', col=alpha(3,1), lwd=wd)
lines(tmp["R30","DR.4","err","C4",], type='l', col=alpha(4,1), lwd=wd)
abline(h=apply(tmp["R30","DR.4","err",,],1,median), lty=2, lwd=2, col=1:4)
abline(h=bl[2], lty=2, lwd=2)
text(x=30, y=bl[2], labels=paste("HOOI error:",bl[2]), pos=1)


# CONVERGENCE CHECKS ------------------------------------------------------

Rhat = function(X){
  n = dim(X)[2]/2
  m = dim(X)[1]*2
  XX = cbind(t(X[,1:n]), t(X[,-(1:n)])) # XX is n*m
  # Rhat
  B = n/(m-1)*sum((apply(XX,2,mean)-mean(XX))^2)
  W = (1/m)*sum(apply(XX,2,var))
  V = (n-1)*W/n + B/n
  R = sqrt(V/W)
  # n_eff
  Vt = sapply(0:(n-1), function(t) sum((XX[(t+1):n,]-XX[1:(n-t),])^2)/m/(n-t))
  rhot = 1-Vt/(2*V)
  TT = min(which(sapply(1:(n-2), function(t) t%%2==1 && rhot[t+1]+rhot[t+2]<0)))
  n_eff = m*n/(1+2*sum(rhot[1:TT]))
  return(list(R = R, n_eff = n_eff))
}

face_tab = data.frame(R=c(), DR=c(), quant=c(), Rhat=c(), neff=c(), err=c(), time=c())
face_rat = matrix(NA, nrow=0, ncol=nrow(face_res_raw)*dd[4])
face_ratt = matrix(NA, nrow=0, ncol=nrow(face_res_raw)*dd[4])
for(ir in 1:dd[1])
  for(im in 1:dd[2])
    for(iq in setdiff(1:dd[3],2)){
      iv = c(ir,im,iq)
      Rn = Rhat(face_res[ir,im,iq,,])
      face_tab = rbind(face_tab, data.frame(R=dn[[1]][ir], DR=dn[[2]][im], quant=dn[[3]][iq],
                                            Rhat=Rn$R, n_eff=Rn$n_eff,
                                            err=median(face_res[ir,im,"err",,]), time=median(face_res[ir,im,"time",,])))
      face_rat = rbind(face_rat, c(face_res[ir,im,"err",,]/face_res[ir,1,"err",,]))
      face_ratt = rbind(face_ratt, c(face_res[ir,im,"time",,]/face_res[ir,1,"time",,]))
    }

range(face_tab$Rhat)
range(face_tab$n_eff)
# face_tab used in generating the following plots, but see below for a better tabular summary of each setting


# RELATIVE ERROR/TIME PLOTS -----------------------------------------------

ff = function(x) (x/x[1])
ff_arrows = function(dd,ii,cc,ww){
  tt = apply(dd[which(ii),],1,quantile,probs=c(0.25,0.5,0.75))
  #print(tt)
  for(i in 2:5){
    x0=(6-i)/5
    y0 = tt[2,i]
    arrows(x0=x0, y0=y0, x1=x0, y1=tt[1,i], angle = 90, length=0.05, col=cc, lwd=ww/2)
    arrows(x0=x0, y0=y0, x1=x0, y1=tt[3,i], angle = 90, length=0.05, col=cc, lwd=ww/2)
  }
}

par(mfrow=c(2,1), mar=c(3,4,1.5,0.2), oma=c(0,0,0,0), las=1, cex.lab=1.15, cex.axis=1.1, cex.main=1.25)

plot(ff(face_tab[face_tab$quant=="err" & face_tab$R=="R30",]$err)~c(1,0.8,0.6,0.4,0.2), type='o', ylim=c(0,3), col="royalblue", pch=20, lwd=3, ylab="Error (relative to DR = 1)", xlab="")
lines(ff(face_tab[face_tab$quant=="err" & face_tab$R=="R15",]$err)~c(1,0.8,0.6,0.4,0.2), type='o', pch=20, col="royalblue", lwd=2)
lines(ff(face_tab[face_tab$quant=="err" & face_tab$R=="R5",]$err)~c(1,0.8,0.6,0.4,0.2), type='o', pch=20, col="royalblue", lwd=1)
ff_arrows(face_rat, face_tab$quant=="err" & face_tab$R=="R30", "royalblue", 3)
ff_arrows(face_rat, face_tab$quant=="err" & face_tab$R=="R15", "royalblue", 2)
ff_arrows(face_rat, face_tab$quant=="err" & face_tab$R=="R5", "royalblue", 1)
abline(h=1, lty=2)
title("Face")
title(xlab="DR", line = 2.15)

plot(ff(face_tab[face_tab$quant=="err" & face_tab$R=="R30",]$time)~c(1,0.8,0.6,0.4,0.2), type='o', ylim=c(0,1.5), col="royalblue", pch=20, lwd=3, ylab="Time (relative to DR = 1)", xlab="")
lines(ff(face_tab[face_tab$quant=="err" & face_tab$R=="R15",]$time)~c(1,0.8,0.6,0.4,0.2), type='o', pch=20, col="royalblue", lwd=2)
lines(ff(face_tab[face_tab$quant=="err" & face_tab$R=="R5",]$time)~c(1,0.8,0.6,0.4,0.2), type='o', pch=20, col="royalblue", lwd=1)
ff_arrows(face_ratt, face_tab$quant=="err" & face_tab$R=="R30", "royalblue", 3)
ff_arrows(face_ratt, face_tab$quant=="err" & face_tab$R=="R15", "royalblue", 2)
ff_arrows(face_ratt, face_tab$quant=="err" & face_tab$R=="R5", "royalblue", 1)
abline(h=1, lty=2)
title(xlab="DR", line = 2.15)

legend('bottomright', legend = c("R = 5","R = 15","R = 30"), col=c(1,1,1), lwd=c(1,2,3), cex=1)


# OTHER -------------------------------------------------------------------

face_tab = data.frame(R=c(), DR=c(), quant=c(), Rhat=c(), neff=c(), err=c(), time=c())
for(ir in 1:dd[1])
  for(im in 1:dd[2])
    for(iq in setdiff(1:dd[3],2)){
      iv = c(ir,im,iq)
      Rn = Rhat(face_res[ir,im,iq,,])
      face_tab = rbind(face_tab, data.frame(R=dn[[1]][ir], DR=dn[[2]][im], quant=dn[[3]][iq],
                                            Rhat=Rn$R, n_eff=Rn$n_eff,
                                            err=paste0(format(round(median(face_res[ir,im,"err",,]),1),nsmall=1), " (",format(round(IQR(face_res[ir,im,"err",,]),1),nsmall=1),")"),
                                            time=paste0(format(round(median(face_res[ir,im,"time",,]),2),nsmall=2), " (",format(round(IQR(face_res[ir,im,"time",,]),2),nsmall=2),")")))
    }
face_tab
