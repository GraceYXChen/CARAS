n1 <- 28 / 2
n2 <- 24 / 2
n3 <- 20 / 2
p1 = .5; p2 = .68
R <- 36
D <- 36 + 8
ana_time <- c(13, 25, 44)
baseParam=c(1,1); cen_rate= .0016; p = 2 

a1 <- 3.71; a2 <- 2.511; a3 <- 1.993
b1 <- -0.236; b2 <- 1.170

library(survival)
library(mvtnorm)
library(gsDesign)
library(survminer)
library(rpact)
design <- getDesignGroupSequential(
  # assume one-sided test, equal information increment
  # OF alpha and beta spending function
  sided = 1, alpha = 0.025, beta = 0.2,
  informationRates = c(1/3, 2/3, 1), kMax=3,
  typeOfDesign = "asOF",typeBetaSpending = "bsOF"
)

sampleSizeResult <- getSampleSizeSurvival(
  design = design, lambda2 =1, lambda1 = 0.5,
  dropoutRate1 = .0016, dropoutRate2 = .0016, dropoutTime = 1,
  followUpTime = 8, # the whole duration is the maximum accrual time 3
  accrualTime = c(0,36)
)

sampleSizeResult

out <- NULL
rej = NULL; early_stop_f <- NULL;dif01 <- NULL;early_stop_s <- NULL;dif_ms <- NULL;
ss = NULL; H = NULL; ms <- NULL; fail=NULL
ar = NULL; d_1 = NULL; d_2 = NULL; d_3 = NULL; d_4 = NULL
ms11 <- NULL
ms10 <- NULL
ms21 <- NULL
ms20 <- NULL

ms1n1 <- NULL
ms1n0 <- NULL
ms2n1 <- NULL
ms2n0 <- NULL

ms1 <- NULL; ms0 <- NULL
method="CARA3"
set.seed(2024)
for(iter in 1:1000){
  try({
    n = n1;  X1 = rbinom(n*2,1,p1); X2 = rbinom(n*2,1,p2);
    G = rep(0:1,each=n)
    XX = cbind(X1,X2,G)
    r = runif((n1+n2+n3)*2,0,R) # accrual time
    r1 <- r[which(r<=sort(r)[(n1*2)])]
    r2 <- r[which(r<=sort(r)[((n1+n2)*2)] & r>sort(r)[(n1*2)])]
    r3 <- r[which(r>sort(r)[((n1+n2)*2)])]
    XParam = c(0,0)
    GParam = 0.86
    IntParam = c(-1.45,-1.45)
    data_1 = DataModel(X=XX[,1:p], XParam=XParam,r=r1,
                       G = G, GParam = GParam, IntParam = IntParam,
                       n=n*2, p=p, baseParam=baseParam, sigma=sigma,
                       R=ana_time[1], D=D, cenParam=cen_rate,mod="cox")
    tt1 = data_1$tc
    trc_data1 = data_1$data
    trc_data1$r = r1
    trc_data1 <- trc_data1[trc_data1$r<ana_time[1],]
    
    trc_data1$y = ifelse((trc_data1$r+ trc_data1$y>ana_time[1]), ana_time[1]-trc_data1$r,
                         trc_data1$y)
    trc_data1$del = ifelse(trc_data1$r+ trc_data1$y>ana_time[1], 0,
                           trc_data1$del)
    
    trc_data1_1 = trc_data1[trc_data1$G==0,-6]
    trc_data1_2 = trc_data1[trc_data1$G==1,-6]
    
    Z.stat <- NULL; pv = NULL;
    Z.stat[1] <- owstat(trc_data1)
    if(Z.stat[1] > a1){
      rej <- c(rej,1);early_stop_s <- c(early_stop_s,1)
      early_stop_f <- c(early_stop_f,0); dif01 = c(dif01,0)
      ss <- c(ss, n1*2)
      H <- c(H,exp(-mean(trc_data1$del==0 & trc_data1$y+trc_data1$r>=ana_time[1])))
      ms <- c(ms,NA)
      ms1 <- c(ms1,NA);ms0 <- c(ms0,NA)
      ar <- c(ar,NA)
      dif_ms <- c(dif_ms, NA)
      ms10 <- c(ms10,NA);ms11 <- c(ms11,NA);ms20 <- c(ms20,NA);ms21 <- c(ms21,NA)
      ms1n0 <- c(ms1n0,NA);ms1n1 <- c(ms1n1,NA);ms2n0 <- c(ms2n0,NA);ms2n1 <- c(ms2n1,NA)
      fail <- c(fail,NA)
      d_1 <- c(d_1,NA)
      d_2 <- c(d_2,NA)
      d_3 <- c(d_3,NA)
      d_4 <- c(d_4,NA)

      next}else if(Z.stat[1]<b1){
        rej <- c(rej,0);early_stop_s <- c(early_stop_s,0)
        early_stop_f <- c(early_stop_f,1); dif01 = c(dif01,0)
        ss <- c(ss, n1*2)
        H <- c(H,exp(-mean(trc_data1$del==0 & trc_data1$y+trc_data1$r>=ana_time[1])))
        ms <- c(ms,NA)
        ar <- c(ar,NA)
        ms1 <- c(ms1,NA);ms0 <- c(ms0,NA)
        ms10 <- c(ms10,NA);ms11 <- c(ms11,NA);ms20 <- c(ms20,NA);ms21 <- c(ms21,NA)
        ms1n0 <- c(ms1n0,NA);ms1n1 <- c(ms1n1,NA);ms2n0 <- c(ms2n0,NA);ms2n1 <- c(ms2n1,NA)
        dif_ms <- c(dif_ms, NA)
        fail <- c(fail,NA)
        d_1 <- c(d_1,NA)
        d_2 <- c(d_2,NA)
        d_3 <- c(d_3,NA)
        d_4 <- c(d_4,NA)
     
        next}
    else if((Z.stat[1]<a1)&(Z.stat[1]>b1)){
      # Interim analysis 2
      # generate covariates of the new data
      n = n2
      X1 = rbinom(n,1,p1); X2 = rbinom(n,1,p2);XX1 = cbind(X1,X2) # Group 1 covariates
      X1 = rbinom(n,1,p1); X2 = rbinom(n,1,p2);XX2 = cbind(X1,X2) # Group 2 covariates
      XX = rbind(XX1,XX2)
      trt <- NULL
      tau = min(max(trc_data1_1[trc_data1_1$del==1,-5]$y),
                max(trc_data1_2[trc_data1_2$del==1,-5]$y))
      
      if(method=="RAR"){
        A1 <- Acoef_RAR(trc_data1_1[,-(3:5)],tau)
        A2 <- Acoef_RAR(trc_data1_2[-(3:5)],tau)
        pi <- allocationRatio("Neyman", A1,A2,trc_data1_1[,-(3:5)],trc_data1_2[,-(3:5)],newx=data.frame(X1=0,X2=0))
        trt <- rbinom(n*2,1,1-pi)
      }else if(method=="Trad"){
        trt = rep(c(0:1),n)
      }else{
        try({A1 = c(Acoef(trc_data1_1[,-5],data.frame(X1=0,X2=0), tau),
                    Acoef(trc_data1_1[,-5],data.frame(X1=0,X2=1), tau),
                    Acoef(trc_data1_1[,-5],data.frame(X1=1,X2=0), tau),
                    Acoef(trc_data1_1[,-5],data.frame(X1=1,X2=1), tau))
        A2 = c(Acoef(trc_data1_2[,-5],data.frame(X1=0,X2=0), tau),
               Acoef(trc_data1_2[,-5],data.frame(X1=0,X2=1), tau),
               Acoef(trc_data1_2[,-5],data.frame(X1=1,X2=0), tau),
               Acoef(trc_data1_2[,-5],data.frame(X1=1,X2=1), tau))},silent=T)
        for(i in 1:(2*n)){
          newx = data.frame(t(XX[i,]))
          ind = which(list(c(0,0),c(0,1),c(1,0),c(1,1))%in%list(as.numeric(newx)))
          if(method == "CARA1"){
            try({pi <- allocationRatio("Neyman",A1[ind] ,A2[ind],trc_data1_1[,-5],trc_data1_2[,-5],newx)},T)
          }else if(method == "CARA2"){
            try({pi <- allocationRatio("Survival",A1[ind] ,A2[ind],trc_data1_1[,-5],trc_data1_2[,-5],newx)},T)
          }else if(method=="CARA3"){
            try({pi <- allocationRatio("Hazard",A1[ind] ,A2[ind],trc_data1_1[,-5],trc_data1_2[,-5],newx)},T)
          }
          trt[i] <- rbinom(1,1,1-pi)  # 1 indicates group 1, 0 indicates group 0
        }
      }
      
      data_2 = DataModel(X=XX[,1:p], XParam=XParam,r=r2,
                         G = trt, GParam=GParam, IntParam=IntParam,
                         n=n*2, p=p, baseParam=baseParam, mod = "cox",
                         sigma=sigma, R=12, D=D, cenParam=cen_rate)
      tt2 = data_2$tc
      data_2only = data_2$data
      data_2 = rbind(data_1$data,data_2$data)
      trc_data2 = data_2
      trc_data2$r = c(r1,r2)
      trc_data2 = trc_data2[trc_data2$r<ana_time[2],]
      trc_data2$y = ifelse(trc_data2$r+trc_data2$y>ana_time[2], ana_time[2]-trc_data2$r,
                           trc_data2$y)
      trc_data2$del = ifelse(trc_data2$r+trc_data2$y>ana_time[2], 0,
                             trc_data2$del)
      
      trc_data2_1 = trc_data2[trc_data2$G==0,-6]
      trc_data2_2 = trc_data2[trc_data2$G==1,-6]
      
      
      Z.stat[2] <-  owstat(trc_data2)     
      dif = sum(trt==1)-sum(trt==0)
      
      if(Z.stat[2] > a2){
        rej <- c(rej,1);early_stop_s <- c(early_stop_s,1)
        early_stop_f <- c(early_stop_f,0); dif01 = c(dif01,dif)
        ss <- c(ss, (n1+n2)*2)
        
        fitn <- coxph(Surv(y,del) ~., data=data_2only)
        # baseline survival function
        obj0 <- survfit(fitn,newdata=data.frame(X1=0,X2=0,G=0))
        S0_tau <- summary(obj0,extend=T)$surv
        expn <- exp(predict(fitn, type="lp")) 
        H <- c(H,mean(exp(-S0_tau^(expn))))
        
        ms <- c(ms,dif_median_surv(data_2only))
        ms1 <- c(ms1,dif_median_surv(data_2only[data_2only$G==1,]))
        ms0 <- c(ms0,dif_median_surv(data_2only[data_2only$G==0,]))
        ms10 <- c(ms10,dif_median_surv(data_2only[data_2only$G==0&data_2only$X1==1,]));
        ms11 <- c(ms11,dif_median_surv(data_2only[data_2only$G==1&data_2only$X1==1,]));
        ms20 <- c(ms20,dif_median_surv(data_2only[data_2only$G==0&data_2only$X2==1,]));
        ms21 <- c(ms21,dif_median_surv(data_2only[data_2only$G==1&data_2only$X2==1,]));
        
        ms1n0 <- c(ms1n0,dif_median_surv(data_2only[data_2only$G==0&data_2only$X1==0,]));
        ms1n1 <- c(ms1n1,dif_median_surv(data_2only[data_2only$G==1&data_2only$X1==0,]));
        ms2n0 <- c(ms2n0,dif_median_surv(data_2only[data_2only$G==0&data_2only$X2==0,]));
        ms2n1 <- c(ms2n1,dif_median_surv(data_2only[data_2only$G==1&data_2only$X2==0,]));
        
        fail = c(fail, sum(trc_data2$del==1))
        
        dif_ms <- c(dif_ms,mean(data_2only$y<2 & data_2only$del==1,na.rm=T))
        ar <- c(ar,sum(data_2$G)/sum(data_2$G==0))
        d_1 <- c(d_1,sum(data_2[data_2$X1==1,]$G==1)-sum(data_2[data_2$X1==1,]$G==0))
        d_2 <- c(d_2,sum(data_2[data_2$X1==0,]$G==1)-sum(data_2[data_2$X1==0,]$G==0))
        d_3 <- c(d_3,sum(data_2[data_2$X2==1,]$G==1)-sum(data_2[data_2$X2==1,]$G==0))
        d_4 <- c(d_4,sum(data_2[data_2$X2==0,]$G==1)-sum(data_2[data_2$X2==0,]$G==0))
        
        
        
        next}else if(Z.stat[2]<b2){
          rej <- c(rej,0);early_stop_s <- c(early_stop_s,0)
          early_stop_f <- c(early_stop_f,1); dif01 = c(dif01,dif)
          
          ss <- c(ss, (n1+n2)*2)
          
          fitn <- coxph(Surv(y,del) ~., data=data_2only)
          # baseline survival function
          obj0 <- survfit(fitn,newdata=data.frame(X1=0,X2=0,G=0))
          S0_tau <- summary(obj0,extend=T)$surv
          expn <- exp(predict(fitn, type="lp")) 
          H <- c(H,mean(exp(-S0_tau^(expn))))
          
          ms <- c(ms,dif_median_surv(data_2only))
          ms1 <- c(ms1,dif_median_surv(data_2only[data_2only$G==1,]))
          ms0 <- c(ms0,dif_median_surv(data_2only[data_2only$G==0,]))
          ms10 <- c(ms10,dif_median_surv(data_2only[data_2only$G==0&data_2only$X1==1,]));
          ms11 <- c(ms11,dif_median_surv(data_2only[data_2only$G==1&data_2only$X1==1,]));
          ms20 <- c(ms20,dif_median_surv(data_2only[data_2only$G==0&data_2only$X2==1,]));
          ms21 <- c(ms21,dif_median_surv(data_2only[data_2only$G==1&data_2only$X2==1,]));
          
          ms1n0 <- c(ms1n0,dif_median_surv(data_2only[data_2only$G==0&data_2only$X1==0,]));
          ms1n1 <- c(ms1n1,dif_median_surv(data_2only[data_2only$G==1&data_2only$X1==0,]));
          ms2n0 <- c(ms2n0,dif_median_surv(data_2only[data_2only$G==0&data_2only$X2==0,]));
          ms2n1 <- c(ms2n1,dif_median_surv(data_2only[data_2only$G==1&data_2only$X2==0,]));
          
          d_1 <- c(d_1,sum(data_2[data_2$X1==1,]$G==1)-sum(data_2[data_2$X1==1,]$G==0))
          d_2 <- c(d_2,sum(data_2[data_2$X1==0,]$G==1)-sum(data_2[data_2$X1==0,]$G==0))
          d_3 <- c(d_3,sum(data_2[data_2$X2==1,]$G==1)-sum(data_2[data_2$X2==1,]$G==0))
          d_4 <- c(d_4,sum(data_2[data_2$X2==0,]$G==1)-sum(data_2[data_2$X2==0,]$G==0))
          
          
          fail = c(fail, sum(trc_data2$del==1))
          ar <- c(ar,sum(data_2$G)/sum(data_2$G==0))
          
          dif_ms <- c(dif_ms,mean(data_2only$y<2 & data_2only$del==1,na.rm=T))
          next}else if((Z.stat[2]<a2)&(Z.stat[2]>b2)){
            early_stop_s <- c(early_stop_s,0)
            early_stop_f <- c(early_stop_f,0)
            # Final Analysis
            # generate covariates of the new data
            n = n3
            X1 = rbinom(n,1,p1); X2 = rbinom(n,1,p2);XX1 = cbind(X1,X2) # Group 1 covariates
            X1 = rbinom(n,1,p1); X2 = rbinom(n,1,p2);XX2 = cbind(X1,X2) # Group 2 covariates
            XX = rbind(XX1,XX2)
            trt <- NULL
            tau = min(max(trc_data2_1[trc_data2_1$del==1,-5]$y),
                      max(trc_data2_2[trc_data2_2$del==1,-5]$y))
            
            if(method=="RAR"){
              A1 <- Acoef_RAR(trc_data2_1[,-(3:5)],tau)
              A2 <- Acoef_RAR(trc_data2_2[-(3:5)],tau)
              pi <- allocationRatio("Neyman", A1,A2,trc_data2_1[,-(3:5)],trc_data2_2[,-(3:5)],newx=data.frame(X1=0,X2=0))
              trt <- rbinom(n*2,1,1-pi)
            }else if(method=="Trad"){
              trt = rep(c(0:1),n)
            }else{
              try({A1 = c(Acoef(trc_data2_1[,-5],data.frame(X1=0,X2=0), tau),
                          Acoef(trc_data2_1[,-5],data.frame(X1=0,X2=1), tau),
                          Acoef(trc_data2_1[,-5],data.frame(X1=1,X2=0), tau),
                          Acoef(trc_data2_1[,-5],data.frame(X1=1,X2=1), tau))
              A2 = c(Acoef(trc_data2_2[,-5],data.frame(X1=0,X2=0), tau),
                     Acoef(trc_data2_2[,-5],data.frame(X1=0,X2=1), tau),
                     Acoef(trc_data2_2[,-5],data.frame(X1=1,X2=0), tau),
                     Acoef(trc_data2_2[,-5],data.frame(X1=1,X2=1), tau))},silent=T)
              
              for(i in 1:(2*n)){
                newx = data.frame(t(XX[i,]))
                ind = which(list(c(0,0),c(0,1),c(1,0),c(1,1))%in%list(as.numeric(newx)))
                if(method == "CARA1"){
                  try({pi <- allocationRatio("Neyman",A1[ind] ,A2[ind],trc_data2_1[,-5],trc_data2_2[,-5],newx)},T)
                }else if(method == "CARA2"){
                  try({pi <- allocationRatio("Survival",A1[ind] ,A2[ind],trc_data2_1[,-5],trc_data2_2[,-5],newx)},T)
                }else if(method=="CARA3"){
                  try({pi <- allocationRatio("Hazard",A1[ind] ,A2[ind],trc_data2_1[,-5],trc_data2_2[,-5],newx)},T)
                }
                trt[i] <- 1-rbinom(1,1,pi)  # 1 indicates group 1, 0 indicates group 0
              }
            }
            dif = dif + sum(trt==1)-sum(trt==0)
            
            data_3 = DataModel(X=XX[,1:p], XParam=XParam,r=r3,
                               G = trt, GParam=GParam, IntParam=IntParam,
                               n=n*2, p=p, baseParam=baseParam, sigma =sigma,
                               R=19, D=D, cenParam=cen_rate, mod="cox")      
            tt3 = data_3$tc
            data_23 = rbind(data_2only,data_3$data)
            data_3 = rbind(data_2,data_3$data)
            
            trc_data3 = data_3
            
            Z.stat[3] <-  owstat(trc_data3)   
            #dif_ms <- c(dif_ms,dif_median_surv(data_3))
            #dif_ms <- c(dif_ms, mean((data_3$del==1 & tt>=1)))
            #dif_ms <- c(dif_ms,RISCA:::rmst(data_3$y,trc_data3[trc_data3$G==1,]$del,0.5,"s"))
            #dif_ms <- c(dif_ms, mean(c(tt2,tt3)>3,na.rm=T))
            
            dif_ms <- c(dif_ms,mean(data_23$y<2 & data_23$del==1,na.rm=T))
            
            fitn <- coxph(Surv(y,del) ~., data=data_23)
            # baseline survival function
            obj0 <- survfit(fitn,newdata=data.frame(X1=0,X2=0,G=0))
            S0_tau <- summary(obj0,extend=T)$surv
            expn <- exp(predict(fitn, type="lp")) 
            H <- c(H,mean(exp(-S0_tau^(expn))))
            
            dif01 <- c(dif01,dif)
            ss <- c(ss, NA)
            
            ms <- c(ms,dif_median_surv(data_23))
            ms1 <- c(ms1,dif_median_surv(data_23[data_23$G==1,]))
            ms0 <- c(ms0,dif_median_surv(data_23[data_23$G==0,]))
            ms10 <- c(ms10,dif_median_surv(data_23[data_23$G==0&data_23$X1==1,]));
            ms11 <- c(ms11,dif_median_surv(data_23[data_23$G==1&data_23$X1==1,]));
            ms20 <- c(ms20,dif_median_surv(data_23[data_23$G==0&data_23$X2==1,]));
            ms21 <- c(ms21,dif_median_surv(data_23[data_23$G==1&data_23$X2==1,]));
            
            ms1n0 <- c(ms1n0,dif_median_surv(data_23[data_23$G==0&data_23$X1==0,]));
            ms1n1 <- c(ms1n1,dif_median_surv(data_23[data_23$G==1&data_23$X1==0,]));
            ms2n0 <- c(ms2n0,dif_median_surv(data_23[data_23$G==0&data_23$X2==0,]));
            ms2n1 <- c(ms2n1,dif_median_surv(data_23[data_23$G==1&data_23$X2==0,]));
            
            ar <- c(ar,sum(data_3$G)/sum(data_3$G==0))
            d_1 <- c(d_1,sum(data_3[data_3$X1==1,]$G==1)-sum(data_3[data_3$X1==1,]$G==0))
            d_2 <- c(d_2,sum(data_3[data_3$X1==0,]$G==1)-sum(data_3[data_3$X1==0,]$G==0))
            d_3 <- c(d_3,sum(data_3[data_3$X2==1,]$G==1)-sum(data_3[data_3$X2==1,]$G==0))
            d_4 <- c(d_4,sum(data_3[data_3$X2==0,]$G==1)-sum(data_3[data_3$X2==0,]$G==0))
            
            fail = c(fail, NA)
            if(Z.stat[3]>a3){rej=c(rej,1)}else{rej=c(rej,0)}
          }
    }},silent=T)
}

print(data.frame(rej=mean(rej),early_stop_s=mean(early_stop_s),
                 early_stop_f=mean(early_stop_f),dif01=mean(dif01,na.rm=T),
                 deathr=mean(dif_ms,na.rm=T),H=mean(H,na.rm=T),
                 ms = mean(ms,na.rm=T),fail=mean(fail,na.rm=T),
                 ms1 =  mean(ms1,na.rm=T),ms0 =  mean(ms0,na.rm=T),
                 msX1G0 = mean(ms10,na.rm=T), msX1G1 = mean(ms11,na.rm=T),
                 msX2G0 = mean(ms20,na.rm=T), msX2G1 = mean(ms21,na.rm=T),
                 msX1nG0 = mean(ms1n0,na.rm=T), msX1nG1 = mean(ms1n1,na.rm=T),
                 msX2nG0 = mean(ms2n0,na.rm=T), msX2nG1 = mean(ms2n1,na.rm=T)))
mean(ar,na.rm=T)
mean(d_1,na.rm=T)
mean(d_2,na.rm=T)
mean(d_3,na.rm=T)
mean(d_4,na.rm=T)
