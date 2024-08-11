beta = .2; alpha = .05; beta_H1 = -0.6
alpha_vec = c(0.0007,0.0164,0.0500)
beta_vec = c(0.0264,0.1165,0.2000)
info_time = c(1/3,2/3,1)
n1 = 33; n2 = 19; n3=17

# probability of positive biomarker for X1 and X2
p1 = .5; p2 = .5

# generated data from the model
mod = "AFT"
mod = "cox"

GParam = 0
XParam = c(0,0)
IntParam = c(0,0)
sigma=1
GParam_all = c(rep(0,4),rep(-0.6,5),rep(0,3),rep(-0.6,6))
XParam_all = matrix(c(0,0, 0.8,0, 1.6,0,
                      2,0, 0,0, 0.8,0, 1.6,0, 
                      0,0, 0.8,0, 0,0,
                      1.5,0.8, 2,1, 0,0,
                      0.5,0.3, 0.8,0.5, 0,0, 0,0, .3,.5),byrow = T,ncol=2)
IntParam_all = matrix(c(rep(0,14),-.4,0,-.4,0,
                        rep(0,12),-.2,-0.3, -.5,-0.8,
                        -.5,-.8),byrow=T,ncol=2)
XParam_all = - XParam_all
GParam_all = - GParam_all
IntParam_all = - IntParam_all


baseParam=c(1,1); cen_rate= .1; p = 2; R = 1  # recruitment duration
D=3  # maximum duration
shape = c(0.2,0.5,1,2,5)

method = "CARA1"
mod="cox"
out = NULL
p1=0.7
lamb=c(1/2,1,2,4,8)
for(sc in c(18:18)){
  print(out)
  GParam = GParam_all[sc]
  XParam = XParam_all[sc,]
  IntParam = IntParam_all[sc,]
  # IntParam=c(-0.2,0);XParam=c(0.8,0);GParam=-0.6
  # baseParam=c(1,lamb[sc])
  rej = NULL; early_stop_f <- NULL;dif01 <- NULL;early_stop_s <- NULL;dif_ms <- NULL
  for(iter in 1:500){
    print(c(mean(rej),mean(early_stop_f),mean(early_stop_s),
            mean(dif01,na.rm=T),mean(dif_ms,na.rm=T)))
    n = n1;  X1 = rbinom(n*2,1,p1); X2 = rbinom(n*2,1,p2);
    #G = rep(1,2*n); G[sample(1:(2*n),n,replace=F)]=0
    G = rep(0:1,n)
    
    # G = rep(1,2*n)
    # G[X1==0&X2==0][1:(length(which(X1==0&X2==0))/2)]=0
    # G[X1==1&X2==0][1:(length(which(X1==1&X2==0))/2)]=0
    # G[X1==1&X2==1][1:(length(which(X1==1&X2==1))/2)]=0
    # G[X1==0&X2==1][1:(length(which(X1==0&X2==1))/2)]=0
    XX = cbind(X1,X2,G)
   # r = pmin(rexp((n1+n2+n3)*2,1),D)
     r = runif((n1+n2+n3)*2,0,D)
    # r = rbeta((n1+n2+n3)*2,0.5,0.5)*D   # Accrual time
    data_1 = DataModel(X=XX[,1:p], XParam=XParam,r=sort(r)[1:(n1*2)],
                       G = G, GParam = GParam, IntParam = IntParam,
                       n=n*2, p=p, baseParam=baseParam, sigma=sigma,
                       R=1.47, D=D, cenParam=cen_rate,mod=mod)
    #data_1 = AFTModel(XX[,1:p], XParam=XParam, G=G, GParam, IntParam, n=n*2, p, sigma, R, D, cenParam=cen_rate)
    t_c = 1.47; tc1 = data_1$tc # calendar time of survival 
    data_1 = data_1$data;  trc_data1 = data_1
    trc_data1$y = ifelse(tc1>t_c, data_1$y-(tc1 - t_c),data_1$y)
    trc_data1$del = ifelse(tc1>t_c,0,1)
    
    trc_data1_1 = trc_data1[trc_data1$G==0,]
    trc_data1_2 = trc_data1[trc_data1$G==1,]
    
    
    #fit1 <- coxph(Surv(y,del)~G,data=trc_data1)
    Z.stat <- NULL; pv = NULL;
    Z.stat[1] <- owstat(trc_data1)
    #Z.stat[1] <-  summary(fit1)$coefficients[1,2]; se1 <- summary(fit1)$coefficients[1,3]
    # Z.stat
    #try({Z.stat[1] <- owstat(trc_data1)},silent=T)
    #a1 <- 0.283; a2 <- 0.55; a3=0.68; b1=1.178; b2=0.774
    a1 = 3.2; a2 = 2.141; a3=1.695; b1=-0.415; b2=0.917
    # a1 <- -qnorm(1-alpha_vec[1])
    # b1 = -qnorm(beta_vec[1]) + beta_H1 / se1
    #a1 = -3.2; a2=-2.141; a3=-1.695; b1=0.415; b2=-0.917
    if(Z.stat[1] > a1){
      rej <- c(rej,1);early_stop_s <- c(early_stop_s,1)
      early_stop_f <- c(early_stop_f,0); dif01 = c(dif01,0)
      dif_ms <- c(dif_ms,dif_median_surv(trc_data1)); 
      #a = survfit(Surv(y,del)~1,trc_data1)
      #dif_ms <- c(dif_ms,a$surv[which(a$time<1)[length(which(a$time<1))]]);
      # res <- summary(survfit(Surv(y,del)~1,trc_data1))
      # dif_ms <- c(dif_ms,rmst(res$time,res$surv,type="l",max.time=1))
      next}else if(Z.stat[1]<b1){rej <- c(rej,0);early_stop_s <- c(early_stop_s,0)
      early_stop_f <- c(early_stop_f,1); dif01 = c(dif01,0)
      dif_ms <- c(dif_ms,dif_median_surv(trc_data1));
      #a = survfit(Surv(y,del)~1,trc_data1)
      #dif_ms <- c(dif_ms,a$surv[which(a$time<1)[length(which(a$time<1))]]);
      # res <- summary(survfit(Surv(y,del)~1,trc_data1))
      # dif_ms <- c(dif_ms,rmst(res$time,res$surv,type="l",max.time=1))
      next}else if((Z.stat[1]<a1)&(Z.stat[1]>b1)){
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
          trt = rep(0:1,each=n)
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
        
        data_2 = DataModel(X=XX[,1:p], XParam=XParam,r=sort(r)[(2*n1+1):(2*n1+2*n2)],
                           G = trt, GParam=GParam, IntParam=IntParam,
                           n=n*2, p=p, baseParam=baseParam, mod = mod,
                           sigma=sigma, R=0.81, D=D, cenParam=cen_rate)
        #data_2 = AFTModel(XX[,1:p], XParam, G=trt, GParam, IntParam, n=n*2, p, sigma, R, D, cenParam=cen_rate)
        tc2 = data_2$tc
        data_2 = data_2$data
        t_c2 = 0.81    # follow the second cohort for 1 year and do interim analysis
        
        trc_data2 = data_2
        trc_data2$y = ifelse(tc2>2.28, data_2$y-(tc2 - 2.28),data_2$y)
        trc_data2$del = ifelse(tc2>2.28,0,trc_data2$del)
        
        trc_data2 <- rbind(trc_data2,
                           data.frame(y=ifelse(tc1>t_c2+t_c,data_1$y + (t_c2+t_c-tc1),data_1$y),
                                      del = ifelse(tc1>t_c2+t_c,0,data_1$del),
                                      X1 = data_1$X1,X2=data_1$X2,G=data_1$G))
        trc_data2_1 = trc_data2[trc_data2$G==0,]
        trc_data2_2 = trc_data2[trc_data2$G==1,]
        
        #fit <- coxph(Surv(y,del)~G,data=trc_data2)
        #
        Z.stat[2] <- owstat(trc_data2)
        #Z.stat[2] <-  summary(fit)$coefficients[1,2];se2 <- summary(fit)$coefficients[1,3]
        #try({Z.stat[2] <- owstat(trc_data2)},silent=T)
        # rho12 <- sqrt(1/2)
        # sigma <- matrix(c(1,rho12,rho12,1),nrow=2)
        # f <- function(t){
        #   pmvnorm(c(a1,-Inf),c(b1,t),c(0,0),sigma=sigma)+alpha_vec[1]
        # }
        # fb <- function(t){
        #   pmvnorm(c(a1-beta_H1/se1,t-beta_H1/se2),c(b1-beta_H1/se1,Inf),c(0,0),sigma=sigma)+beta_vec[1]
        # }
        # a2 <- uniroot(function(t){f(t)[1]-alpha_vec[2]},c(-5,5))$root
        # b2 <- uniroot(function(t){fb(t)[1]-beta_vec[2]},c(-5,5))$root
        
        dif = sum(trt==1)-sum(trt==0)
        
        if(Z.stat[2] > a2){
          rej <- c(rej,1);early_stop_s <- c(early_stop_s,1)
          early_stop_f <- c(early_stop_f,0); dif01 = c(dif01,dif)
          dif_ms <- c(dif_ms,dif_median_surv(trc_data2));
          #a = survfit(Surv(y,del)~1,trc_data2)
          #dif_ms <- c(dif_ms,a$surv[which(a$time<1)[length(which(a$time<1))]]);
          # res <- summary(survfit(Surv(y,del)~1,trc_data2))
          # dif_ms <- c(dif_ms,rmst(res$time,res$surv,type="l",max.time=1))
          next}else if(Z.stat[2]<b2){rej <- c(rej,0);early_stop_s <- c(early_stop_s,0)
          early_stop_f <- c(early_stop_f,1); dif01 = c(dif01,dif)
          dif_ms <- c(dif_ms,dif_median_surv(trc_data2));
          #a = survfit(Surv(y,del)~1,trc_data2)
          
          #dif_ms <- c(dif_ms,a$surv[which(a$time<1)[length(which(a$time<1))]]);
          # res <- summary(survfit(Surv(y,del)~1,trc_data2))
          # dif_ms <- c(dif_ms,rmst(res$time,res$surv,type="l",max.time=1))
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
              trt = rep(0:1,each=n)
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
            
            data_3 = DataModel(X=XX[,1:p], XParam=XParam,r=sort(r)[(2*n1+2*n2+1):length(r)],
                               G = trt, GParam=GParam, IntParam=IntParam,
                               n=n*2, p=p, baseParam=baseParam, sigma=sigma,
                               R=0.72, D=D, cenParam=cen_rate, mod=mod)
            #data_3 = AFTModel(XX[,1:p], XParam, G=trt, GParam, IntParam, n=n*2, p, sigma, R, D, cenParam=cen_rate)
            tc3 = data_3$tc;data_3 = data_3$data
            t_c3 = 0.72
            
            trc_data3 = data_3
            trc_data3$y = ifelse(tc3>3,data_3$y - (tc3-3),trc_data3$y)
            trc_data3$del = ifelse(tc3>3,0,data_3$del)
            
            trc_data3 <- rbind(trc_data3,
                               data.frame(y=ifelse(tc1>D,data_1$y - (tc1-D),data_1$y),
                                          del = ifelse(tc1>D,0,data_1$del),
                                          X1 = data_1$X1,X2=data_1$X2,G=data_1$G),
                               data.frame(y=ifelse(tc2>D,data_2$y - (tc2-D),data_2$y),
                                          del = ifelse(tc2>D,0,data_2$del),
                                          X1 = data_2$X1,X2=data_2$X2,G=data_2$G))
            
            fit <- coxph(Surv(y,del)~G,data=trc_data3)
            # summary(fit)
            Z.stat[3] <- owstat(trc_data3)
            #Z.stat[3] <-  summary(fit)$coefficients[1,2]
            #try(Z.stat[3] <- owstat(trc_data3),silent=T)
            
            #se3 <- summary(fit)$coefficients[1,3]
            
            # rho13 <- sqrt(1/3); rho23 <- sqrt(2/3)
            # sigma <- matrix(c(1,rho12,rho13,rho12,1,rho23,rho13,rho23,1),nrow=3,byrow=T)
            # #
            # f <- function(t){
            #   pmvnorm(c(a1,a2,-Inf),c(b1,b2,t),c(0,0,0),sigma=sigma)+alpha_vec[2]
            # }
            #a3 <- uniroot(function(t){f(t)[1]-alpha_vec[3]},c(-10,10))$root
            dif_ms <- c(dif_ms,dif_median_surv(trc_data3))
            #a = survfit(Surv(y,del)~1,trc_data3)
            #dif_ms <- c(dif_ms,a$surv[which(a$time<1)[length(which(a$time<1))]])
            # res <- summary(survfit(Surv(y,del)~1,trc_data3))
            # dif_ms <- c(dif_ms,rmst(res$time,res$surv,type="l",max.time=1))
            dif01 <- c(dif01,dif)
            if(Z.stat[3]>a3){rej=c(rej,1)}else{rej=c(rej,0)}
          }
      }
  }
  iter
  out <- rbind(out,data.frame(sc=sc,rej=mean(rej),
                              early_stop_s=mean(early_stop_s),
                              early_stop_f=mean(early_stop_f),
                              dif01=mean(dif01,na.rm=T),
                              dif_ms=mean(dif_ms,na.rm=T)))
  print(out)
}



mean(rej); mean(early_stop_s)
mean(early_stop_f); mean(dif01)
iter