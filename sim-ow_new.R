D = 4
ana_time = c(1.7,2.7,D)
a1 = 3.2; a2 = 2.141; a3 = 1.695
b1 = -0.415; b2 = 0.917
n1 = 31; n2 = 18; n3 = 5
p1 = .5; p2 = .5
baseParam=c(1,1); cen_rate= .1; p = 2; R = 1 

method="CARA1"
GParam_all = c(rep(0,4),rep(-0.6,5),rep(0,3),rep(-0.6,6))
XParam_all = matrix(c(0,0, -1,0, -2,0, -5,0, 
                      0,0, -1,0, -2,0, 0,0, 
                      -2,0, -2,-1, -2,-2,
                      -2,-5, -1,-1,  -1,-2,
                      0,0, -1,-2),byrow = T,ncol=2)
IntParam_all = matrix(c(rep(0,14),-2,0,-2,0,
                        rep(0,10),rep(-2,4)),byrow=T,ncol=2)

sc=5
set.seed(2023)
out <- NULL
for(sc in c(1:16)[-c(1,5,9)]){
  GParam = GParam_all[sc]
  XParam = XParam_all[sc,]
  IntParam = IntParam_all[sc,]
  rej = NULL; early_stop_f <- NULL;dif01 <- NULL;early_stop_s <- NULL;dif_ms <- NULL;
  ss = NULL; H = NULL; ms <- NULL; fail=NULL
  for(iter in 1:1000){
    try({
    # print(c(mean(rej),mean(early_stop_f),mean(early_stop_s),
    #                  mean(dif01,na.rm=T),mean(ms,na.rm=T)))
    n = n1;  X1 = rbinom(n*2,1,p1); X2 = rbinom(n*2,1,p2);
    G = rep(0:1,each=n)
    XX = cbind(X1,X2,G)
    
    # r = pmin(rexp((n1+n2+n3)*2,1),3)
    r = runif((n1+n2+n3)*2,0,3)
    # r = rbeta((n1+n2+n3)*2,0.5,0.5)*3   # Accrual time
    
    r1 <- r[which(r<=sort(r)[(n1*2)])]
    r2 <- r[which(r<=sort(r)[((n1+n2)*2)] & r>sort(r)[(n1*2)])]
    r3 <- r[which(r>sort(r)[((n1+n2)*2)])]
    
    data_1 = DataModel(X=XX[,1:p], XParam=XParam,r=r1,
                       G = G, GParam = GParam, IntParam = IntParam,
                       n=n*2, p=p, baseParam=baseParam, sigma=sigma,
                       R=1.74, D=D, cenParam=cen_rate,mod=mod)
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
      #dif_ms <- c(dif_ms,RISCA:::rmst(trc_data1[trc_data1$G==1,]$y,trc_data1[trc_data1$G==1,]$del,0.5,"s"))
      #dif_ms <- c(dif_ms, mean((data_1$data$del==1 & data_1$tc>=1)))
      #dif_ms <- c(dif_ms,dif_median_surv(data_1$data)); 
      
      ss <- c(ss, n1*2)
      H <- c(H,exp(-mean(trc_data1$del==0 & trc_data1$y+trc_data1$r>=ana_time[1])))
      ms <- c(ms,NA)
      dif_ms <- c(dif_ms, NA)
      fail <- c(fail,NA)
      next}else if(Z.stat[1]<b1){
        rej <- c(rej,0);early_stop_s <- c(early_stop_s,0)
        early_stop_f <- c(early_stop_f,1); dif01 = c(dif01,0)
        #dif_ms <- c(dif_ms, mean((data_1$data$del==1 & data_1$tc>=1)))
        #dif_ms <- c(dif_ms,RISCA:::rmst(trc_data1[trc_data1$G==1,]$y,trc_data1[trc_data1$G==1,]$del,0.5,"s"))
        #dif_ms <- c(dif_ms,dif_median_surv(data_1$data))
        ss <- c(ss, n1*2)
        H <- c(H,exp(-mean(trc_data1$del==0 & trc_data1$y+trc_data1$r>=ana_time[1])))
        ms <- c(ms,NA)
        dif_ms <- c(dif_ms, NA)
        fail <- c(fail,NA)
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
                         n=n*2, p=p, baseParam=baseParam, mod = mod,
                         sigma=sigma, R=1, D=D, cenParam=cen_rate)
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
        fail = c(fail, sum(trc_data2$del==1))

        dif_ms <- c(dif_ms,mean(data_2only$y<2 & data_2only$del==1,na.rm=T))

 
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
          fail = c(fail, sum(trc_data2$del==1))
       
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
                               R=1.26, D=D, cenParam=cen_rate, mod=mod)      
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
            fail = c(fail, NA)
            if(Z.stat[3]>a3){rej=c(rej,1)}else{rej=c(rej,0)}
          }
    }},silent=T)
  }

  out <- rbind(out,data.frame(sc=sc,rej=mean(rej),
                              early_stop_s=mean(early_stop_s),
                              early_stop_f=mean(early_stop_f),
                              dif01=mean(dif01,na.rm=T),
                              deathr=mean(dif_ms,na.rm=T),
                              H=mean(H,na.rm=T),
                              ms = mean(ms,na.rm=T),
                              fail=mean(fail,na.rm=T)))
  print(out)
}
out <- out[order(out$sc),]

write.csv(out,"D://CARAOW207.csv")
