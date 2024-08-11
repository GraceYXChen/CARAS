library(survival)
library(mvtnorm)
library(gsDesign)
library(survminer)
library(rpact)
design <- getDesignGroupSequential(
  # assume one-sided test, equal information increment
  # OF alpha and beta spending function
  sided = 1, alpha = 0.05, beta = 0.2,
  informationRates = c(1/3, 2/3, 1), kMax=3,
  typeOfDesign = "asOF",typeBetaSpending = "bsOF"
)

sampleSizeResult <- getSampleSizeSurvival(
  design = design, lambda2 =1, lambda1 = exp(-0.5),
  dropoutRate1 = .1, dropoutRate2 = .1, dropoutTime = 1,
  followUpTime = 1, # the whole duration is the maximum accrual time 3
  accrualTime = c(0,3)
)

sampleSizeResult <- getSampleSizeSurvival(
  design = design, lambda2 =1, lambda1 = exp(-0.5),
  dropoutRate1 = .1, dropoutRate2 = .1, dropoutTime = 1,
  followUpTime = 1, # the whole duration is the maximum accrual time 3
  accrualTime = c(0,R)
)

kable(sampleSizeResult)  # study design
kable(summary(sampleSizeResult))
  
AFTModel = function(X, XParam, G, GParam, IntParam, n, p, sigma, R, D, cenParam){
  #### simulate Cox model data with Weibull baseline hazard #####
  ## XParam: vector of beta param in Cox model
  ## G: treatment indicator
  ## GParam: parameter of treatment indicator
  ## IntParam: parameter of interaction between G and X
  ## X: n by p matrix
  ## baseParam: a vector of a and lambda in weibull model
  ## R: length of recruitment 
  ## D: overall trial duration
  ## cenParam: exponential censoring rate
  ## return the data ordered by the recruit time
  ######################################################
  W = rnorm(n)
  t = exp(X %*% XParam + G * GParam + apply(X,2,function(u){u*G}) %*% IntParam + sigma*W)

  r = runif(n,0,R) # patient entry calendar time
  cen = rexp(n, cenParam)
  
  y = pmin(cen,t)
  del = ((cen<t) + (r+y>D)==0)
  y = ifelse(r+y>D,D-r,y)
  #y[which(r + y > D)] = D-r
  y = y[order(r,decreasing=F)]
  del = del[order(r,decreasing=F)]
  X = X[order(r,decreasing=F),]
  G = G[order(r,decreasing=F)]
  data = data.frame(cbind(y,del,X,G))
  return(list(data=data,tc=r+y))
}

DataModel = function(X, XParam, r, G, GParam, IntParam, 
                    n, p, baseParam=NULL, R, D, cenParam,
                    mod=c("cox","AFT"), sigma=NULL){
  #### simulate Cox model data with Weibull baseline hazard #####
  ## XParam: vector of beta param in Cox model
  ## G: treatment indicator
  ## GParam: parameter of treatment indicator
  ## IntParam: parameter of interaction between G and X
  ## X: n by p matrix
  ## baseParam: a vector of a and lambda in weibull model
  ## R: length of recruitment 
  ## D: overall trial duration
  ## cenParam: exponential censoring rate
  ## mod: Cox PH model or non PH AFT model 
  ## return the data ordered by the recruit time
  ######################################################
  if(mod=="cox"){
    a = baseParam[1] # shape
    lamb = baseParam[2] # scale 
    U = runif(n)
     t = lamb * (-log(1-U) / as.numeric(exp(
       X %*% XParam + G * GParam + apply(X,2,function(u){u*G}) %*% IntParam)))^(1/a)
    # log logistic basedline hazard (a=1,b=lamb)
    # t = (exp(-log(1-U) / as.numeric(exp(
    #   X %*% XParam + G * GParam + apply(X,2,function(u){u*G}) %*% IntParam)))-1)^(1/lamb)
   }else if(mod=="AFT"){
    W = rnorm(n)
    t = exp(X %*% XParam + G * GParam + apply(X,2,function(u){u*G}) %*% IntParam + sigma*W)
  }
  cen = rexp(n, cenParam)
  
  y = pmin(cen,t) # time to failure
  del = ((cen<t) + (r+y>D)==0) # censor at D or not
  y = ifelse(r+y>D,D-r,y)
  #y[which(r + y > D)] = D-r
  # y = y[order(r,decreasing=F)]
  # del = del[order(r,decreasing=F)]
  # X = X[order(r,decreasing=F),]
  # G = G[order(r,decreasing=F)]
  data = data.frame(cbind(y,del,X,G))
  return(list(data=data,tc=r+y))
}


coxScore <- function(X, y, del, nullCoef=NULL){
  ################## Score vector in Cox model ################
  ## X: n*p covariate matirx
  ## y: observed time
  ## del: censoring indicator
  ## nullCoef: a vector of the index of the null coefficients under H0
  ##           e.g. c(1,2) indicates H0: beta1 = beta2 = 0 
  #############################################################
  if(length(nullCoef)==dim(X)[2]){
    nX = X
    z = X
    e = rep(1,dim(X)[1])
  }else{
    nX = X[,-nullCoef]
    z = as.matrix(X[,nullCoef])
    e = exp(nX %*% matrix(coxph(Surv(y,del)~nX)$coef))
  }
  nd = length(y[del==1])
  u = rep(0,length(nullCoef))
  I = matrix(0,nrow = length(nullCoef),ncol=length(nullCoef))
  for(i in 1:nd){
    ti = sort(y[del==1])[i]
    a = which(y==ti)[1]
    r1 = which(y>=ti)
    u = u + del[a] * (z[a,] - apply(as.matrix((z[r1,]*e[r1])/sum(e[r1])),2,sum))
    w = e[r1] / sum(e[r1])
    for(g in 1:length(nullCoef)){
      for(l in 1:length(nullCoef)){
        I[g,l] = I[g,l] + sum(z[r1,g] * z[r1,l] * w) - sum(z[r1,g]*w) * sum(z[r1,l]*w)
      }
    }
  }
  return(list(score=u, var=I))
}


Acoef_RAR <- function(data, tau){
  ############### Compute A in Var(S(t)) ############
  ### data: a dataframe contains y,del
  ########################################################
  fit <- coxph(Surv(y,del) ~., data=data)
  # baseline survival function
  obj <- survfit(fit,newdata=data.frame(X1=0,X2=0))
  S0_tau <- summary(obj,time=tau,extend=T)$surv
  exp_linear_pred <- exp(predict(fit, type="lp"))  # exp(b'newx)
  exp_linear <- exp(predict(fit, type="lp"))   # exp(b'x)
  S_tau <- S0_tau^(exp_linear_pred)   
  
  km_fit <- survfit(Surv(y, del) ~ 1, data)
  km_fit <- summary(km_fit)
  if(sum(km_fit$time<=tau)==0){return(S_tau^2)}
  else{
  l_time <- km_fit$time[1:which(km_fit$time<=tau)[length(which(km_fit$time<=tau))]]  # distinct failure times before tau

  prod_term = 0
  for(t in l_time){
    w = sum(exp_linear[data$y>=t])
    pi = exp(-exp_linear_pred / w)
    qi = 1 - pi
    S_ti = km_fit$surv[which(km_fit$time==t)]
    prod_term = prod_term + ifelse(S_ti==0,0,qi / (S_ti * pi))
  }
  # t_lplus1 <- km_fit$time[which(km_fit$time>tau)[1]]
  # 
  # if(is.na(t_lplus1)==F & t!=t_lplus1 & tau>t){
  # 
  # p_lplus1 = exp(-exp_linear_pred * (tau-t)/(t_lplus1-t) / 
  #                  sum(exp_linear[data$y >= t_lplus1]))
  # prod_term = prod_term + ifelse(km_fit$surv[which(km_fit$time==t_lplus1)]!=0,(1-p_lplus1) / 
  #                                  (km_fit$surv[which(km_fit$time==t_lplus1)]*p_lplus1),0)
  # }
  A = S_tau[1]^2 * prod_term[1]
  return(A)
  }
}

Acoef <- function(data, newx, tau){
  ############### Compute A in Var(S(t|newx)) ############
  ### data: a dataframe contains y,del,X (includes X1 X2 and X3)
  ### newx: a data frame with X1, X2, X3
  ########################################################
  fit <- coxph(Surv(y,del) ~., data=data)
  # baseline survival function
  obj <- survfit(fit,newdata=data.frame(X1=0,X2=0))
  S0_tau <- summary(obj,time=tau,extend=T)$surv
  exp_linear_pred <- exp(predict(fit, type="lp",newdata = newx))  # exp(b'newx)
  exp_linear <- exp(predict(fit, type="lp"))   # exp(b'x)
  S_tau <- S0_tau^(exp_linear_pred)   # survival probability at tau given newx
  
  km_fit <- survfit(Surv(y, del) ~ 1, data)
  km_fit <- summary(km_fit)
  if(sum(km_fit$time<=tau)==0){return(S_tau^2)}
  else{
  z = 0
  dataX = data[,-which(colnames(data)%in%c("y","del"))]
  
  prod_term = 0
  l_time <- km_fit$time[1:which(km_fit$time<=tau)[length(which(km_fit$time<=tau))]]  # distinct failure times before tau
  for(t in l_time){
    term_1 = apply((exp_linear * dataX)[which(data$y>=t),],2,sum)
    w = sum(exp_linear[data$y>=t])
    z = z + term_1 / w^2 - as.matrix(newx) / w
    pi = exp(-exp_linear_pred / w)
    qi = 1 - pi
    #S_ti = km_fit$surv[which(km_fit$time==t)]
    S_ti = summary(obj,time=t,extend=T)$surv^(exp_linear_pred)
    if(pi>0){
    prod_term = prod_term + ifelse(S_ti==0,0,qi / (S_ti * pi))
    }
  }
  z = as.vector(z)

  v = 0
  for(t in data$y){
    w = sum(exp_linear[data$y>=t])
    ind = which(data$y>=t)
    v = v + (t(as.matrix(dataX[ind,]))%*%(as.matrix(dataX[ind,])*exp_linear[ind])/w -
      t(as.matrix(dataX[ind,]*exp_linear[ind])) %*% (as.matrix(dataX[ind,]*exp_linear[ind]))/ w^2)
  }

  inv_vbar = solve(v / length(data$y))
  
  eps1 = mean(data$del)
  # 
  # prod_term = 0
  # for(t in l_time){
  #   w = sum(exp_linear[data$y>=t])
  #   pi = exp(-exp_linear_pred / w)
  #   qi = 1 - pi
  #   S_ti = km_fit$surv[which(km_fit$time==t)]
  #   prod_term = prod_term + qi / (S_ti * pi)
  # }
# 
  # try comment the following and improve the reuslt
  # t_lplus1 <- km_fit$time[which(km_fit$time>tau)[1]]
  # if(is.na(t_lplus1)==F & t!=t_lplus1 & tau>t){
  # 
  #   t <- km_fit$time[which(km_fit$time>tau)[1]-1]
  #   p_lplus1 = exp(-exp_linear_pred * (tau-t)/(t_lplus1-t) /
  #                  sum(exp_linear[data$y >= t_lplus1]))
  #    prod_term = prod_term + ifelse(km_fit$surv[which(km_fit$time==t_lplus1)]!=0 & p_lplus1!=0,(1-p_lplus1) /
  #       (km_fit$surv[which(km_fit$time==t_lplus1)]*p_lplus1),0)
  # }

  A = S_tau^2 * (prod_term +
                   exp_linear_pred^2 * z %*% inv_vbar %*% z / eps1)
  # A = ifelse(S_tau==1, 1, A = S_tau^2 * (prod_term +
  #   exp_linear_pred^2 * z %*% inv_vbar %*% z / eps1))
  #A = ifelse(A==0,0.01,A)
  return(A)
  }
}

predict_surv <- function(tau, newx, data){
  ################################################
  # tau: analysis time
  # newx: new covariates of the predicted data (data frame, same as data)
  # data: a data frame containing column 'y', 'del', 'X1','X2',...'Xp'
  ################################################
  fit <- coxph(Surv(y,del) ~., data=data)
  obj <- survfit(coxph(Surv(y,del) ~1, data=data),newdata = newx)
  S0_tau <- summary(obj,time=tau)$surv
  S <- S0_tau^(exp(predict(fit, type="lp",newdata = newx)))  
  return(S) 
}




allocationRatio <- function(type=c("Neyman","Survival","Hazard"),
                            A1, A2, data1, data2, newx){
  ################### allocation ratio ###############
  ## A1, A2: The denominator of the variance of survival function
  ####################################################
  if(type=="Neyman"){
    return( (sqrt(A1)) / (sqrt(A1)+sqrt(A2)) )
  }else if(type=="Survival"){
    fit1 <- coxph(Surv(y,del) ~., data=data1)
    fit2 <- coxph(Surv(y,del) ~., data=data2)
    obj1 <- survfit(coxph(Surv(y,del) ~1, data=data1),newdata = newx)
    obj2 <- survfit(coxph(Surv(y,del) ~1, data=data2),newdata = newx)
    S01_tau <- summary(obj1,time=tau)$surv
    S02_tau <- summary(obj2,time=tau)$surv
    S1 <- S01_tau^(exp(predict(fit1, type="lp",newdata = newx)))  
    S2 <- S02_tau^(exp(predict(fit2, type="lp",newdata = newx))) 
    return((sqrt((1-S2)*A1)) / 
             (sqrt((1-S2)*A1)+sqrt(A2*(1-S1))))
  }else{
    fit1 <- coxph(Surv(y,del) ~., data=data1)
    fit2 <- coxph(Surv(y,del) ~., data=data2)
    obj1 <- survfit(coxph(Surv(y,del) ~1, data=data1),newdata = newx)
    obj2 <- survfit(coxph(Surv(y,del) ~1, data=data2),newdata = newx)
    S01_tau <- summary(obj1,time=tau)$surv
    S02_tau <- summary(obj2,time=tau)$surv
    S1 <- S01_tau^(exp(predict(fit1, type="lp",newdata = newx)))  
    S2 <- S02_tau^(exp(predict(fit2, type="lp",newdata = newx)))
    H1 <- -log(S1)
    H2 <- -log(S2)
    return((sqrt(A1*H1*H2)) / 
             (sqrt(A1*H1*H2)+H1*sqrt(A2)))
  }
}

dif_median_surv <- function(data){
  ####### data contains y, del columns
  fit <- surv_fit(Surv(y,del) ~ 1, data)
  r = surv_median(fit, combine = FALSE)
  return(r$median)
}

oneyrsurv <- function(data){
  fit <- coxph(Surv(y,del) ~., data=data)
  # baseline survival function
  obj <- survfit(fit)
  S0_tau <- summary(obj,time=1,extend=T)$surv
  exp_linear_pred <- exp(predict(fit, type="lp"))  # exp(b'newx)
  exp_linear <- exp(predict(fit, type="lp"))   # exp(b'x)
  S_tau <- S0_tau^(exp_linear_pred)   
  return(S_tau)
}



owstat <- function(trc_data1){
  time = trc_data1$y
  group = trc_data1$G
  event = trc_data1$del
  
  n <- length(time)
  ng <- table(group)
  group <- factor(group)
  mod = glm(G~X1+X2,trc_data1,family="binomial")
  wi = as.numeric(predict(mod,type="response"))
  wi[trc_data1$G==1] = 1-wi[trc_data1$G==1]
  
  Ag <- aggregate(event, by = list(time = time, group = group), 
                  FUN = sum, drop = FALSE)
  Ag$x <- ifelse(is.na(Ag$x), 0, Ag$x)
  tab <- data.frame(time = Ag$time[Ag$group == levels(group)[1]], 
                    event1 = Ag$x[Ag$group == levels(group)[1]], 
                    event2 = Ag$x[Ag$group == levels(group)[2]])
  Agw = aggregate(event*wi, by = list(time = time, group = group), 
                      FUN = sum, drop = FALSE)
  Agw$x <- ifelse(is.na(Agw$x), 0, Agw$x)
  tab$event1w <- Agw$x[Agw$group==levels(group)[1]]
  tab$event2w <- Agw$x[Agw$group==levels(group)[2]]
  
  w <- aggregate(wi, by = list(time = time, group = group), 
            FUN = sort, drop = FALSE)

  if(class(w$x)=="list"){
    w$x = lapply(w$x,function(u){ifelse(is.na(u),0,u)})
    w$x <- unlist(lapply(w$x,function(u){u[1]}))
  }
  w$x <- ifelse(is.na(w$x),0,w$x)
  tab$w1 <- w$x[w$group==levels(group)[1]]
  tab$w2 <- w$x[w$group==levels(group)[2]]
  
  tab$atrisk1 <- NA; tab$atrisk1w <- NA
  tab$atrisk2 <- NA; tab$atrisk2w <- NA
  for (i in 1:dim(tab)[1]) {
    tab$atrisk1[i] <- sum(time[group == levels(group)[1]] >= 
                            tab$time[i])
    tab$atrisk2[i] <- sum(time[group == levels(group)[2]] >= 
                            tab$time[i])
    tab$atrisk1w[i] <- sum(wi[group == levels(group)[1]][time[group == levels(group)[1]] >= tab$time[i]])
    tab$atrisk2w[i] <-sum(wi[group == levels(group)[2]][time[group == levels(group)[2]] >= tab$time[i]])
  }
  
  nz <- dim(tab)[1]
  tab$atrisk <- tab$atrisk1 + tab$atrisk2
  tab$event <- tab$event1 + tab$event2
  tab$eventw <- tab$event1w + tab$event2w
  tab$atriskw <- tab$atrisk1w + tab$atrisk2w
  
  D <- tab[tab$event > 0, ]
  
  Gw <- sum(D$event1w - D$atrisk1w * D$eventw / D$atriskw)
  varGw <- 0
  for(j in 1:(dim(D)[1]-1)){
    varGw <- varGw + (D$event*(D$atrisk-D$event)/D$atrisk/(D$atrisk-1))[j] * 
      (((D$atrisk1w/D$atriskw)^2)[j]*sum((tab$w1^2)[1:(D$atrisk[j])],na.rm=T) + ((D$atrisk2w/D$atriskw)^2)[j]*sum((tab$w2^2)[1:(D$atrisk[j])],na.rm=T))
  }
  return(Gw/sqrt(varGw))
}
owstat(trc_data1)


Z.stat=NULL
for(i in 1:1000){
  n = n1;  X1 = rbinom(n*2,1,p1); X2 = rbinom(n*2,1,p2);
  #G = rep(1,2*n); G[sample(1:(2*n),n,replace=F)]=0
  G = rep(0:1,n)
  XX = cbind(X1,X2,G)
  data_1 = DataModel(X=XX[,1:p], XParam=XParam,
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
  Z.stat[i] <- owstat(trc_data1)
}
hist(Z.stat)
hist(rnorm(1000),add=T,col="red")




owstat_mok <- function(trc_data1){
  time = trc_data1$y
  group = trc_data1$G
  event = trc_data1$del
  
  n <- length(time)
  ng <- table(group)
  group <- factor(group)
  mod = glm(G~Z,trc_data1,family="binomial")
  wi = as.numeric(predict(mod,type="response"))
  wi[trc_data1$G==1] = 1-wi[trc_data1$G==1]
  
  Ag <- aggregate(event, by = list(time = time, group = group), 
                  FUN = sum, drop = FALSE)
  Ag$x <- ifelse(is.na(Ag$x), 0, Ag$x)
  tab <- data.frame(time = Ag$time[Ag$group == levels(group)[1]], 
                    event1 = Ag$x[Ag$group == levels(group)[1]], 
                    event2 = Ag$x[Ag$group == levels(group)[2]])
  Agw = aggregate(event*wi, by = list(time = time, group = group), 
                  FUN = sum, drop = FALSE)
  Agw$x <- ifelse(is.na(Agw$x), 0, Agw$x)
  tab$event1w <- Agw$x[Agw$group==levels(group)[1]]
  tab$event2w <- Agw$x[Agw$group==levels(group)[2]]
  
  w <- aggregate(wi, by = list(time = time, group = group), 
                 FUN = sort, drop = FALSE)
  
  if(class(w$x)=="list"){
    w$x = lapply(w$x,function(u){ifelse(is.na(u),0,u)})
    w$x <- unlist(lapply(w$x,function(u){u[1]}))
  }
  w$x <- ifelse(is.na(w$x),0,w$x)
  tab$w1 <- w$x[w$group==levels(group)[1]]
  tab$w2 <- w$x[w$group==levels(group)[2]]
  
  tab$atrisk1 <- NA; tab$atrisk1w <- NA
  tab$atrisk2 <- NA; tab$atrisk2w <- NA
  for (i in 1:dim(tab)[1]) {
    tab$atrisk1[i] <- sum(time[group == levels(group)[1]] >= 
                            tab$time[i])
    tab$atrisk2[i] <- sum(time[group == levels(group)[2]] >= 
                            tab$time[i])
    tab$atrisk1w[i] <- sum(wi[group == levels(group)[1]][time[group == levels(group)[1]] >= tab$time[i]])
    tab$atrisk2w[i] <-sum(wi[group == levels(group)[2]][time[group == levels(group)[2]] >= tab$time[i]])
  }
  
  nz <- dim(tab)[1]
  tab$atrisk <- tab$atrisk1 + tab$atrisk2
  tab$event <- tab$event1 + tab$event2
  tab$eventw <- tab$event1w + tab$event2w
  tab$atriskw <- tab$atrisk1w + tab$atrisk2w
  
  D <- tab[tab$event > 0, ]
  
  Gw <- sum(D$event1w - D$atrisk1w * D$eventw / D$atriskw)
  varGw <- 0
  for(j in 1:(dim(D)[1]-1)){
    varGw <- varGw + (D$event*(D$atrisk-D$event)/D$atrisk/(D$atrisk-1))[j] * 
      (((D$atrisk1w/D$atriskw)^2)[j]*sum((tab$w1^2)[1:(D$atrisk[j])],na.rm=T) + ((D$atrisk2w/D$atriskw)^2)[j]*sum((tab$w2^2)[1:(D$atrisk[j])],na.rm=T))
  }
  return(Gw/sqrt(varGw))
}
