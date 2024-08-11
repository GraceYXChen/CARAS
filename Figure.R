# Figure: 800*300
d =readxl::read_xlsx("D://fig.xlsx",sheet=1)
d$rej[is.na(d$rej)]=0
library(ggplot2)
library(RColorBrewer)
color = brewer.pal(8,"Paired")[1:8]
color[4] = "purple"

###### rejection prob ######
d =readxl::read_xlsx("D://fig.xlsx",sheet=1)
color = c(colorspace::sequential_hcl(n = 7, palette = "Purples2")[1:2],
          colorspace::sequential_hcl(n = 7, palette = "Reds")[3],
          colorspace::sequential_hcl(n = 7, palette = "Peach")[1],
          colorspace::sequential_hcl(n = 7, palette = "Terrain")[1],
          colorspace::sequential_hcl(n = 7, palette = "Reds")[4],
          colorspace::sequential_hcl(n = 7, palette = "Peach")[2],
          colorspace::sequential_hcl(n = 7, palette = "Terrain")[2])
par(mfrow=c(1,3))
plot(d$rej[d$method=="Trad"],type="o", lty=1,pch=1, ylim=c(0,0.15),cex=0.6,
     cex.axis=0.7,cex.lab=0.7,main="Scenario 1",cex.main=0.8,
     ylab="Rejection Probability",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$rej[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$rej[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$rej[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")
leg = unique(d$method)
leg[3] = "CARAS1"
leg[4] = "CARAS2"
leg[5] = "CARAS3"
leg[6] = "CARASOW1"
leg[7] = "CARASOW2"
leg[8] = "CARASOW3"
legend(c(0.1,0.15),leg,lty=c(2,2,1,3,4,1,3,4),box.lty = 0,cex=0.8,box.lwd = 1,
       col=c("black","blue","red","orange","magenta","red","orange","magenta"),pch=1:8)


d =readxl::read_xlsx("d://fig.xlsx",sheet=2)
plot(d$rej[d$method=="Trad"],type="o", lty=1,pch=1, ylim=c(0.4,0.9),cex=0.6,
     cex.axis=0.7,cex.lab=0.7,main="Scenario 5",cex.main=0.8,
     ylab="Rejection Probability",xlab="Model")
lines(d$rej[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$rej[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$rej[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")




d =readxl::read_xlsx("d://fig.xlsx",sheet=3)
plot(d$rej[d$method=="Trad"],type="o", lty=1,pch=1, ylim=c(0.5,1),cex=0.6,
     cex.axis=0.7,cex.lab=0.7,main="Scenario 9",cex.main=0.8,
     ylab="Rejection Probability",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$rej[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$rej[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$rej[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")



d =readxl::read_xlsx("d://fig.xlsx",sheet=4)
plot(d$rej[d$method=="Trad"],type="o", lty=1,pch=1, ylim=c(0,0.2),cex=0.6,
     cex.axis=0.7,cex.lab=0.7,main="Scenario 3",cex.main=0.8,
     ylab="Rejection Probability",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$rej[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$rej[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$rej[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$rej[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$rej[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")




############# Median Survival Time #################
d =readxl::read_xlsx("d://fig.xlsx",sheet=1)
par(mfrow=c(1,3))
plot(d$MS[d$method=="Trad"],type="o", lty=1,pch=1, cex=0.6,
     cex.axis=0.7,cex.lab=0.7,main="Scenario 1",cex.main=0.8,ylim=c(0.4,1),
     ylab="Median Survival Time",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$MS[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$MS[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$MS[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")


d =readxl::read_xlsx("d://fig.xlsx",sheet=2)
plot(d$MS[d$method=="Trad"],type="o", lty=1,pch=1,cex=0.6,ylim=c(0.4,1.1),
     cex.axis=0.7,cex.lab=0.7,main="Scenario 5",cex.main=0.8,
     ylab="Median Survival Time",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$MS[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$MS[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$MS[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")




d =readxl::read_xlsx("d://fig.xlsx",sheet=3)
plot(d$MS[d$method=="Trad"],type="o", lty=1,pch=1,cex=0.6,ylim=c(0.2,1),
     cex.axis=0.7,cex.lab=0.7,main="Scenario 9",cex.main=0.8,
     ylab="Median Survival Time",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$MS[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$MS[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$MS[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")


d =readxl::read_xlsx("d://fig.xlsx",sheet=4)
plot(d$MS[d$method=="Trad"],type="o", lty=1,pch=1,cex=0.6,ylim=c(0,1),
     cex.axis=0.7,cex.lab=0.7,main="Scenario 3",cex.main=0.8,
     ylab="Median Survival Time",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$MS[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$MS[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$MS[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$MS[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$MS[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")



############# Difference #################
d =readxl::read_xlsx("d://fig.xlsx",sheet=1)
par(mfrow=c(1,3))
plot(d$dif[d$method=="Trad"],type="o", lty=1,pch=1, cex=0.6,
     cex.axis=0.7,cex.lab=0.7,main="Scenario 1",cex.main=0.8,ylim=c(-1,6),
     ylab="Difference of Number of Patients in Control and Experiment",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$dif[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$dif[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="Red")
lines(d$dif[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$dif[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="red")
lines(d$dif[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")


d =readxl::read_xlsx("d://fig.xlsx",sheet=2)
plot(d$dif[d$method=="Trad"],type="o", lty=1,pch=1,cex=0.6,ylim=c(-1,20),
     cex.axis=0.7,cex.lab=0.7,main="Scenario 5",cex.main=0.8,
     ylab="Difference of Number of Patients in Control and Experiment",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$dif[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$dif[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$dif[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$dif[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="Red")
lines(d$dif[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")




d =readxl::read_xlsx("d://fig.xlsx",sheet=3)
plot(d$dif[d$method=="Trad"],type="o", lty=1,pch=1,cex=0.6,ylim=c(-1,20),
     cex.axis=0.7,cex.lab=0.7,main="Scenario 9",cex.main=0.8,
     ylab="Difference of Number of Patients in Control and Experiment",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$dif[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$dif[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$dif[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$dif[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="Red")
lines(d$dif[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")


d =readxl::read_xlsx("d://fig.xlsx",sheet=4)
plot(d$dif[d$method=="Trad"],type="o", lty=1,pch=1,cex=0.6,ylim=c(-1,5),
     cex.axis=0.7,cex.lab=0.7,main="Scenario 3",cex.main=0.8,
     ylab="Difference of Number of Patients in Control and Experiment",xlab="Model",mgp=c(1.5,0.5,0))
lines(d$dif[d$method=="RAR"],type="o",lty=2,pch=2,cex=0.6,col="blue")
lines(d$dif[d$method=="CARA1"],type="o",lty=1,pch=3,cex=0.6,col="red")
lines(d$dif[d$method=="CARA2"],type="o",lty=1,pch=4,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3"],type="o",pch=5,lty=1,cex=0.6,col="magenta")
lines(d$dif[d$method=="CARA1*"],type="o",pch=6,lty=2,cex=0.6,col="Red")
lines(d$dif[d$method=="CARA2*"],type="o",pch=7,lty=2,cex=0.6,col="orange")
lines(d$dif[d$method=="CARA3*"],type="o",pch=8,lty=2,cex=0.6,col="magenta")






