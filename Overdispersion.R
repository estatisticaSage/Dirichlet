library(dplyr)
library(ggplot2)
library(DirichletReg)
library(tidyr)
library(ggrepel)
library(rstan)
#### FUNCTIONS ####
s<-function(alpha,y){
  alpha_0<- sum(alpha)
  digamma(alpha_0)-digamma(alpha)+log(y)
}


# overdispersion em cada observacao
D_ij<-function(alpha,y){
  alpha_0<- sum(alpha)
  n<-length(y)
  Dij<-trigamma(alpha_0)-trigamma(alpha)+s(alpha,y)^2
  return(Dij)
}


D.Aitchison <- function(y1,y2){
  r <- y1/y2
  d <-  sqrt(sum((log(r) - mean(log(r)))^2))
  return(d)
}

g<-function(y){
  C<-dim(y)[2]
  p<-NA
  n<-dim(y)[1]
  for(i in 1:C){
    p[i]<-(prod(y[,i]))^(1/n)
  }
  return(p)
}


#### OVERDISPERSION ####

#CALCULO NOS DADOS SIMULADOS #
load("~/Sim_r.rda")
theta_estimado<-Par.Estimates[1,]$mean

betas_estimados<- matrix(c(Par.Estimates[2:4,]$mean,0,Par.Estimates[5:7,]$mean,0),ncol=2,nrow=4)

E.A<-exp(betas_estimados[,1])*theta_estimado
E.B<-exp(betas_estimados[,1] + betas_estimados[,2])*theta_estimado


Dij.A<-as.data.frame(D_ij(E.A,Y[1:250,]))
Dij.B<-as.data.frame(D_ij(E.B,Y[251:500,]))


Dij.A%>%str
apply(Dij.A,2,sum)
apply(Dij.A,1,sum) # soma as linhas para ver qual individuo tem mais overdispersion

Dij.A%>%dplyr::select(V3)%>%ggplot()+
  geom_density(aes(V3))

dado <- Dij.A%>%dplyr::select(V3)%>%
  mutate(Id = 1:n())%>%  
  dplyr::filter(V3 > quantile(Dij.A$V3,0.90))

dado1 <- Dij.A%>%dplyr::select(V3)%>%
  mutate(Id = 1:n())%>%  
  dplyr::filter(V3 > quantile(Dij.A$V3,0.975))

dado$id.text <- paste(dado$Id,sep = "")

p <- dado%>%
  ggplot()+
  geom_segment(aes(x=Id,xend=Id,y=0,yend=V3),size=1.5,alpha=0.5)+
  geom_text_repel(data = dado,aes(x=Id, y=V3, 
                                  label=id.text))+
  geom_point(aes(x=Id,y=V3),size=2.5)+
  geom_point(data=dado1,aes(x=Id,y=V3),size=2.5,color="red")+
  geom_segment(data=dado1,aes(x=Id,xend=Id,y=0,yend=V3),size=1.5,alpha=0.5,color="red")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,angle=0,hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=20),
        strip.text =  element_text(size=18))+
  labs(x="Possible influential observations",
       y="Overdispersion",subtitle = "A")
p


### Sites ####

Dij.A <-as.data.frame(D_ij(E.A,Y[1:250,]))
Dij.B <-as.data.frame(D_ij(E.A,Y[251:500,]))


Dij.A<-Dij.A%>%mutate(Site = "A")
Dij.B<-Dij.B%>%mutate(Site = "B")



Dij.A<-Dij.A%>%dplyr::mutate(Id=1:n())
Dij.B<-Dij.B%>%dplyr::mutate(Id=1:n())

Dij.A$id.text <- paste(Dij.A$Id,sep = "")
Dij.B$id.text<-paste(Dij.B$Id,sep='')

Dij.Dado<-rbind(Dij.A,Dij.B)

dado <- Dij.Dado%>%dplyr::select(V3,Site,Id,id.text)%>%
  dplyr::filter(V3 > quantile(Dij.Dado$V3,0.95))

dado1 <- Dij.Dado%>%dplyr::select(V3,Site,Id,id.text)%>%
  dplyr::filter(V3 > quantile(Dij.Dado$V3,0.975))


#ajeitei o grafico para saur sempre de 0 a 250 obs
p <- dado%>%
  ggplot()+
  geom_segment(aes(x=Id,xend=Id,y=0,yend=V3),size=1.5,alpha=0.5)+
  geom_text_repel(data = dado,aes(x=Id, y=V3, 
                                  label=id.text))+
  geom_point(aes(x=Id,y=V3),size=2.5)+
  geom_point(data=dado1,aes(x=Id,y=V3),size=2.5,color="red")+
  geom_segment(data=dado1,aes(x=Id,xend=Id,y=0,yend=V3),size=1.5,alpha=0.5,color="red")+
  facet_wrap(~Site)+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,angle=0,hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=20),
        strip.text =  element_text(size=18))+
  labs(x="Possible influential observations",
       y="Overdispersion")
p



# CALCULO NO RESULTADO DA NOSSA REGRESSAO #
library(rstan)
library(DirichletReg)
library(dplyr)
load("C:/Users/Mariana/Downloads/BayesDirReg (1).rda")
fit1

betas<-extract(fit1)$beta
betas
#extraindo os betas
theta<-extract(fit1)$exptheta

# media a posteriori deles
betas.E<-apply(betas[,,1],MAR=2,FUN=mean)
theta.E<-mean(theta)

for(i in 2:6){
  x<-apply(betas[,,i],MAR=2,FUN=mean)
  betas.E<-rbind(betas.E,x)
}
B.E<-as.matrix(betas.E)
B.E # guardados as medias a posteriori dos betas


# calculando os alphas estimados
A.E<-matrix(NA,nrow=9,ncol=6)


#PENSAR NUMA FORMA MAIS PRATICA DE FAZER ISSO 
A.E[1,]<-exp(B.E[,1])*theta.E
A.E[2,]<-exp(B.E[,1]+B.E[,2])*theta.E
A.E[3,]<-exp(B.E[,1]+B.E[,2]+B.E[,3])*theta.E
A.E[4,]<-exp(B.E[,1]+B.E[,2]+B.E[,3]+B.E[,4])*theta.E
A.E[5,]<-exp(B.E[,1]+B.E[,2]+B.E[,3]+B.E[,4]+B.E[,5])*theta.E
A.E[6,]<-exp(B.E[,1]+B.E[,2]+B.E[,3]+B.E[,4]+B.E[,5]+B.E[,6])*theta.E
A.E[7,]<-exp(B.E[,1]+B.E[,2]+B.E[,3]+B.E[,4]+B.E[,5]+B.E[,6]+B.E[,7])*theta.E
A.E[8,]<-exp(B.E[,1]+B.E[,2]+B.E[,3]+B.E[,4]+B.E[,5]+B.E[,6]+B.E[,7]+B.E[,8])*theta.E
A.E[9,]<-exp(B.E[,1]+B.E[,2]+B.E[,3]+B.E[,4]+B.E[,5]+B.E[,6]+B.E[,7]+B.E[,8]+B.E[,9])*theta.E



load("~/IC/Dataset.rda")
options(scipen=999)
RD_SR[,7:15] <- RD_SR[,7:15]/100


#Overdispersion sitio PAB2 ficou bem estranho
dados<-RD_SR%>%filter(Site=='PAB2')
dados<-as.matrix(DR_data(dados[,7:15]))
str(dados)

(D_ij(A.E[,1],dados))
#valores muito altoos muito louco

# Overdispersion sitio PAB3
dados<-RD_SR%>%filter(Site=='PAB3')
dados<-DR_data(dados[,7:15])


D_ij(A.E[,2],dados)

#Overdispersion sitio TIM
dados<-RD_SR%>%filter(Site=='TIM')
dados<-DR_data(dados[,7:15])

D_ij(A.E[,6],dados)


#### R2 MEASURE based  on sum of squares ####
center.Y<-g(Y)

# erro por causa do 0 
CSST<-sum((D.Aitchison(Y,center.Y))^2)

predito<-matrix(as.numeric(Predictive$mean),ncol=4,nrow=500)

CSSE<-sum((D.Aitchison(Y,predito))^2)

R2<- 1 - (CSSE/CSST)

# tudo 0 PORQUE OS DADOS TEM 0 , muitos problemas com os 0
RD_SR[,7:15]<-RD_SR[,7:15]/100
pab2<-RD_SR%>%dplyr::filter(Site=='PAB2')
g(DR_data(pab2[,7:15]))

# nao tem como fazer o CSSE PORQUE NAO FIZEMOS O Y PREDITO

#### R2 BASED ON TOTAL VARIABILITY ####

#install.packages("devtools")
library(devtools)
#install_github("pchiroque/benthic")
library(benthic)

# fazer log ratio usando o pacote da pamela :)

V.L<-var(logratio(Y))
V.L.E<-vat(logratio(predito))
# para fazer a variabilidade total precisamos somar o traco da matriz

tot.var<- sum(diag(V.L))


tot.var.E<-sum(diag(V.L.E))

# parece estar bem ruimm
R2<-tot.var.E/tot.var

V.L<-var(logratio(DR_data(RD_SR[,7:15])))
tot.var<-sum(diag(V.L))
