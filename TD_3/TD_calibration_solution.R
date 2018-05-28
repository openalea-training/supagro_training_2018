###__________________________________________________________________________________________##
###_______________________ECUE2:DEMARCHE DE MODELISATION_________________________##
###_______________________TD définir une fonction et calibrer ses paramètres_________________##
###__________________________________________________________________________________________##
#Raphael PEREZ & Christian FOURNIER , Mai 2018


###libraries
require(ggplot2)
require(cowplot)

###environnement de travail
#setwd(dir='/home/perez/Documents/Enseignement/supagro_training_2018/TD_3/')


####____________Etape 1: Définir une fonction_____________####

###importer les données
donLW=read.csv(file = 'dataLW_Maize.csv')

###identifier chaque plante
donLW$plant_id=paste0('plante_',as.numeric(as.factor((paste(donLW$year,donLW$genotype,donLW$Nmax,sep='_')))))

###identifier chaque feuille
donLW$leaf_id=paste('leaf',donLW$N,donLW$plant_id,sep='_')

###Représenter une feuille à partir des données disponibles ####
leaf=12
ggplot(data=donLW[donLW$leaf_id==unique(donLW$leaf_id)[leaf],],aes(x=l,y=w))+
  geom_point()+
  theme_classic()

###Programmer la fonction polynomiale sur R#####
polym=function(x,w0,lm){
  a0=w0
  c0=(w0-1)/(lm**2)
  b0=-2*c0*lm


  c1=-1/(1-lm)**2
  b1=-2*c1*lm
  a1=-b1-c1
  y=ifelse(x<=lm,a0+b0*x+c0*x**2,a1+b1*x+c1*x**2)
  return(y)
}



#### Représenter une feuille avec les différentes valeurs des paramètres####
w0=0.8
lm=0.01
x_sim=seq(from = 0,to = 1,by = 0.02)
y_sim=polym(x = x_sim,w0 = w0,lm = lm)

ggplot()+
  geom_path(aes(y=y_sim,x=x_sim))+
  geom_path(aes(y=-y_sim,x=x_sim))+
  ggtitle(paste('w0=',w0,' lm=',lm))


####____________Etape 2: Calibrer les paramètres d’une fonction_____________####

####Ajuster les valeurs de lm et w0 pour l’ensemble des feuilles mesurées et représenter ces ajustements####

donFit=NULL
pdf('Ajustement_LW.pdf')
for (i in unique(donLW$leaf_id)){
  sub=donLW[donLW$leaf_id==i,]

  fit=nls(data=sub,sub$w~polym(x = sub$l,w0 = w0,lm = lm),start = list(w0=0.5,lm=0.5))
  w0=as.numeric(coef(fit)['w0'])
  lm=as.numeric(coef(fit)['lm'])

  donFit_sub=data.frame(year=unique(sub$year),genotype=unique(sub$genotype),N=unique(sub$N),Nmax=unique(sub$Nmax),lmax=unique(sub$lmax),wmax=unique(sub$wmax),leaf_id=i,w0=w0,lm=lm)
  donFit=rbind(donFit,donFit_sub)

  graph=ggplot()+
    geom_point(data=sub,aes(x=l,y=w))+
    geom_path(aes(x=x_sim,y=polym(x = x_sim,w0=w0,lm=lm),col='red'))+
    theme_classic()+
    theme(legend.position='none')+
    ggtitle(i,paste0('w0=',round(w0,2),' lm=',round(lm,2)))

  print(graph)
}
dev.off()

####Représenter le domaine de variation des 2 paramètres####
par(mfcol=c(1,2))
boxplot(donFit$lm,ylab='lm',las=1)
boxplot(donFit$w0,ylab='w0',las=1)

####estimer leur coefficient de variation cv=sd/mean####
sd(donFit$lm)/mean(donFit$lm)
sd(donFit$w0)/mean(donFit$w0)

####correlation entre les 2 paramètres?####
model=lm(donFit$lm~donFit$w0)
summary(model)
ggplot(data=donFit,aes(x=lm,y=w0))+
  geom_point()+
    geom_smooth(method='lm',se = F)


####Estimer la surface de chaque feuille avec le modele et comparer avec les donnees  ####

###estimation des surfaces reeles par la méthode des trapèzes
###valeurs observées
leaf_area_estim=function(l,lmax,wmax){
    W=w*wmax
    L=l*lmax
    t=NULL
    for (i in 1:(length(W)-1)){
        t[i]=(L[i]-L[i+1])*(W[i]+W[i+1])/2
    }
    area=sum(t)
    return(area)
}
####valeurs prédites
polym_integral=function(w0,lm,lmax,wmax){
  a0=w0
  c0=(w0-1)/(lm**2)
  b0=-2*c0*lm
  
  c1=-1/(1-lm)**2
  b1=-2*c1*lm
  a1=-b1-c1
  
  int=a0*lm + b0/2*lm**2 + c0/3*lm**3 + a1 + b1/2 + c1/3 - a1*lm - b1/2 * lm**2 - c1/3 * lm**3

  return(int*wmax*lmax)
}

donArea=NULL
for (leaf in 1:length(unique(donLW$leaf_id))){
  # leaf=1
  sub=donLW[donLW$leaf_id==unique(donLW$leaf_id)[leaf],]
  l=sub$l
  w=sub$w
  lmax=unique(sub$lmax)
  wmax=unique(sub$wmax)

  w0=donFit[donFit$leaf_id==unique(donLW$leaf_id)[leaf],]$w0
  lm=donFit[donFit$leaf_id==unique(donLW$leaf_id)[leaf],]$lm



  obs = leaf_area_estim(l=l,lmax=lmax,wmax=wmax)
  sim = polym_integral(lmax=lmax,wmax=wmax,w0=w0,lm=lm)

  donArea_sub = data.frame(genotype=unique(donLW[donLW$leaf_id==unique(donLW$leaf_id)[leaf],]$genotype),
                           plant_id=unique(donLW[donLW$leaf_id==unique(donLW$leaf_id)[leaf],]$plant_id),
                           leaf_id=unique(donLW[donLW$leaf_id==unique(donLW$leaf_id)[leaf],]$leaf_id),
                           N=unique(donLW[donLW$leaf_id==unique(donLW$leaf_id)[leaf],]$N),
                             leaf_area_estim=obs,
                             leaf_area_predict=sim)
  donArea=rbind(donArea,donArea_sub)

}


####Root mean square error
f.rmse=function(Y_obs,Y_estim){
  
  nb=length(Y_obs)
  SCR = sum((Y_estim - Y_obs)**2)
  RMSE=sqrt(1/nb * SCR)
  
  return(round(RMSE,2))
}

###Bias
f.bias=function(Y_obs,Y_estim){
  
  nb=length(Y_obs)
  SCR = sum((Y_estim - Y_obs)**2)
  biais = 1/nb * sum(Y_estim - Y_obs)
  
  return(round(biais,2))
}


ggplot(data=donArea,aes(x=leaf_area_estim,y=leaf_area_predict))+
  geom_point()+
  ylab('Simulated leaf area (cm2)')+
  xlab('Observed leaf area (cm2)')+
  geom_abline(slope = 1,intercept = 0)

model=lm(donArea$leaf_area_estim~donArea$leaf_area_predict)
summary(model)

f.rmse(Y_obs = donArea$leaf_area_estim,Y_estim = donArea$leaf_area_predict)
f.bias(Y_obs = donArea$leaf_area_estim,Y_estim = donArea$leaf_area_predict)


####____________Etape 3: Modélisation des effets ontogeniques et genetiques _____________####


###Représenter les boxplots des valeurs de paramètres pour chaque génotype####
ggplot(data=donFit,aes(x=genotype,y=lm,col=genotype))+
  geom_boxplot()+
  theme_classic()

ggplot(data=donFit,aes(x=genotype,y=w0,col=genotype))+
  geom_boxplot()+
  theme_classic()

###Représenter les des valeurs de paramètres en fonction du rang de la feuille (N), du nombre de feuilles total (Nmax), et du de la longeur des feuilles (lmax)####

###N
g1=ggplot(data=donFit,aes(x=N,y=w0))+
  geom_point()
g2=ggplot(data=donFit,aes(x=N,y=lm))+
  geom_point()

###Nmax
g3=ggplot(data=donFit,aes(x=Nmax,y=w0))+
  geom_point()
g4=ggplot(data=donFit,aes(x=Nmax,y=lm))+
  geom_point()

##lmax
g5=ggplot(data=donFit,aes(x=lmax,y=w0))+
  geom_point()
g6=ggplot(data=donFit,aes(x=lmax,y=lm))+
  geom_point()

plot_grid(g1,g2,g3,g4,g5,g6,ncol=2)


####Réaliser une analyse de variances sur lm et w0####
anova_w0=aov(data=donFit,w0~genotype+N+Nmax+lmax)
summary(anova_w0)

anova_lm=aov(data=donFit,lm~genotype+N+Nmax+lmax)
summary(anova_lm)

###verifier les hypothèses####
par(mfcol=c(2,2))
plot(anova_lm)
plot(anova_w0)


####Représenter les longeurs/largeurs de feuilles en fonction de leur rang pour chaque plante observée####
###subset des données pour avoir 1 lmax et 1 wmax par feuille
donL=donLW[donLW$N>3 & donLW$w==0,]

ggplot(data=donL,aes(y=wmax,x=N,col=plant_id))+
  geom_point()+
  geom_line()+
  theme_classic()

ggplot(data=donL,aes(y=lmax,x=N,col=plant_id))+
  geom_point()+
  geom_line()+
  theme_classic()


### Allometrie w / l####
donL$wl_ratio=donL$wmax/donL$lmax
ggplot(data=donL,aes(y=wl_ratio,x=N,col=genotype,pch=plant_id))+
  geom_point()+
  theme_classic()

ggplot(data=donL,aes(y=wmax/lmax,x=genotype,col=genotype))+
  geom_boxplot()+
  theme_classic()

###anova
anova_lw=aov(data=donL,wl_ratio~genotype*N)
summary(anova_lw)

par(mfcol=c(2,2))
plot(anova_lw)

#### Modelisation des effets

lm_wl=lm(data=donL,wl_ratio~N:genotype)
lm_w0=lm(data=donFit,w0~N:genotype)
lm_lm=lm(data=donFit,lm~N:genotype)

summary(lm_wl)

## TODO : plot des fits

### Estimation des surfaces a partir de lmax, rank et des modele ci dessus

##plot obs / sim

##
## inversion en imposant smax et en deformant lmax
