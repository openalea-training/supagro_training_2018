###__________________________________________________________________________________________##
###_______________________ECUE2:DEMARCHE DE MODELISATION_________________________##
###_______________________TD définir une fonction et calibrer ses paramètres_________________##
###__________________________________________________________________________________________##
#Raphael PEREZ & Christian FOURNIER , Mai 2018


###libraries
require(ggplot2)
require(cowplot)

###environnement de travail
setwd(dir='/home/perez/Documents/Enseignement/supagro_training_2018/TD_3/')


####____________Etape 1: Définir une fonction_____________####

###importer les données
donLW=read.csv(file = 'dataLW_Maize.csv')

###identifier chaque plante
donLW$plant_id=paste(donLW$year,donLW$genotype,donLW$Nmax,sep='_')

###identifier chaque feuille
donLW$leaf_id=paste(donLW$plant,donLW$genotype,donLW$lmax,donLW$N,sep='_')

###Représenter une feuille à partir des données disponibles ####

###Programmer la fonction polynomiale sur R#####
polym=function(x,w0,lm){
  a0=
  c0=
  b0=
    
    
  c1=
  b1=
  a1=
  y=ifelse(x<=lm,a0+b0*x+c0*x**2,a1+b1*x+c1*x**2)
  return(y)
}

#### Représenter une feuille avec les valeurs ####
w0=0.8
lm=0.01
x_sim=seq(from = 0,to = 1,by = 0.02)
y_sim=


####____________Etape 2: Calibrer les paramètres d’une fonction_____________####

####Ajuster les valeurs de lm et w0 pour l’ensemble des feuilles mesurées et représenter ces ajustements#### 



####Représenter le domaine de variation des 2 paramètres####

####estimer leur coefficient de variation cv=sd/mean####

####correlation entre les 2 paramètres?####

####____________Etape 3: Les paramètres sont-ils génétiques? _____________####


###Représenter les boxplots des valeurs de paramètres pour chaque génotype####


###Représenter les des valeurs de paramètres en fonction du rang de la feuille (N), du nombre de feuilles total (Nmax), et du de la longeur des feuilles (lmax)####

###N


###Nmax


##lmax


####Réaliser une analyse de variances sur lm et w0####

###verifier les hypothèses####


####____________Changer d’échelle et intégrer plusieurs relations allométriques_____________####

####Représenter les longeurs/largeurs de feuilles en fonction de leur rang pour chaque plante observée####
###subset des données pour avoir 1 lmax et 1 wmax par feuille



####Le rapport wmax/lmax est-il constant avec le rang et génétiquement dépendent?####


###anova


####Estimer la surface de chaque feuille et représenter les surfaces de feuilles d’un plante en fonction du rang ####
