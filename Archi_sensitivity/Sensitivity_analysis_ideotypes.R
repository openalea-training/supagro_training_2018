#----------------------------------------------------------------------------------#
#-----------------------SENSITIVITY ANALYSIS and IDEOTYPES IDENTIFICATION----------##
#----------------------------------------------------------------------------------#

##Raphael PEREZ, Christian FOURNIER, June 2018


#load package for sensitivity
# install.packages(c('sensitivity'))
library(sensitivity)

#load package for graphics
# install.packages(c('ggplot2','cowplot'))
require(ggplot2)
require(cowplot)

# install.packages(c('ggrepel','plotrix','stringr'))
require(ggrepel)
require(plotrix)
require(stringr)

#Directory
# Directory='/home/perez/Documents/Enseignement/supagro_training_2018/TD_5/'
Directory=paste(dirname(rstudioapi::getActiveDocumentContext()$path),sep='')

#===================#
# Theme ggplot
#===================#

myTheme <- theme(
  panel.background=element_rect(fill="transparent", color=NA),  
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.line=element_line(colour="black"), 
  axis.title=element_text(size=16),
  axis.text.y=element_text(size=14, colour="black"), 
  axis.text.x=element_text(size=14, colour="black", angle=0, hjust=0.5),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_line(colour="grey90", size=0.2), 
  legend.position="right", 
  legend.text=element_text(size=14),
  legend.title=element_text(size=14)
)

myThemeBar <- theme(
  panel.background=element_rect(fill="transparent", color=NA),  
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.line=element_line(colour="black"), 
  axis.title=element_text(size=16),
  axis.text.y=element_text(size=14, colour="black"), 
  axis.text.x=element_text(size=14, colour="black", angle=0, hjust=0.5),
  panel.grid.minor = element_blank(), 
  panel.grid.major =element_blank(),
  legend.position="bottom", 
  legend.text=element_text(size=14),
  legend.title=element_text(size=14)
)

###color palette   http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=7
palette7=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628')
paletteRaph=alpha(colour = palette7,alpha = 0.7)


###______________________________________________________________________________####
####----------------------------FUNCTIONS-----------------------####
####______________________________________________________________________________####

f.k=function(area,k){
  y=1-exp(-k*area)
  return(y)
}

# area_sim=seq(0,1,0.1)
# plot(y=f.k(area=area_sim,k=1),x=area_sim,ylim=c(0,1))
# lines(y=f.k(area=area_sim,k=6),x=area_sim)


####______________________________________________________________________________####
####----------------------------METHODE DE MORRIS--------------------------------####
####______________________________________________________________________________####

####----------------------------Design the plan-------------------------------####

##------fixed parameters per sensitivity analysis

#plant_area: plant area at final stage (cm2)
#plant_orientation: plant orientation in azimtuh (degree) (0°= inter-row)

###-----generate variation range of model parameters----####

#plant_height: plant height (cm)

#rmax: relative position of the leaf with highest area

#skew: skewness of leaf area distribution 

#wl_slp: ratio width/length of leaf with leaf rank

#w0_slp: evolution of leaflet width at insertion with leaf rank

#lm_slp: evolution of the position of leaflet maximum width with leaf rank

#incli_top: insertion inclination of the highest leaf on the stem (degree)

#incli_base: insertion inclination of the lowest leaf on the stem (degree)

#infl: curvature coefficient; inflexion point of the logistic relation between insertion angle and leaf tip angle (degree)

#phyllotactic_angle: difference angle between two vertical rows of leaves on the stem (degree) (0° = superposed, 180° = oposed)

#phyllotactic_deviation: leaf azimuth deviation from phyllotactic_angle (degree)

parameters=c('plant_height','rmax','skew','wl_slp','w0_slp','lm_slp','incli_top','incli_base','infl','phyllotactic_angle','phyllotactic_deviation')

Pvar=data.frame(parameter=parameters,min=NA,max=NA) #generate a data frame
rownames(Pvar)=parameters

average_plant=c('plant_height'=200,'rmax'=0.5,'skew'=0.1,'wl_slp'= 0.007,'w0_slp'=0.01,'lm_slp'=-0.02,'incli_top'= 20,'incli_base'=75,'infl'=40,'phyllotactic_angle'=179.9,'phyllotactic_deviation'=30)

Pvar$min=average_plant*0.7
Pvar$max=average_plant*1.3

###----function arguments----####

#number of factors (parameters)
nFact= parameters
#number of trajectories
r=30
#inf limit
binf=Pvar[,'min']; names(binf)=parameters
#sup limit 
bsup=Pvar[,'max'];names(bsup)=parameters
#discretisation levels
Q=5
#discretisation step
step=2

###----create Morris trajectories----####

#set random seed
RNGkind(kind="L'Ecuyer-CMRG")
set.seed(1) 

#plan
etude.morris=morris(model=NULL,factors=as.character(nFact),r=r,design=list(type='oat',levels=Q, grid.jump=step),scale=T,binf= binf,bsup=bsup)

###----save the design----####
planMorris=etude.morris$X

###visualisation of the distribution of sampling values per parameter
plan=as.data.frame(planMorris)
colnames(plan)


####______________________________________________________________________________####
####----------------------------SENSITIVITY ANALYSIS--------------------------------####
####______________________________________________________________________________####

####____________________________load and merge the data______________________________________####

# ####list of output files per group in the folder 
# outputs=list.files(path = paste0(Directory,'/outputs'),pattern ='planMorris')
# 
# ###loop to load and estimate sensitivity indices
# 
# ###dataframe initialization
# 
# ###dataframe of outputs
# don=NULL
# 
# ###dataframe of sensitivity indices
# donS=NULL
# 
# 
# 
# ###loops on output files
# for (i in outputs){
#   print(i)
#   # i=outputs[1]
#   
#   sub=read.csv(paste0(Directory,'/outputs/',i))
#   
#   ###add the simu index
#   sub$simu=paste0('p',sprintf("%03d", as.numeric(sub$X)+1))
#   
#   
#   pdf(file = paste0(Directory,'/outputs/Adjusted_k_area',unique(sub$plant_area),'_d',unique(sub$density),'.pdf'),onefile = T)
#   
#   ###get plant outputs for SA
#   subS=do.call('rbind',lapply(split(sub, sub$simu, drop=TRUE), function(xS) {
#     # xS=sub[sub$simu=='p053',]
#     xS$leaf_rank=1:nrow(xS)
#     xS$relative_rank=xS$leaf_rank/max(xS$leaf_rank)
#     
#     xS$I=xS$area*xS$Ei
#     xS$cumulArea=cumsum(xS$area)/sum(xS$area)
#     xS$cumulI=cumsum(xS$I)/sum(xS$I)
#     
#     res=unique(xS[,c('simu','density','plant_area',parameters)])
#     
#     res$plant_leaf_area=sum(xS$area)
#     res$plant_I=sum(xS$I)
#     res$plant_Ei=res$plant_I/res$plant_leaf_area
#     res[,c(paste0('Ei_leaf_',1:16))]=xS$Ei
#     
#     ###divide the plant in 3 layers
#     res$area_low=sum(xS[2:6,'area'])
#     res$area_mid=sum(xS[7:11,'area'])
#     res$area_high=sum(xS[12:16,'area'])
#     
#     res$I_low=sum(xS[2:6,'area']*xS[2:6,'Ei'])
#     res$I_mid=sum(xS[7:11,'area']*xS[7:11,'Ei'])
#     res$I_high=sum(xS[12:16,'area']*xS[12:16,'Ei'])
#     
#     res$Ei_low=res$I_low/res$area_low
#     res$Ei_mid=res$I_mid/res$area_mid
#     res$Ei_high=res$I_high/res$area_high
#     
#     
#     ###adjust k
#     fit=nls(data=xS,formula = 1-cumulI~f.k(area = 1-cumulArea,k=k),start = list(k=6))
#     don_fit=data.frame(cumulArea=1-xS$cumulArea,cumulI=predict(fit))
#     k_estim=coef(fit)['k']
#    
#     print(unique(res$simu))
#     
#      # graphK=ggplot()+
#      #  geom_point(data=xS,aes(x=1-cumulArea,y=1-cumulI))+
#      #  geom_path(data=xS,aes(x=1-cumulArea,y=1-cumulI))+
#      #  geom_path(data=don_fit,aes(x=cumulArea,y=cumulI),col='red')+
#      #  myTheme+
#      #  ggtitle(unique(res$simu),paste0('k = ',round(k_estim,2)))+
#      #  xlab('Cumulated leaf area from top')+
#      #  ylab('Cumulated interception')
#      # 
#      # print(graphK)
#      # 
#     
#      res$k=k_estim
#      
#      ##approx rh50
#      res$rh50=approx(y =xS[,'relative_rank'] ,x = xS[,'cumulArea'],xout=0.5)$y
#      
#      # ggplot(data=xS,aes(x=xS[,'cumulArea'],y=xS[,'relative_rank']))+
#        # geom_point()
#      
#      
#     return(res)
#     }))
#     
#   ###estimate LAI & rie
#   subS$LAI=subS$plant_leaf_area*subS$density
#   subS$RIE=subS$LAI*subS$plant_I
#   
#   ###merge the data
#   don=rbind(don,subS)
#   
#   ###loop to estimate sensitivity indices for each variable
#   for (var in c('LAI','RIE','plant_Ei',paste0('Ei_leaf_',1:16),'area_low','area_mid','area_high','I_low','I_mid','I_high','Ei_low','Ei_mid','Ei_high','k','rh50')){
#     
#     ###merge morris plan with output to estimate indices
#     out=tell(etude.morris,y= subS[,var])
#     # print(etude.morris)
#     
#     result=data.frame(t(out$ee))
#     
#     don_out=data.frame(var=var,parameter=parameters,density=unique(subS$density),plant_area=unique(subS$plant_area),mu=apply(X=result,MARGIN = 1,mean),mu_star=apply(X=abs(result),MARGIN = 1,mean),sd=apply(X=result,MARGIN = 1,sd))
#     
#     ####merge the data
#     donS=rbind(donS,don_out)
#   }
#   dev.off()
# }
# 
# save(x=donS,file =paste0(Directory,'/Archi_sensitivity.RData'))

load(file =paste0(Directory,'/Archi_sensitivity.RData'))


####graphics####

###function for ggplot color gradients
gg_color_hue <- function(x,n=3) {
  # n=length(unique(x))
  hues = seq(15, 375, length = n + 1)
  col=hcl(h = hues, l = 65, c = 100)[x]
  return(col)
}

####create factors for visualisation
donS$f_area=paste0('area = ',round(donS$plant_area/10000,2),' m2')
donS$f_density=paste0('density = ',sprintf("%g",donS$density),' plt.m-2')


###relatin rh50 rmax
don$f_area=paste0('area = ',round(don$plant_area/10000,2),' m2')
don$f_density=paste0('density = ',sprintf("%g",don$density),' plt.m-2')

# ggplot(data=don, aes(x=rmax,y=rh50))+
#   geom_point()+
#   geom_smooth()+
#   facet_grid(f_density~f_area)

####sensitivity indices for a given combination of factors####

###select var, latitude, density, plant_area

# var='RIE'
# plant_area=c(8000)
# density=c(5.5)
# 
# ggplot(data=donS[donS$var==var &  donS$density==density  & donS$plant_area==plant_area,],aes(x=mu_star,y=sd,label=parameter))+
#   geom_point()+
#   geom_label_repel()+
#   xlab(expression(mu^'*'))+
#   ylab(expression(sigma))+
#   facet_grid(f_density~f_area)+
#   ggtitle(paste0('Morris indices for ', var)) 
# 
# 
# var='Ei_leaf_9'
# plant_area=c(6000,8000,10000)
# density=c(5.5,7,9,11)
# 
# ggplot(data=donS[donS$var==var &  donS$density %in% density  & donS$plant_area %in% plant_area,],aes(x=as.factor(density),y=parameter,fill=mu_star))+
#   geom_tile()+
#   ylab('')+
#   xlab(expression(density(plants.m^-2)))+
#   facet_grid(~f_area)+
#   scale_fill_gradient2(name = expression(mu^'*'),low='darkblue',mid = 'white',high="darkred")+
#   ggtitle(paste0('Morris indices for ', var)) +
#   myTheme+
#   theme(legend.position='right')




####effect of factor on sensitivity indices####

###tranform to relative indice per environment 
as_rel=donS
as_rel$mu_star_rel=NA
for (d in unique(donS$density)){
  for (a in unique(donS$plant_area)){
    for (v in unique(donS$var)){
      as_rel[as_rel$density==d  & as_rel$plant_area==a & as_rel$var==v,'mu_star_rel']=(as_rel[as_rel$density==d  & as_rel$plant_area==a & as_rel$var==v,'mu_star']-min(as_rel[as_rel$density==d  & as_rel$plant_area==a & as_rel$var==v,'mu_star'],na.rm=T))/(max(as_rel[as_rel$density==d  & as_rel$plant_area==a & as_rel$var==v,'mu_star'],na.rm=T)-min(as_rel[as_rel$density==d  & as_rel$plant_area==a & as_rel$var==v,'mu_star'],na.rm=T))
    }
  }
}

# var='LAI'
# plant_area=c(6000,8000,10000)
# density=c(5.5,7,9,11)

# ggplot(data=as_rel[as_rel$var==var &  as_rel$density %in% density ,],aes(x=density,y=mu_star_rel,col=parameter,pch=parameter,fill=parameter))+
#   geom_point(size=3)+
#   geom_line()+
#   scale_shape_manual(name='',labels=levels(as_rel[as_rel$var==var &  as_rel$density==density,]$parameter),values= c(15:25))+
#   scale_color_manual(name='',labels=levels(as_rel[as_rel$var==var &  as_rel$density==density,]$parameter),values= gg_color_hue(x =c(1:11),n=11))+
#   scale_fill_manual(name='',labels=levels(as_rel[as_rel$var==var &  as_rel$density==density,]$parameter),values= gg_color_hue(x =c(1:11),n=11))+
#   ylab(expression(mu[rel]^'*'))+
#   facet_grid(~f_area)+
#   ggtitle(paste0('Morris indices for ', var)) 


#####color map of sensitivity indices####
plant_area=sort(unique(donS$plant_area))
density=sort(unique(donS$density))

pdf(file = paste0(Directory,'/outputs/Map_AS_mu.pdf'),onefile = T)
for (v in unique(as_rel$var)){
  # var='Ei_plant'

  gmap=ggplot(data=as_rel[as_rel$var==v &  as_rel$density %in% density  & as_rel$plant_area %in% plant_area,],aes(x=as.factor(density),y=parameter,fill=mu_star_rel))+
    geom_tile()+
    ylab('')+
    xlab(expression(density(plants.m^-2)))+
    facet_grid(~f_area)+
    scale_fill_gradient2(name = expression(mu[rel]^'*'),low='darkblue',mid = 'white',high="darkred")+
    ggtitle(paste0('Morris indices for ', v)) +
    myTheme+
    theme(legend.position='right')
  
  print(gmap)
}
dev.off()

pdf(file = paste0(Directory,'/outputs/Map_AS_sigma.pdf'),onefile = T)
for (v in unique(as_rel$var)){
  # var='Ei_plant'
  
  gmap=ggplot(data=as_rel[as_rel$var==v &  as_rel$density %in% density  & as_rel$plant_area %in% plant_area,],aes(x=as.factor(density),y=parameter,fill=sd))+
    geom_tile()+
    ylab('')+
    xlab(expression(density(plants.m^-2)))+
    facet_grid(~f_area)+
    scale_fill_gradient2(name = expression(mu[rel]^'*'),low='darkblue',mid = 'white',high="darkred")+
    ggtitle(paste0('Morris indices for ', v)) +
    myTheme+
    theme(legend.position='right')
  
  print(gmap)
}
dev.off()











####______________________________________________________________________________####
####----------------------------IDEOTYPES--------------------------------####
####______________________________________________________________________________####
#load package for experimanental design
library(lhs)
#####______________________LHS design____________________####

####function to generate the LHS plan
RandomLHS=function(factors,distribParameters,size,preserveDraw=FALSE){
  set.seed(1)
  nbf=length(factors)
  design=randomLHS(n=size,k=nbf,preserveDraw = preserveDraw)
  for (i in 1:nbf){
    design[,i]=distribParameters[[1]][[i]]+design[,i]*(distribParameters[[2]][[i]]-distribParameters[[1]][[i]])
  }
  colnames(design)=factors
  resultats=as.data.frame(design)
  return(resultats)
}

####----------------------------Design the plan-------------------------------####

parameters=c('plant_height','rmax','skew','wl_slp','w0_slp','lm_slp','incli_top','incli_base','infl','phyllotactic_angle','phyllotactic_deviation')

parametersLHS=c('plant_height','rmax','incli_base','infl','phyllotactic_angle')

fix_parameters=parameters[!(parameters %in% parametersLHS)]

PvarLHS=data.frame(parameter=parametersLHS,min=NA,max=NA) #generate a data frame
rownames(PvarLHS)=parametersLHS

average_plant=c('plant_height'=200,'rmax'=0.5,'skew'=0.1,'wl_slp'= 0.007,'w0_slp'=0.01,'lm_slp'=-0.02,'incli_top'= 20,'incli_base'=75,'infl'=40,'phyllotactic_angle'=179.9,'phyllotactic_deviation'=30)

PvarLHS$min=average_plant[parametersLHS]*0.5
PvarLHS$max=average_plant[parametersLHS]*1.5


size=5**5
distribParameters=PvarLHS[,c('min','max')]

planLHS=RandomLHS(factors = parametersLHS,size=size,distribParameters = distribParameters,preserveDraw=FALSE)

for (f in fix_parameters){
  planLHS[f]=average_plant[f]
}

planLHS_10000=planLHS
planLHS_10000$plant_area=10000

planLHS_6000=planLHS
planLHS_6000$plant_area=6000
###----save the design----####
filename='planLHS_6000'
write.csv(x=planLHS_6000,file =paste(Directory,'/',filename,'.csv',sep=''),row.names = F)

filename='planLHS_10000'
write.csv(x=planLHS_10000,file =paste(Directory,'/',filename,'.csv',sep=''),row.names = F)

####list of output files per group in the folder 
outputsLHS=list.files(path = paste0(Directory,'/outputs'),pattern ='planLHS')

###loop to load and estimate sensitivity indices

###dataframe initialization

###dataframe of outputs
donLHS=NULL

###loops on output files
for (i in outputsLHS){
  print(i)
  # i=outputsLHS[1]
  
  sub=read.csv(paste0(Directory,'/outputs/',i))
  
  ###add the simu index
  sub$simu=paste0('p',sprintf("%03d", as.numeric(sub$X)+1))
  
  
  # pdf(file = paste0(Directory,'/outputs/Adjusted_k_area',unique(sub$plant_area),'_d',unique(sub$density),'.pdf'),onefile = T)
  
  ###get plant outputs for SA
  subS=do.call('rbind',lapply(split(sub, sub$simu, drop=TRUE), function(xS) {
    # xS=sub[sub$simu=='p053',]
    xS$leaf_rank=1:nrow(xS)
    xS$relative_rank=xS$leaf_rank/max(xS$leaf_rank)
    
    xS$I=xS$area*xS$Ei
    xS$cumulArea=cumsum(xS$area)/sum(xS$area)
    xS$cumulI=cumsum(xS$I)/sum(xS$I)
    
    res=unique(xS[,c('simu','density','plant_area',parameters)])
    
    res$plant_leaf_area=sum(xS$area)
    res$plant_I=sum(xS$I)
    res$plant_Ei=res$plant_I/res$plant_leaf_area
    res[,c(paste0('Ei_leaf_',1:nrow(xS)))]=xS$Ei
    
    lim=ceiling(nrow(xS)/3)
    ###divide the plant in 3 layers
    res$area_low=sum(xS[2:lim,'area'])
    res$area_mid=sum(xS[(lim+1):(2*lim-1),'area'])
    res$area_high=sum(xS[(2*lim):nrow(xS),'area'])
    
    res$I_low=sum(xS[2:lim,'area']*xS[2:lim,'Ei'])
    res$I_mid=sum(xS[(lim+1):(2*lim-1),'area']*xS[(lim+1):(2*lim-1),'Ei'])
    res$I_high=sum(xS[(2*lim):nrow(xS),'area']*xS[(2*lim):nrow(xS),'Ei'])
    
    res$Ei_low=res$I_low/res$area_low
    res$Ei_mid=res$I_mid/res$area_mid
    res$Ei_high=res$I_high/res$area_high
    
    
    ###adjust k
    fit=nls(data=xS,formula = 1-cumulI~f.k(area = 1-cumulArea,k=k),start = list(k=6))
    don_fit=data.frame(cumulArea=1-xS$cumulArea,cumulI=predict(fit))
    k_estim=coef(fit)['k']
    
    print(unique(res$simu))
    
    # graphK=ggplot()+
    #  geom_point(data=xS,aes(x=1-cumulArea,y=1-cumulI))+
    #  geom_path(data=xS,aes(x=1-cumulArea,y=1-cumulI))+
    #  geom_path(data=don_fit,aes(x=cumulArea,y=cumulI),col='red')+
    #  myTheme+
    #  ggtitle(unique(res$simu),paste0('k = ',round(k_estim,2)))+
    #  xlab('Cumulated leaf area from top')+
    #  ylab('Cumulated interception')
    # 
    # print(graphK)

    
    res$k=k_estim
    
    ##approx rh50
    res$rh50=approx(y =xS[,'relative_rank'] ,x = xS[,'cumulArea'],xout=0.5)$y
    
    # ggplot(data=xS,aes(x=xS[,'cumulArea'],y=xS[,'relative_rank']))+
    # geom_point()
    
    
    return(res)
  }))
  
  ###estimate LAI & rie
  subS$LAI=subS$plant_leaf_area*subS$density
  subS$RIE=subS$LAI*subS$plant_I
  
  ###merge the data
  donLHS=rbind(donLHS,subS)
  # dev.off()
}


####create factors for visualisation

donLHS$f_area=paste0('area = ',round(donLHS$plant_area/10000,2),' m2')
donLHS$f_density=as.factor(donLHS$density)
levels(donLHS$f_density)=list('d = 5.5 plt.m-2'='5.5','d = 7 plt.m-2'='7','d = 9 plt.m-2'='9','d = 11 plt.m-2'='11','d = 12.5 plt.m-2'='12.5')

###----fit the polynomial metamodel----####
###select the variable
var='Ei_leaf_9'

#complete model
for (d in unique(donLHS$density)){
  for (a in unique(donLHS$plant_area)){
    MM_poly_total=lm(formula=donLHS[donLHS$density==d & donLHS$plant_area==a,var]~polym(plant_height,rmax,incli_base,infl,phyllotactic_angle,degree=3),data=donLHS[donLHS$density==d & donLHS$plant_area==a,])
    
    summary(lm(formula=donLHS[donLHS$density==d & donLHS$plant_area==a,var]~predict(MM_poly_total)))
    
    r2_total=summary(MM_poly_total)$adj.r.squared
    print(paste(var,'  r2=',r2_total))
    
    #model adjustment
    
    plot(predict(MM_poly_total)~donLHS[donLHS$density==d & donLHS$plant_area==a,var],xlab='Meta model polynome',ylab='model')
    abline(a=0,b=1,col=2)
  }
}


#estimate principal effect and interaction
tableMM=NULL

for (d in unique(donLHS$density)){
  for (a in unique(donLHS$plant_area)){
    donLHS_sub=donLHS[donLHS$density==d & donLHS$plant_area==a,]
    V=donLHS_sub[,parametersLHS]
    for (i in 1:ncol(V)){
      v=paste(colnames(V)[i])
      
      model_seul=lm(formula=donLHS_sub[,var]~polym(V[,i],degree=3),data=donLHS_sub)
      r2_seul=summary(model_seul)$adj.r.squared
      
      vecteur=data.frame(r=c(1:ncol(V)),c=c(1:ncol(V)))
      w=which(vecteur$c!=i)  
      
      model_sauf=lm(formula=donLHS_sub[,var]~polym(V[,w[1]],V[,w[2]],V[,w[3]],V[,w[4]],degree=3),data=donLHS_sub)
      r2_sauf=summary(model_sauf)$adj.r.squared
      
      r2_spe=r2_total-r2_sauf
      r2_int=r2_spe-r2_seul
      
      tableMM_sub=data.frame(Alone=r2_seul,Specific=r2_spe,Total=r2_spe+abs(r2_int), Interaction=r2_int,row.names=v)
      tableMM_sub$density=d
      tableMM_sub$plant_area=a
      tableMM=rbind(tableMM,tableMM_sub)
    }
    
    assign(paste('table',var),tableMM)
    
    print(tableMM)
  }
}



###----graphs----####

# ggplot(data=donLHS,aes(x=rh50,y=rmax))+
#   geom_point()+
#   geom_smooth(method='lm')

#variance decomposition and parameter effects
# barplot(t(tableMM[,c('Specific','Interaction')]),las=2,main=var)

data_bar=data.frame(density=rep(tableMM[,'density'],2),plant_area=rep(tableMM[,'plant_area'],2),r2=c(tableMM[,'Specific'],tableMM[,'Interaction']),Effect=c(rep('Specific',nrow(tableMM)),rep('Interaction',nrow(tableMM))),Parameter=rep(parametersLHS,2))

data_bar$f_area=paste0('area = ',round(data_bar$plant_area/10000,2),' m2')
data_bar$f_density=as.factor(data_bar$density)
levels(data_bar$f_density)=list('d = 5.5 plt.m-2'='5.5','d = 7 plt.m-2'='7','d = 9 plt.m-2'='9','d = 11 plt.m-2'='11','d = 12.5 plt.m-2'='12.5')

graphDecomp=ggplot() + 
  geom_hline(yintercept=r2_total,lty=2)+
  geom_bar(data=data_bar, aes(fill=Effect,y=r2, x=Parameter), stat="identity",col=1,lwd=0.1)+
  ylab('')+
  ylim(0,1)+
  coord_flip()+
  ggtitle(var)+
  facet_grid(f_area~f_density)+
  myThemeBar+
  scale_fill_manual(name='Effect',labels=c(expression(Interaction),expression(Specific)),values =paletteRaph[c(1,3)])+
  theme(legend.position='right')

print(graphDecomp)




###subset the best and worst plants architecture for var####
donLHS$type='intermed'

##"loop to extract best and worst per environment
don_ideo=NULL
for (d in unique(donLHS$density)){
  for (a in unique(donLHS$plant_area)){
    donLHS[donLHS$density==d  & donLHS$plant_area==a & donLHS[donLHS$density==d & donLHS$plant_area==a,var]<quantile(x =donLHS[donLHS$density==d &  donLHS$plant_area==a,var],probs = 0.05),]$type='worst'
    donLHS[donLHS$density==d & donLHS$plant_area==a & donLHS[donLHS$density==d & donLHS$plant_area==a,var]>quantile(x =donLHS[donLHS$density==d & donLHS$plant_area==a,var],probs = 0.95),]$type='best'
    
    sub_max=donLHS[donLHS$density==d  & donLHS$plant_area==a & donLHS[donLHS$density==d & donLHS$plant_area==a,var]==max(donLHS[donLHS$density==d &  donLHS$plant_area==a,var]),]
    don_ideo=rbind(don_ideo,sub_max)
  }
}

####graphic####
ggplot(data=donLHS,aes_string(x='rmax',y='Ei_mid',col='area_mid'))+
  geom_point()+
  scale_color_gradient(name = expression(area(m^2)),low='yellow',high="darkred")+
  myTheme+
  facet_grid(f_area~f_density)

ggplot(data=donLHS,aes_string(x='rmax',y='Ei_mid',col='type'))+
  geom_point()+
  myTheme+
  facet_grid(f_area~f_density)

ggplot(data=donLHS,aes_string(x='area_mid',y='I_mid',col='rmax'))+
  geom_point()+
  scale_color_gradient(name = expression(rmax),low='yellow',high="darkred")+
  myTheme+
  facet_grid(f_area~f_density)


ggplot(data=donLHS[donLHS$type=='best',],aes_string(x='density',y='I_mid',col='simu',group='simu'))+
  geom_point()+
  geom_line()+
  myTheme+
  theme(legend.position='none')

ggplot(data=donLHS[donLHS$simu %in% unique(don_ideo$simu),],aes_string(x='density',y='I_mid',col='simu',group='simu'))+
  geom_point()+
  geom_line()+
  myTheme


###difference in parameter values between best and worst architecture####
donLHS_type=do.call('data.frame',aggregate(x = donLHS[,parametersLHS],by = list(density=donLHS$density,plant_area=donLHS$plant_area,type=donLHS$type),FUN=function(x){c(mean = mean(x,na.rm=T),sd=sd(x,na.rm=T))}))


donLHS_type[donLHS_type$type=='best' & donLHS_type$density==12.5,]


ggplot(data=donLHS_type[donLHS_type$type %in% c('best','worst'),],aes(x=density,y=rmax.mean,col=type))+
  geom_pointrange(aes(ymin=rmax.mean-rmax.sd,ymax=rmax.mean+rmax.sd))+
  geom_line()+
  facet_wrap(~plant_area)


###tranform to relative value of parameter
donLHS_rel=donLHS
for (p in parametersLHS){
    donLHS_rel[,p]=(donLHS[,p]-min(donLHS[,p],na.rm=T))/(max(donLHS[,p],na.rm=T)-min(donLHS[,p],na.rm=T))
}

donLHS_rel_type=do.call('data.frame',aggregate(x = donLHS_rel[,parametersLHS],by = list(density=donLHS_rel$density,plant_area=donLHS_rel$plant_area,type=donLHS_rel$type),FUN=function(x){c(mean = mean(x,na.rm=T),sd=sd(x,na.rm=T))}))

###only keep parameter mean values vor best and worst archi
donLHS_rel_type=droplevels(donLHS_rel_type[donLHS_rel_type$type %in% c('best','worst'),])

rownames(donLHS_rel_type)=paste(donLHS_rel_type$density,donLHS_rel_type$plant_area,donLHS_rel_type$type,sep='_')


####graphic####

###select environements

don_rad=donLHS_rel_type[donLHS_rel_type$density==5.5,]

radial.plot(don_rad[,paste0(parametersLHS,'.mean')],rp.type='p',labels=parametersLHS,poly.col=alpha(palette7[1:2],alpha=0.3),line.col=palette7[1:2],lwd=2,main=paste0('area=',unique(don_rad$plant_area)),label.prop=1.2,mar=c(4,4,8,0))
legend(legend=c(paste('best',var),paste('worst',var)),fill=alpha(gg_color_hue(x = c(2,1), n=2),alpha=0.3),bty='n',horiz=F,x=0.6,y=1.5,cex=0.8)




###



