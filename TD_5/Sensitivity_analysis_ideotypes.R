#----------------------------------------------------------------------------------#
#-----------------------SENSITIVITY ANALYSIS and IDEOTYPES IDENTIFICATION----------##
#----------------------------------------------------------------------------------#

##Raphael PEREZ, Christain FOURNIER, June 2018


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
Directory='/home/perez/Documents/Enseignement/supagro_training_2018/TD_5/'

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


####______________________________________________________________________________####
####----------------------------SENSITIVITY ANALYSIS--------------------------------####
####______________________________________________________________________________####

####group
group='groupe2'

####____________________________load and merge the data______________________________________####

####list of output files per group in the folder 
outputs=list.files(path = paste0(Directory,'outputs'),pattern =group)

###loop to load and estimate sensitivity indices

###dataframe initialization

###dataframe of outputs
don=NULL

###dataframe of sensitivity indices
donS=NULL


###loops on output files
for (i in outputs){
  ###loop on output variables
  sub=read.csv(paste0(Directory,'outputs/',i))
  
  ###add the simu index
  sub$simu=paste0('p',sprintf("%03d", as.numeric(rownames(sub))))
  
  ###estimate LAI & rie
  sub$LAI=sub$Area*sub$density
  sub$RIE=sub$LAI*sub$Ei
  
  ###estimate light intercepted (relatively to incident radiation)
  sub$I=sub$Ei_leaf*sub$Area_leaf
  
  ###merge the data
  don=rbind(don,sub)
  
  ###loop estimate sensitivity indices for each variable
  for (var in c('Ei_leaf','I','RIE')){
    
    ###merge morris plan with output to estimate indices
    out=tell(etude.morris,y= sub[,var])
    # print(etude.morris)
    
    res=data.frame(t(out$ee))
    
    don_out=data.frame(var=var,parameter=parameters,latitude=unique(sub$latitude),density=unique(sub$density),plant_area=unique(sub$plant_area),mu=apply(X=res,MARGIN = 1,mean),mu_star=apply(X=abs(res),MARGIN = 1,mean),sd=apply(X=res,MARGIN = 1,sd))
    
    ####merge the data
    donS=rbind(donS,don_out)
  }
}



####graphics####

###function for ggplot color gradients
gg_color_hue <- function(x,n=3) {
  # n=length(unique(x))
  hues = seq(15, 375, length = n + 1)
  col=hcl(h = hues, l = 65, c = 100)[x]
  return(col)
}

####create factors for visualisation
donS$f_latitude=paste0('lat = ',sprintf("%02d",donS$latitude),' °')
donS$f_area=paste0('area = ',round(donS$plant_area/10000,2),' m2')
donS$f_density=paste0('density = ',sprintf("%02d",donS$density),' plt.m-2')


####sensitivity indices for a given combination of factors####

###select var, latitude, density, plant_area

var='RIE'
latitude=c(0)
plant_area=c(8000)
density=c(4)

ggplot(data=donS[donS$var==var &  donS$density==density & donS$latitude==latitude & donS$plant_area==plant_area,],aes(x=mu_star,y=sd,label=parameter))+
  geom_point()+
  geom_label_repel()+
  xlab(expression(mu^'*'))+
  ylab(expression(sigma))+
  facet_grid(f_density~f_latitude)+
  ggtitle(paste0('Morris indices for ', var)) 


####effect of factor on sensitivity indices####

###tranform to relative indice per environment 
as_rel=donS
as_rel$mu_star_rel=NA
for (d in unique(donS$density)){
  for (l in unique(donS$latitude)){
    for (a in unique(donS$plant_area)){
      for (v in unique(donS$var)){
        as_rel[as_rel$density==d & as_rel$latitude==l & as_rel$plant_area==a & as_rel$var==v,'mu_star_rel']=(as_rel[as_rel$density==d & as_rel$latitude==l & as_rel$plant_area==a& as_rel$var==v,'mu_star']-min(as_rel[as_rel$density==d & as_rel$latitude==l & as_rel$plant_area==a & as_rel$var==v,'mu_star'],na.rm=T))/(max(as_rel[as_rel$density==d & as_rel$latitude==l & as_rel$plant_area==a & as_rel$var==v,'mu_star'],na.rm=T)-min(as_rel[as_rel$density==d & as_rel$latitude==l & as_rel$plant_area==a & as_rel$var==v,'mu_star'],na.rm=T))
      }
    }
  }
}

var='Ei_leaf'
plant_area=c(8000)
density=c(4,9,13)
latitude=c(0,45)


ggplot(data=as_rel[as_rel$var==var &  as_rel$density %in% density & as_rel$latitude %in% latitude,],aes(x=density,y=mu_star_rel,col=parameter,pch=parameter,fill=parameter))+
  geom_point(size=3)+
  geom_line()+
  scale_shape_manual(name='',labels=levels(as_rel[as_rel$var==var &  as_rel$density==density,]$parameter),values= c(15:25))+
  scale_color_manual(name='',labels=levels(as_rel[as_rel$var==var &  as_rel$density==density,]$parameter),values= gg_color_hue(x =c(1:11),n=11))+
  scale_fill_manual(name='',labels=levels(as_rel[as_rel$var==var &  as_rel$density==density,]$parameter),values= gg_color_hue(x =c(1:11),n=11))+
  ylab(expression(mu[rel]^'*'))+
  facet_grid(~f_latitude)+
  ggtitle(paste0('Morris indices for ', var)) 


####______________________________________________________________________________####
####----------------------------IDEOTYPES--------------------------------####
####______________________________________________________________________________####

####create factors for visualisation
don$f_latitude=paste0('lat = ',sprintf("%02d",don$latitude),' °')
don$f_area=paste0('area = ',round(don$plant_area/10000,2),' m2')
don$f_density=paste0('density = ',sprintf("%02d",don$density),' plt.m-2')

###select the variable
var='Ei_leaf'

###subset the best and worst plants architecture for var####
don$type='intermed'

##"loop to extract best and worst per environment
for (d in unique(don$density)){
  for (l in unique(don$latitude)){
    for (a in unique(don$plant_area)){
      don[don$density==d & don$latitude==l & don$plant_area==a & don[don$density==d & don$latitude==l & don$plant_area==a,var]<quantile(x =don[don$density==d & don$latitude==l & don$plant_area==a,var],probs = 0.05),]$type='worst'
      don[don$density==d & don$latitude==l & don$plant_area==a & don[don$density==d & don$latitude==l & don$plant_area==a,var]>quantile(x =don[don$density==d & don$latitude==l & don$plant_area==a,var],probs = 0.95),]$type='best'
      
    }
  }
}

####graphic####
ggplot(data=don,aes_string(x='f_latitude',y=var,col='type'))+
  geom_boxplot()+
  facet_grid(f_area~f_density)


###difference in parameter values between best and worst architecture####

###tranform to relative value of parameter
don_rel=don
for (p in parameters){
    don_rel[,p]=(don[,p]-min(don[,p],na.rm=T))/(max(don[,p],na.rm=T)-min(don[,p],na.rm=T))
}

don_rel_type=do.call('data.frame',aggregate(x = don_rel[,parameters],by = list(latitude=don_rel$latitude,density=don_rel$density,plant_area=don_rel$plant_area,type=don_rel$type),FUN=function(x){c(mean = mean(x,na.rm=T))}))

###only keep parameter mean values vor best and worst archi
don_rel_type=droplevels(don_rel_type[don_rel_type$type %in% c('best','worst'),])

rownames(don_rel_type)=paste(don_rel_type$latitude,don_rel_type$density,don_rel_type$plant_area,don_rel_type$type,sep='_')


####graphic####

###select environements
latitude=0
density=9
plant_area=8000

radial.plot(don_rel_type[don_rel_type$latitude==latitude & don_rel_type$plant_area==plant_area & don_rel_type$density==density,parameters],rp.type='p',labels=parameters,poly.col=alpha(gg_color_hue(x = c(2,1), n=2),alpha=0.3),line.col=gg_color_hue(x = c(2,1), n=2),radial.labels=c(NA,NA,NA,NA,NA,1),radial.lim=c(0,1),lwd=2,main=paste0('lat=',latitude,' density=',density,' area=',plant_area),label.prop=1.2,mar=c(4,4,8,0))
legend(legend=c(paste('best',var),paste('worst',var)),fill=alpha(gg_color_hue(x = c(2,1), n=2),alpha=0.3),bty='n',horiz=F,x=0.6,y=1.5,cex=0.8)


###difference in parameter values between environments for the best architectures####

####graphic####

###select environements
latitude=0
plant_area=8000

don_env=don_rel_type[don_rel_type$plant_area==plant_area & don_rel_type$latitude==latitude & don_rel_type$type=='best',parameters]

radial.plot(don_env,rp.type='p',labels=parameters,line.col=c('yellow','orange','darkred'),radial.labels=c(NA,NA,NA,NA,NA,1),radial.lim=c(0,1),lwd=3,main=paste0('lat=',latitude,' area=',plant_area),label.prop=1.2,mar=c(4,4,8,0))
legend(title='density (plt.m-2)',legend=c('4','9','13'),col=c('yellow','orange','darkred'),lwd=3,bty='n',horiz=F,x=0.6,y=1.6,cex=0.8)


####______________________________________________________________________________####
####----------------------------IDEOTYPES VISUALISATION--------------------------------####
####______________________________________________________________________________####
###slect var
var='RIE'

###select environement
latitude=0
plant_area=8000
density=13

don[don$latitude==latitude & don$density==density & don$plant_area==plant_area & don[don$latitude==latitude & don$density==density & don$plant_area==plant_area,var]==max(don[don$latitude==latitude & don$density==density & don$plant_area==plant_area,var]),c('simu',parameters,var)]


