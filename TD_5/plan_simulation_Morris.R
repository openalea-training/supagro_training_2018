#----------------------------------------------------------------------------------#
#-----------------------------DESIGN MORRIS PLAN FOR SENSITIVITY ANALYSIS-----------------------------##
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


###visualisation of the distribution of sampling values per parameter
plan=as.data.frame(planMorris)
colnames(plan)

hist(plan$incli_top)

plot(plan[,1:4],pch=19,col=alpha(c='red',alpha=0.1))

#####____________________________Groupe1______________________________________####
group='groupe1'

###fixed parameters

###developemental stage (number of leaves)
stages=c(4,8,12)

####date
daydate='1997-07-09'

###final plant area (cm2) (=16 leaves)
plant_area=6000

###density
density=8

###variable parameters

####latitude

for (latitude in c(0,15,30,45,60)){
  
  ###save the file
  filename=paste0('planMorris_',group,'_area',plant_area,'_d',density,'_lat',latitude,'_stade')
  plan=data.frame(as.data.frame(planMorris),plant_area=rep(plant_area,length(stages)*nrow(planMorris)),stage=rep(stages,each=nrow(planMorris)))
  
  write.csv(x=plan,file =paste0(Directory,filename,'.csv'),row.names = F)
  
}


#####____________________________Groupe2______________________________________####
group='groupe2'

###fixed parameters

###developemental stage (number of leaves)
stages=c(4,8,12)

####date
daydate='2018-06-01'

###final plant area (cm2) (=16 leaves)
plant_area=8000

###density
for (density in c(4,9,13)){
  
  
  ###variable parameters
  
  ####latitude
  
  for (latitude in c(0,45)){
    
    ###save the file
    filename=paste0('planMorris_',group,'_area',plant_area,'_d',density,'_lat',latitude,'_stade')
    plan=data.frame(as.data.frame(planMorris),plant_area=rep(plant_area,length(stages)*nrow(planMorris)),stage=rep(stages,each=nrow(planMorris)))
    
    write.csv(x=plan,file =paste0(Directory,filename,'.csv'),row.names = F)
  }
}

#####____________________________Groupe3______________________________________####
group='groupe3'

###fixed parameters

###developemental stage (number of leaves)
stages=c(4,8,12)

####date
daydate='1975-07-18'

###latitude 
latitude=43

###variable parameters

###final plant area (cm2) (=16 leaves)
for (plant_area in c(4000,10000)){
  ###density
  for (density in c(3,7,10)){
    
    ###save the file
    filename=paste0('planMorris_',group,'_area',plant_area,'_d',density,'_lat',latitude,'_stade')
    plan=data.frame(as.data.frame(planMorris),plant_area=rep(plant_area,length(stages)*nrow(planMorris)),stage=rep(stages,each=nrow(planMorris)))
    
    write.csv(x=plan,file =paste0(Directory,filename,'.csv'),row.names = F)
  }
}

