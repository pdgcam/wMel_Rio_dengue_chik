setwd("C:/Users/gri2s/OneDrive - University of Cambridge/R_M2/1_Data/Code_paper/Wolbachia_Rio")

library(INLA)
library(FNN)
library(ggplot2)

df.st = read.csv("data_wmel_rio.csv")
source("model_function_rev.R")


#Run main spatial models presented in the manuscript
#"virus" argument indicates the dataset to use (dengue or chikungunya)
#"wmel_proportion" is a boolean indicating whether wMel % variable should be included in the model
#"in.project" is a boolean indicating whether "In.project" variable should be included in the model
#"prediction_type is an integer ranging from 1 to 4 defining how to remove randomly data from the analysis to assess the model performance and predictive power
#1 : Full set - All cases are used for data prediciton (Figure S9A)
#2 : Half set - 50% of cases are removed at random (Figure S9B)
#3 : Spatial clusters of cells are removed (20%) (Figure S9C)

fit_dengue_wmelprop = INLA_run(virus = "dengue",
                               wmel_proportion = T,
                               in.project = F,
                               prediction_type = 1)

fit_chik_wmelprop = INLA_run(virus = "chik",
                             wmel_proportion = T,
                             in.project = F,
                             prediction_type = 1)

fit_dengue_inproject = INLA_run(virus = "dengue",
                                wmel_proportion = F,
                                in.project = T,
                                prediction_type = 1)

fit_chik_inproject = INLA_run(virus = "chik",
                              wmel_proportion = F,
                              in.project = T,
                              prediction_type = 1)

#Get inferred parameters from the fits
list_outputs = list(
  Wolbachia = list(dengue = fit_dengue_wmelprop$summary.random$Wolbachia.bin,
                   chik = fit_chik_wmelprop$summary.random$Wolbachia.bin),
  time = list(dengue = fit_dengue_wmelprop$summary.random$t,
              chik = fit_chik_wmelprop$summary.random$t),
  space = list(dengue = fit_dengue_wmelprop$summary.random$spatial.field,
               chik = fit_chik_wmelprop$summary.random$spatial.field),
  In.project = list(dengue = fit_dengue_inproject$summary.fixed,
                    chik = fit_chik_inproject$summary.fixed)
)

#Visualise results - example
ggplot()+
  geom_pointrange(data=list_outputs$Wolbachia$dengue ,aes(x=ID, y=exp(mean), ymin=exp(`0.025quant`),ymax=exp(`0.975quant`)))+
  geom_pointrange(data=list_outputs$Wolbachia$chik ,aes(x=ID, y=exp(mean), ymin=exp(`0.025quant`),ymax=exp(`0.975quant`)), col=2)+
  expand_limits(y=0)+
  scale_y_continuous(name="Risk ratio")+
  scale_x_discrete(limits=levels(as.factor(list_outputs$Wolbachia$dengue$ID))[c(8,7,1:6)])






