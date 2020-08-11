############################################################################
#  SIMULATION STUDY - QUADRATIC LOGISTIC REGRESSION
############################################################################
set.seed(9517)
###required packages
library(lme4)
library (survival)
library(matrixStats)
library(ggplot2)
memory.limit(size=4000000000)
# change working directory to current directory (if working with Rstudio)
if(Sys.getenv('RSTUDIO') == '1')
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions4simulation sub-sampling techniques") # functions for simulations
#####################

#sim.dat <- create_dat(num_cluster = 500, cluster_size = 500)
#data_use <- sampling(samp_method=1,dat=sim.dat,ratio=2)


#############################################################################
# Starting main simulation to compare the random effects and CML estimators
# using different types of within-cluster case-control sub-sampling and sub-sampling across the whole data frame
#
#
# input:
###########
#num_sim is the number of simulations
#samp_mode - see sampling function
# samp_meth: the sampling mechanisim generated the data.
#            one of the following options should be used
#                 0 no sampling - original data
#                 1 is samp_meth1 - random sampling
#                 2 is samp_meth2 - case control ignoring clusters
#                 3 is samp_meth3 - case control prop to cluster size
#                 4 is samp_meth4 - case control prop to no. cluster cases
# es_method: 1 - random effects using glmer or 2 - cml using clogit 
# ratio is the ratio between the number of cases to the number of controls used by the sampling method
# Output:
###########
#Estimators either by random effects or conditional logistic regression

################################################################

mainfun<-function(num_sim,samp_mode,es_mode,ratio=4)
{ ptm <- proc.time()

res <-matrix(NA, num_sim,4) 
for (u in (1:num_sim))
{print(u)
  
  data_org<-create_dat()

    if(samp_mode==0)
    {
      data_use<-data_org 
    }
    
    if(samp_mode==1)
    {
      data_use<-sampling(samp_method=1,dat=data_org,ratio)
    }
    if(samp_mode==2)
    {
      data_use<-sampling (samp_method=2,dat=data_org,ratio) 
    }
    if(samp_mode==3)
    {
      data_use<-sampling (samp_method=3,dat=data_org,ratio) 
    }
    if(samp_mode==4)
    {
      data_use<-sampling (samp_method=4,dat=data_org,ratio)
    }
    
    if(es_mode==1)
    {
      estimator<-est(da=data_use,es_method=1,samp_meth=samp_mode) 
    }
    if(es_mode==2)
    {
      estimator<-est(da=data_use,es_method=2,samp_meth=samp_mode) 
    }
    res[u,] <- estimator

  
}

time<-as.vector(proc.time() - ptm)[3]
return(list(res,time))
}

#crearing figure 6 in the article
out_1_0_1<-mainfun(300,samp_mode=0,es_mode=1,ratio=6) 
out_1_1_1<-mainfun(300,samp_mode=1,es_mode=1,ratio=6) 
out_1_2_1<-mainfun(300,samp_mode=2,es_mode=1,ratio=6) 
out_1_3_1<-mainfun(300,samp_mode=3,es_mode=1,ratio=6) 
out_1_4_1<-mainfun(300,samp_mode=4,es_mode=1,ratio=6) 

out_1_0_2<-mainfun(300,samp_mode=0,es_mode=2,ratio=6) 
out_1_1_2<-mainfun(300,samp_mode=1,es_mode=2,ratio=6) 
out_1_2_2<-mainfun(300,samp_mode=2,es_mode=2,ratio=6) 
out_1_3_2<-mainfun(300,samp_mode=3,es_mode=2,ratio=6) 
out_1_4_2<-mainfun(300,samp_mode=4,es_mode=2,ratio=6) 


out <- list()
out[[1]] <- out_1_0_1
out[[2]] <- out_1_1_1
out[[3]] <- out_1_2_1
out[[4]] <- out_1_3_1
out[[5]] <- out_1_4_1

out[[6]] <-out_1_0_2
out[[7]] <-out_1_1_2
out[[8]] <-out_1_2_2
out[[9]] <-out_1_3_2
out[[10]] <-out_1_4_2


###########presenting results###########
estimators <- list()
V_estimators <- list()
bias <- list()
MSE <- list()
time<- list()
beta_1_r<-2
beta_2_r<- -2
real<-c(beta_1_r,beta_2_r)
for (i in 1:10)
{ if(i<6)
  {
    all <- out[[i]][[1]][,2:3]
}
  if(i>5)
  {
    all <- out[[i]][[1]][,1:2]
  }
  
estimators [[i]]<- colMeans(all)
V_estimators[[i]]<-colVars(all)
bias[[i]]<-estimators[[i]]-real
MSE[[i]]<-(bias[[i]])^2+V_estimators[[i]]
time[[i]]<-out[[i]][[2]]

}

#MSE plot x
MSE_plot_x<-numeric()
for(i in 1:10)
{ 
  MSE_plot_x<-c(MSE_plot_x,abs(MSE[[i]][1]))
  
}

names<-c("(0)r","(i)r","(ii)r","(iii)r","(iv)r","(0)c","(i)c","(ii)c","(iii)c","(iv)c")

MSE_plot_x<-round(MSE_plot_x,3)
M_plotdat_x<-cbind(names,MSE_plot_x)
M_plotdat_x<-data.frame(M_plotdat_x)
M_plotdat_x$names <- factor(M_plotdat_x$names, levels = M_plotdat_x$names)
ggplot(M_plotdat_x) + geom_bar(aes(x = M_plotdat_x$names, y = M_plotdat_x$MSE_plot_x), stat = "identity")+xlab("Sampling and estimation method") +ylab("MSE of the X coefficient")

MSE_plot_x2<-numeric()
for(i in 1:10)
{ 
  MSE_plot_x2<-c(MSE_plot_x2,abs(MSE[[i]][2]))
  
}


MSE_plot_x2<-round(MSE_plot_x2,3)
M_plotdat_x2<-cbind(names,MSE_plot_x2)
M_plotdat_x2<-data.frame(M_plotdat_x2)
M_plotdat_x2$names <- factor(M_plotdat_x2$names, levels = M_plotdat_x2$names)
ggplot(M_plotdat_x2) + geom_bar(aes(x = M_plotdat_x2$names, y = M_plotdat_x2$MSE_plot_x2), stat = "identity")+xlab("Sampling and estimation method") +ylab("MSE of the X^2 coefficient")


