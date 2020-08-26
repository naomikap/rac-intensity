############################################################################
#  SIMULATION STUDY - QUADRATIC LOGISTIC REGRESSION
############################################################################
#set.seed(9517)
set.seed(7314)
###required packages
library(lme4)
library (survival)
library(matrixStats)
library(ggplot2)
library(scales)
memory.limit(size=400000000)
# change working directory to current directory (if working with Rstudio)
if(Sys.getenv('RSTUDIO') == '1')
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions4simulation sub-sampling techniques") # functions for simulations
#####################
#############################################################################
# Starting main simulation to compare the random effects and CML estimators
# using different types of within-cluster case-control sub-sampling and sub-sampling across the whole data frame
#
#
# input:
###########
#num_sim is the number of simulations
#samp - see sampling function
# samp: vector of the sampling mechanisims that generated the data.
#            a vector containing some of the following options 
#                 0 no sampling - original data
#                 1 is samp_meth1 - random sampling
#                 2 is samp_meth2 - case control ignoring clusters
#                 3 is samp_meth3 - case control prop to cluster size
#                 4 is samp_meth4 - case control prop to no. cluster cases
# es: a vector of the following options:
#                 1 - random effects using glmer 
#                 2 - cml using clogit 
# ratio is the ratio between the number of cases to the number of controls used by the sampling method
# Output:
###########
#res: A list containing an element for each combinations of samp and es 
#      Each element is a matrix where each row contains the estimators for the corresponding iteration
#mode: Each column of this matrix is the combination of samp and es used in the corresponding element of res
################################################################

mainfun<-function(num_sim,samp,es,ratio=4) {
ptm <- proc.time()

mode <- rbind(rep(samp,length(es)),
              as.vector(t(matrix(rep(es,length(samp)),ncol=length(samp),byrow=FALSE))))
n_mode <- ncol(mode)
res <- rep(list(numeric()),n_mode)

for (u in (1:num_sim)) {
  data_org<-create_dat()
  for(i in 1:n_mode) {
    print(u)
    data_use<-sampling(samp_method=mode[1,i],dat=data_org,ratio)
    res[[i]] <-rbind(res[[i]],est(da=data_use,es_method=mode[2,i],samp_meth=mode[1,i]))
  }
}


time<-as.vector(proc.time() - ptm)[3]
return(list(res=res,mode=mode,time=time))
}

#crearing figure 6 in the article
results <- mainfun(300,samp=0:4,es=1:2,ratio=6) 
                                                                                                                                                    

out <- results$res
mode <- results$mode

###########presenting results###########
estimators <- list()
V_estimators <- list()
bias <- list()
MSE <- list()
time <- list()
beta_1_r <-  2
beta_2_r <- -2
real<-c(beta_1_r,beta_2_r)
for (i in 1:ncol(mode)) {
  if(mode[2,i]==1) {
    all <- out[[i]][,2:3]
  }
  if(mode[2,i]==2)
  {
    all <- out[[i]][,1:2]
  }

estimators[[i]] <- colMeans(all)
V_estimators[[i]] <- colVars(all)
bias[[i]] <- estimators[[i]]-real
MSE[[i]] <- (bias[[i]])^2+V_estimators[[i]]
time[[i]] <- out[[i]][[2]]

}

#MSE plot x
MSE_plot_x<-numeric()
for(i in 1:ncol(mode))
{ 
  MSE_plot_x<-c(MSE_plot_x,abs(MSE[[i]][1]))

}

names<-c("(0)r","(i)r","(ii)r","(iii)r","(iv)r","(0)c","(i)c","(ii)c","(iii)c","(iv)c")

M_plotdat_x<-cbind(names,MSE_plot_x)
M_plotdat_x<-data.frame(M_plotdat_x)
M_plotdat_x$names <- factor(M_plotdat_x$names, levels = M_plotdat_x$names)
ggplot(M_plotdat_x) + geom_bar(aes(x = names, y = as.numeric(paste(MSE_plot_x))), stat = "identity")+xlab("Sampling and estimation method") +ylab("MSE of the X coefficient") + coord_cartesian(ylim = c(0.015, 0.025))

MSE_plot_x2<-numeric()
for(i in 1:ncol(mode))
{ 
  MSE_plot_x2<-c(MSE_plot_x2,abs(MSE[[i]][2]))
  
}


M_plotdat_x2<-cbind(names,MSE_plot_x2)
M_plotdat_x2<-data.frame(M_plotdat_x2)
M_plotdat_x2$names <- factor(M_plotdat_x2$names, levels = M_plotdat_x2$names)
ggplot(M_plotdat_x2) + geom_bar(aes(x = names, y = as.numeric(paste(MSE_plot_x2))), stat = "identity")+xlab("Sampling and estimation method") +ylab("MSE of the X^2 coefficient") + coord_cartesian(ylim = c(0.014, 0.025))


