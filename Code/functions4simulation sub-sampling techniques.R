############################################################################
# FUNCTIONS USED FOR THE SIMULATION STUDY - QUADRATIC LOGISTIC REGRESSION
############################################################################

###required packages
library(lme4)
library (survival)
library(matrixStats)
#####################


############################################################################
# Simple sampling 
# INPUT:
# ======
# dat is the data to be sampled from, it has 4 variables:
#   n_acc - binary response 1 for "yes" 0 for "no"
#   x, x_squared - covariates 
#   cluster - indicator for cluster
# ratio -  number of controls for each case to sample 
#         (it is assumed that the number of "0" is much larger than the number of "1") 
#
######################################################################
samp_meth1<-function(dat,ratio){
  n <- nrow(dat)
  total_case<-nrow(dat[dat$n_acc==1,])
  if (n < total_case*(1+ratio)) {stop("sample size larger than n")}
  samp.size <- ceiling((ratio+1)*total_case)
  samp.obs <- sample(n, size=samp.size, replace=FALSE)
  data_samp<- dat[samp.obs,]
  return(data_samp)   
}

############################################################################
# Case control sampling ignoring clusters
# all cases are sampled and controls are sampled from the entire data 
# INPUT:
# ======
# dat is the data to be sampled from, it has 4 variables:
#   n_acc - binary response 1 for "yes" 0 for "no"
#   x, x_squared - covariates 
#   cluster - indicator for cluster
# ratio -  number of controls for each case to sample 
#         (it is assumed that the number of "0" is much larger than the number of "1") 
#
######################################################################
samp_meth2<-function(dat,ratio){
  n <- nrow(dat)
  case_ind <- dat$n_acc==1
  case <- dat[case_ind,]
  ctrl <- dat[!(case_ind),]
  total_case<-nrow(case)
  total_control<-n-total_case
  if (n < total_case*(1+ratio)) {stop("sample size larger than n")}
  samp.ctrl.size <- ceiling(ratio*total_case)
  samp.ctrl <- sample(total_control, size=samp.ctrl.size, replace=FALSE)
  control <- ctrl[samp.ctrl,] 
  # offset - same for all
  off <- rep(-log(samp.ctrl.size/total_control),(total_case+samp.ctrl.size))  
  data_samp <- data.frame(cbind(rbind(case,control),off))
  return(data_samp)     
}


############################################################################
# Case control sampling proportional to cluster size
# all cases  + sample of controls from each cluster proportional to cluster size
# INPUT:
# ======
# dat is the data to be sampled from, it has 4 variables:
#   n_acc - binary response 1 for "yes" 0 for "no"
#   x, x_squared - covariates 
#   cluster - indicator for cluster
# ratio -  number of controls for each case to sample 
#         (it is assumed that the number of "0" is much larger than the number of "1") 
#
######################################################################
samp_meth3<-function(dat,ratio){
  n <- nrow(dat)
  clst.size <- table(dat$cluster)  # total cluster sizes
  clst.id <- as.integer(names(clst.size))
  n.clst <- length(clst.size)  # no. of clusters
  case <- dat[dat$n_acc==1,]
  ctrl <- dat[dat$n_acc==0,]
  total_case <- nrow(case)
  total_control <- nrow(ctrl)
  samp.ctrl.size <- ceiling(ratio*total_case*clst.size/n) # no. to samp from each cluster
  samp <- numeric()
  for(j in 1:n.clst) { # sample controls from each cluster
    id <- clst.id[j]
    ctrl.j <- ctrl[ctrl$cluster==id,]
    case.j <- case[case$cluster==id,]
    if(dim(case.j)[1]>0) { # remove clusters with 0 cases
      size.j <- dim(ctrl.j)[1]
      if (size.j < samp.ctrl.size[j]){
        warning(paste("cluster",id,"has fewer controls than should be sampled - all controls are used"))
      }
      n.samp <- min(size.j,samp.ctrl.size[j])
      s.ctrl <- sample(x = size.j,size = n.samp,replace = FALSE)
      off <- -log(n.samp/size.j)
      case.j$off <- off
      ctrl.j$off <- off
      samp <- rbind(samp,case.j)
      samp <- rbind(samp,ctrl.j[s.ctrl,])
    }
  }
  data_samp <- data.frame(samp)
  return(data_samp)    
}

############################################################################
# Case control sampling proportional to number of cluster-specific cases
# all cases  + sample of controls from each cluster proportional to the number
# of cases in the cluster
# INPUT:
# ======
# dat is the data to be sampled from, it has 4 variables:
#   n_acc - binary response 1 for "yes" 0 for "no"
#   x, x_squared - covariates 
#   cluster - indicator for cluster
# ratio -  number of controls for each case to sample 
#         (it is assumed that the number of "0" is much larger than the number of "1") 
#
######################################################################
samp_meth4<-function(dat,ratio){
  case <- dat[dat$n_acc==1,]
  ctrl <- dat[dat$n_acc==0,]
  total_case <- table(case$cluster) #number of cases in each cluster
  clst.id <- as.integer(names(total_case))
  n.clst <- length(clst.id)  # no. of clusters
  samp.ctrl.size <- ceiling(ratio*total_case) # no. to samp from each cluster
  samp <- numeric()
  for(j in 1:n.clst) { # sample controls from each cluster
    id <- clst.id[j]
    ctrl.j <- ctrl[ctrl$cluster==id,]
    case.j <- case[case$cluster==id,]
    if(dim(case.j)[1]>0) { # remove clusters with 0 cases 
      size.j <- dim(ctrl.j)[1]
      if (size.j < samp.ctrl.size[j]){
        warning(paste("cluster",id,"has fewer controls than should be sampled - all controls are used"))
      }
      n.samp <- min(size.j,samp.ctrl.size[j])
      s.ctrl <- sample(x = size.j,size = n.samp,replace = FALSE)
      off <- -log(n.samp/size.j)
      case.j$off <- off
      ctrl.j$off <- off
      samp <- rbind(samp,case.j)
      samp <- rbind(samp,ctrl.j[s.ctrl,])
    }
  }
  data_samp <- data.frame(samp)
  return(data_samp)    
}


########################################################
# Calling the sampling method.
# samp_methods - take values 1-4 (see methods above)
# other input variables - see methods above
#########################################################

sampling <-function(samp_method,dat,ratio) {
  if  (samp_method==0){data_samp<-dat}
  if  (samp_method==1){data_samp<-samp_meth1(dat,ratio)}
  if  (samp_method==2){data_samp<-samp_meth2(dat,ratio)}
  if  (samp_method==3){data_samp<-samp_meth3(dat,ratio)}
  if  (samp_method==4){data_samp<-samp_meth4(dat,ratio)}
  return(data_samp)    
}

############################################################################
# Estimation either by random effects or conditional logistic regression
# 
# INPUT:
# ======
# da: data to analyze with 4 or 5 variables:
#   n_acc - binary response 1 for "yes" 0 for "no"
#   x, x_squared - covariates 
#   cluster - indicator for cluster
#   off - offset for case-control sampling requiring an offset
# es_method: 1 - random effects using glmer or 2 - cml using clogit 
# samp_meth: the sampling mechanisim generated the data.
#            one of the following options should be used
#                 0 no sampling - original data
#                 1 is samp_meth1 - random sampling
#                 2 is samp_meth2 - case control ignoring clusters
#                 3 is samp_meth3 - case control prop to cluster size
#                 4 is samp_meth4 - case control prop to no. cluster cases
#Intercept:
#In logistic regression, the estimated effect is consistent when using case
#control sampling, but the estimated intercept may not be valid (Agresti, 2013).
        #FALSE-If we are not interested  in the intercept, we do not need to include an offset
        #True- The intercept is adjusted using offset terms
######################################################################
est<-function(da,es_method,samp_meth,intercept=FALSE){
  if (es_method==1){ 
    if (samp_meth==0 | samp_meth==1){
      rand <- glmer(n_acc ~ x + x_square + (1 | cluster), data=da, 
                    family=binomial(link="logit"),
                    control = glmerControl(optimizer = "bobyqa"))
      estimate <- c(rand@beta,rand@theta)
    }
    if (samp_meth==2 | samp_meth==3 | samp_meth==4){
      rand <- glmer(n_acc ~ offset(off) + x + x_square + (1 | cluster), data=da, 
                    family=binomial(link="logit"),
                    control = glmerControl(optimizer = "bobyqa"))
      estimate <- c(rand@beta,rand@theta)
    }
  }
  if(es_method==2){
    cml <- clogit(n_acc ~ x + x_square + strata(cluster), data=da)
    estimate <- coef(cml)
    
    if(intercept==TRUE)
    {
      xx <- cbind(da$x,da$x_square)
      off <- xx%*%estimate
      rand_1 <- glmer(da$n_acc ~ offset(off) + (1 | da$cluster), 
                      family=binomial(link="logit"),
                      control = glmerControl(optimizer = "bobyqa"))
      estimate_temp_2 <- c(rand_1@beta,rand_1@theta)
      estimate <- c(estimate_temp_2[1],estimate[1],estimate[2],estimate_temp_2[2])
    }
    
  }
  return(estimate)  
}

#############################################################################
# Generate data for simulation
# model is random effect logistic regression with x and x^2 as covariates
# and cluster-specific random effect having a zero-mean normal distribution
#
# input:
###########
# num_cluster - numbr of clusters
# cluster_size - mean cluster size
# sigma_real - std of random effect (random effect ~ N(0,sigma_real))
# beta_0_r - intercept of logistic regression
# beta_1_r - coefficient of x
# beta_2_r - coefficient of x^2
#
# Output:
###########
# data set with 4 variables: 
# cluster - indicator of cluster
# x - a deterministic sequence from 1/cluster_size to 1 (location of pixel)
# x_square - x^2 
# n_acc - the binary 0/1 variable (response)
################################################################

create_dat <- function(num_cluster=500,cluster_size=500,sigma_real=0.75,beta_0_r=-3,beta_1_r=2,beta_2_r=-2){
  a <- rnorm(num_cluster,0,sigma_real) 
  cluster <- as.factor(rep(1:num_cluster,each=cluster_size))  
  x <- rep(1:cluster_size,times=num_cluster)/cluster_size
  x_square <- x^2
  prob <- exp(beta_0_r+beta_1_r*x+beta_2_r*x_square+a[cluster])/(1+exp(beta_0_r+beta_1_r*x+beta_2_r*x_square+a[cluster]))
  n_acc <- (runif(num_cluster*cluster_size) < prob )*1
  dat <- data.frame(cbind(n_acc, x, x_square, cluster))
  return(dat)
}


