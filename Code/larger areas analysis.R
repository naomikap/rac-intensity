###################################################
#This is the code for estimation using larger areas
#The code here provides Figures 5, Table 1, testing the hypothesis that the lambda parameters are equal for all j using the random effects model conducted in Section 5.2
#In addition, it provides Figure 7 in section 6.2  
####################################################
# The r file "Organizing data to include subsets (larger areas)" was used to adjust the following files to include 14 subsets and 36 subsets

#1. locations_data.CSV: A data set of RAC locations. 
# We created a new file named locations_data_include_subset.CSV to include a column named sub_area which is the number of the subset the RAC belongs to (1 to 14)
# and the column sub_area2 which is the number of the subset the RAC belongs to in a different partition (1 to 36)
#in addition 3 variables are used in this analysis
#The shoe column which indicates the shoe number, 
#the x column which indicates the x axis of the RAC location 
#the y column which indicates the Y axis of the RAC location.

#2. contacts_data.txt: A data set of the contact surface 
#This is a pixel data where 1 indicates there is a contact surface and 0 otherwise
#There are 307 columns in each shoe and 395 is the number of rows
#The number of shoes is 387 but 386 is the number of shoes with RACs - shoe 127 has no RACS
#There are 307 columns in each shoe and 395 is the number of rows
# We created 2 files from this file chi and chi 2. 
#chi is the contact surface matrix which corresponds to 14 subsets and chi2 corresponds to 36 subsets
# Each row in these matrices corresponds to the shoe and the column indicates the subset 
#such that the cell indicates the number of pixels with contact surface in subset j of shoe i divided by 10000 in order for the number to be between zero and 1

# In addition from contacts_data.txt we create the list cont_use  and allcount matrix which are used to determine the contour of the shoe - this is only used to make figure 5 look nice
# See the "Organizing data to include subsets (larger areas).R" for details

################################################################################################
setwd("C:\\Users\\kapna\\Dropbox\\naomi-micha\\shoe_data\\Codes for JASA")
set.seed(312)

#install.packages("hglm")
#install.packages("fields")
#install.packages("car")
#install.packages("xtable")
#install.packages("ggplot2")
#install.packages("gridExtra")

library(hglm)
library(fields)
library(car)
library(xtable)
library(ggplot2)
library(gridExtra)

data_RAC<-read.csv("locations_data_include_subset.CSV",header=TRUE) 
chi<-as.matrix(read.csv("chi.CSV",header=TRUE))
chi2<-as.matrix(read.csv("chi2.CSV",header=TRUE))
#The  following files are used for Figure 5 below
cont_use<-list()
for(i in 1:14)
{
  cont_use[[i]]<-as.matrix(read.csv(paste("cont_use",i,".CSV" ,sep = "", collapse = NULL),header=TRUE))
  
}
allcont<-as.matrix(read.csv("allcont.CSV",header=TRUE))
#######################################################################################################
#The Naive estimator
# INPUT:
# ======
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
    #each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calcaulated from shoe and subset (if n_ij ==NULL) or can be given  
########################################################################################################

Naive_estimator<-function(shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij=NULL)
{
  if(is.null(n_ij))
  {
    num_shoes<-length(unique(shoe))
    num_areas<-length(unique(subset))
    n_ij <- matrix(as.numeric(table(shoe,subset)),num_shoes,num_areas,byrow=FALSE)
    
  }
  
  n_ij_chi_ij <- n_ij/contact_mat
  Naive <-colMeans(n_ij_chi_ij,na.rm=TRUE)
  return(Naive)
}
Naive<-Naive_estimator()

#######################################################################################################
#The random effects estimator
# INPUT:
# ======
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
  #each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calculated from shoe and subset (if n_ij ==NULL) or can be given  
########################################################################################################

Random_estimator<-function(shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij=NULL)
{ 
  num_shoes<-nrow(contact_mat)
  num_areas<-ncol(contact_mat)
  
  #Adjusting  the matrices for the random effect procedure
  #creating a chi matrix that its dimentions fits the n used in the hglm procedure
  #the rows of original chi matrix  are organized one after the other, creating a long vector
  
  chi_rand1<-contact_mat[1,]
  for (i in 2:num_shoes)
  {
    chi_rand1<-c(chi_rand1,contact_mat[i,])
  }
  chi_rand1<-as.vector(unlist(chi_rand1))
  chi_rand<-chi_rand1[chi_rand1>0] #We are not including observations where contact_mat=0 and thus contact_mat is adjusted
  
  #X - the design matrix - for the hglm procedure (the dimensions are M-num_shoes * J-num_areas)
  x<-numeric()
  for (i in 1:num_shoes)
  {
    x<-rbind(x,diag(num_areas))
  }
  
  x<-x[chi_rand1>0,] #Adjustind x to suit the correction of not including in the analysis observations (n) when chi=0
  # Z matrix (needed for hglm procedure)
  z<-NULL
  for (i in 1:num_shoes)
  {  
    z_temp <-matrix(0,num_areas,num_shoes)
    z_temp[,i]<-1
    z<-rbind(z,z_temp)
  }
  z<-z[chi_rand1>0,] #adjustind z to suit the correction of not including in the analysis observations (n) when chi=0
  
  if(is.null(n_ij))
  {
    n_ij <- matrix(as.numeric(table(shoe,subset)),num_shoes,num_areas,byrow=FALSE)
    
  }
  sum((contact_mat==0)*(n_ij)) #there are no RACs in places that contact_mat is 0 which is good
  #adjustind  n_ij for the random effect analysis - creating n
  n_tem<-0
  for (i in 1:num_shoes)
  {
    n_tem<-c(n_tem,n_ij[i,])  
  }
  n<-n_tem[-1]
  #We are not including in the random effect analysis observations that contact_mat=0
  n_rand<- n[is.na(n)==FALSE]
  n_rand<- n[(chi_rand1>0)]
  # applying Poisson model with Gamma distributed random effects 
  simul.pois<- hglm(y = n_rand, X = x, Z = z, family = poisson(link = log),rand.family = Gamma(link = log),offset=log(chi_rand),vcovmat=TRUE, method="EQL1")
  FixCoefMat<-summary(simul.pois)$FixCoefMat
  random<-exp(simul.pois$fixef)

  return(list(random=random,FixCoefMat=FixCoefMat,res=simul.pois))
}
Random<-Random_estimator()

#######################################################################################################
#The log CML function, is used by CML_estimator function below
# INPUT:
# ======
#l_lambda_minus1 log lambda not including the first lambda
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
    #each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calculated from shoe and subset (if n_ij ==NULL) or can be given  
########################################################################################################

funLogCMLlambda <- function(l_lambda_minus1,shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij=NULL)
{ 
  num_shoes<-nrow(contact_mat)
  num_areas<-ncol(contact_mat)
  if(is.null(n_ij))
  {
    n_ij <- matrix(as.numeric(table(shoe,subset)),num_shoes,num_areas,byrow=FALSE)
    
  }
  
  l_lambda<-c(0,l_lambda_minus1)
    ell<- -sum(n_ij*(
      log(
        contact_mat*(
          t(as.matrix(exp(l_lambda))%*%rep(1,nrow(contact_mat)))/(
            (contact_mat %*% exp(l_lambda))%*%rep(1,ncol(contact_mat))
          )
        )
      )
    ),na.rm = TRUE)


return(ell)
}  

#######################################################################################################
#The CML estimator
# INPUT:
# ======
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
    #each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calculated from shoe and subset (if n_ij ==NULL) or can be given  
########################################################################################################

CML_estimator<-function(shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij=NULL){   
  
  Naive_est<-Naive_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
    
  opt<-optim(log(Naive_est[-1]/Naive_est[1]), funLogCMLlambda,hessian = TRUE,shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)####
  theta<-  opt$par  # these are estimators of log(Lambdaj/Lambda1) j=2,...J
  hes<-opt$hessian
  LAM<-sum(Naive_est)

  l_lambda_cml_1<-log(LAM)-log(1+sum(exp(theta)))
  l_lambda_cml_2_j<-theta+log(LAM)-log(1+sum(exp(theta)))
  l_lambda_cml<-c(l_lambda_cml_1,l_lambda_cml_2_j)  
  cml<-exp(l_lambda_cml)
  return(list(CML=cml,hes=hes,theta=theta))
}
CML<-CML_estimator()

################### Figure 5: the 3 estimators intensities 
#since we used chi/10000 (see "Organizing data to include subsets (larger areas)".R),in the calculating estimators we should divide the estimators by 10000 (then they will be on the same scale as in the pixel analysis)
Naive_estimator_adj <- Naive/10000
random_effects_estimator_adj<-Random$rand/10000
CML_estimator_adj<-CML$CML/10000

#Adjustments for the plot (we use the contacts_data.txt data here to determine the contour of a shoe and create a nice plot)
n_sub<-14
est_adj<-list()
est_adj[[1]]<-CML_estimator_adj
est_adj[[2]]<-Naive_estimator_adj
est_adj[[3]]<-random_effects_estimator_adj
res<-list()
resu<-list()
for(j in 1:3)
{
  res[[j]]<-cont_use[[1]]*est_adj[[j]][1]
  for(i in 2:n_sub)
  {
    res[[j]] <- res[[j]]+ cont_use[[i]]*est_adj[[j]][i] #plot of all sub areas
    
  }
  
  res[[j]][allcont==0] <- NA
  res[[j]][res[[j]]==0] <- NA
  resu[[j]]<-res[[j]][40:355, 80:220]
  resu[[j]]<-resu[[j]][,ncol(resu[[j]]):1]
}
#The final Figure 5 with all 3 intensities
com_3_est<-cbind(resu[[2]],resu[[3]],resu[[1]])
image.plot(t(com_3_est[nrow(com_3_est):1,]),axes=FALSE,xlab='Naive                                          Random                                          CML')
pdf(file ="large_inten_JASA.pdf", height=6, width=6)
image.plot(t(com_3_est[nrow(com_3_est):1,]),axes=FALSE,xlab='Naive                     Random                     CML')
dev.off()

#######################################################################################################
#The Confidence Intervals calculated in Table 1
# INPUT:
# ======
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
    #each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calculated from shoe and subset (if n_ij ==NULL) or can be given  
########################################################################################################
CI<-function(shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij=NULL)
{ # shoe is a vector which indicates for each RAC, which shoe it belongs to
  # subset is a vector which indicates for each RAC, which subset it belongs to
  
  Naive<-Naive_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
  Rand_res<-Random_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
  Random<-Rand_res$random
  CML_res<-CML_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
  CML<-CML_res$cml
  #naive CI
  n_i_dot <-as.numeric(table(shoe))
  sum_j_chi_ij_lam_j <- contact_mat %*% Naive
  var_a<- mean(((n_i_dot^2-n_i_dot)/(sum_j_chi_ij_lam_j^2))-1) 
  
  one_chi_ij <- 1/contact_mat
  one_chi_ij[one_chi_ij==Inf] <- NA
  sum_one_chi_ij <- colSums(one_chi_ij,na.rm=TRUE)
  m_j <- colSums(contact_mat>0)
  var_naive <- var_a*Naive^2/m_j+Naive*sum_one_chi_ij/m_j^2
  var_log_naive<-var_naive/Naive^2
  Low_naive<-exp(log(Naive)-1.96*sqrt(var_log_naive))
  High_naive<-exp(log(Naive)+1.96*sqrt(var_log_naive))
  confidence_naive<-cbind(Low_naive,High_naive)
  confidence_naive_adj<-round(confidence_naive/10000,4)
  Conf_naive<-paste ("(",confidence_naive_adj[,1],",", confidence_naive_adj[,2],")" ,sep = "", collapse = NULL)
  
  #Random  CI
  coe <- Rand_res$FixCoefMat
  confidence_rand<-exp(cbind(coe[,1]-1.96*coe[,2],coe[,1]+1.96*coe[,2]))
  confidence_rand_adj<-round(confidence_rand/10000,4)
  Conf_rand<-paste ("(",confidence_rand_adj[,1],",", confidence_rand_adj[,2],")" ,sep = "", collapse = NULL)
  
  #CML CI
  num_areas<-length(unique(subset))
  hes <- CML_res$hes
  var_theta<-  solve(hes)
  #The reason that we do not have to multiply the hessian by -1 is that all of the evaluation has been done in terms of -1 times the log-likelihood.
  #This means that the hessian that is produced by optim is already multiplied by -1
  
  #applying the delta method to get log(Lambda_j) j=2,...,J variance 
  #creating the g matrix 
  theta<-CML_res$theta
  part1<- 1/(1+sum(exp(theta)))
  part2<-matrix(NA,num_areas,num_areas-1)
  for(j in 1:(num_areas-1))
  {
    part2[,j]<-rep(exp(theta[j]))
  }
  part3<- -part1*part2
  part4<-cbind(rep(0,num_areas),part3)
  diagi<-as.matrix(diag(x = 1, num_areas,num_areas))
  g_tag_t<- diagi+part4
  g_tag<-g_tag_t[,-1]
  var_1<-g_tag%*%var_theta%*%t(g_tag)
  d <- 1.96*sqrt(diag(var_1))
  
  Low_CML <-CML_res$CML*exp(-d)
  High_CML <-CML_res$CML*exp(d)
  
  confidence_CML<-cbind(Low_CML,High_CML)
  confidence_CML_adj<-round((confidence_CML/10000),4)
  Conf_CML<-paste ("(",confidence_CML_adj[,1],",", confidence_CML_adj[,2],")" ,sep = "", collapse = NULL)
  
  CI<-cbind(seq(1,num_areas),Conf_naive,Conf_rand,Conf_CML)
  colnames(CI)<-c("Sub area","Naive","Random","CML")
  return(xtable(CI))
  
}

CI()

#Testing the hypothesis that the lambda parameters are equal for all j using the random effects model conducted in Section 5.2 checking assumption of equal estimators
#######################################################################################################
#The likelihood of the random effects model in order to calculate the variance
# INPUT:
# ======
#params - includes sig^2, the variance of the random effects and the log of the lambda parameters
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
    #each row in corresponds to the shoe and the column indicates the subset 
########################################################################################################

funLogRanlambda <- function(params,shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi)
{
  num_shoes<-length(unique(shoe))
  num_areas<-length(unique(subset))
  n_ij <- matrix(as.numeric(table(shoe,subset)),num_shoes,num_areas,byrow=FALSE)
  
  sig2 <- params[1]
  l_lambda <- params[2:length(params)]
  return(-(
    sum(n_ij*log((t(as.matrix(exp(l_lambda))%*%rep(1,nrow(contact_mat)))*contact_mat)/(
      (sig2+(contact_mat %*% exp(l_lambda))%*%rep(1,ncol(contact_mat)))^(
        ((rowSums(n_ij)+sig2)/rowSums(n_ij))%*%t(rep(1,ncol(contact_mat)))
      )
    )
    )
    ,na.rm = TRUE)
    +
      nrow(n_ij)*(sig2*log(sig2)-log(gamma(sig2)))
    +
      sum(lgamma(rowSums(n_ij)+sig2))
  ) )
}

#######################################################################################################
#hypothsis_test function tests the hypothesis that the lambda parameters are equal for all j using the random effects model 
#conducted in Section 5.2 - checking assumption of equal estimators
# INPUT:
# ======
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
    #each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calculated from shoe and subset (if n_ij ==NULL) or can be given  

########################################################################################################
hypothsis_test<-function(shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij=NULL)
{
  num_areas<-length(unique(subset))
  Naive_est<-Naive_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
  res_rand<-Random_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
    
  n_i_dot <-as.numeric(table(shoe))
  sum_j_chi_ij_lam_j <- contact_mat %*% Naive_est
  var_a<- mean(((n_i_dot^2-n_i_dot)/(sum_j_chi_ij_lam_j^2))-1) 
  
  optrand<-optim(c(var_a,log(Naive_est)), funLogRanlambda, shoe=shoe,subset=subset,contact_mat=contact_mat,gr = NULL, method ="L-BFGS-B", lower = rep(0, num_areas), upper = rep(Inf, num_areas),hessian = TRUE)
  hes_rand <- optrand$hessian
  vc<-  solve(hes_rand)[-1,-1]
  matr <- cbind(rep(1,13),-diag(13))
  co<-optrand$par[-1]
  Hypoth_res<-linearHypothesis(res_rand$res,hypothesis.matrix=matr,rhs=rep(0,nrow(matr)),test=c("Chisq", "F"),vcov.=vc,coef.=co)
  
  return(Hypoth_res)
}

hypothsis_test()


#######################################################################################################
#Comparison of the three estimators 
#Figure 7 in section 6.2 
# INPUT:
# ======
# num_sim The number of simulations 
# shoe - is a vector which indicates for each RAC, which shoe it belongs to
# subset - is a vector which indicates for each RAC, which subset it belongs to
# contact_mat - is the contact surface matrix which  
#each row in corresponds to the shoe and the column indicates the subset 
#n_ij -  There are two options: can be calculated from shoe and subset (if n_ij ==NULL) or can be given  

########################################################################################################
estimateLambda<-function(num_sim=500,shoe=data_RAC$shoe,subset=data_RAC$sub_area,contact_mat=chi,n_ij =NULL)
{ 
  ptm <- proc.time()
  lambda_sim<-Naive_estimator(shoe=shoe,subset=subset,contact_mat=contact_mat,n_ij=n_ij)
  chi_sim<-contact_mat
  n_i_dot <-as.numeric(table(shoe))
  sum_j_chi_ij_lam_j <- contact_mat %*% lambda_sim
  var_a<- mean(((n_i_dot^2-n_i_dot)/(sum_j_chi_ij_lam_j^2))-1) 
 
  numberOfShoes<-nrow(chi_sim)
  numberOfAreasInShoe<-ncol(chi_sim)
  
  lambdasRandom <- matrix(NA,num_sim,numberOfAreasInShoe)
  lambdasCML <- matrix(NA,num_sim,numberOfAreasInShoe)
  lambdasNaive <- matrix(NA,num_sim,numberOfAreasInShoe)
  
  n<-matrix(NA,dim(chi_sim)[1],dim(chi_sim)[2])
  ##generating the number of RACs in pixel j of shoe i
  for (t in 1:num_sim)
  {
    print(t)
    a_sim <- rgamma(numberOfShoes,shape=1/var_a,scale = var_a)
    for (i in 1:numberOfShoes)
    {
      for (j in 1:numberOfAreasInShoe)
      { 
        if(chi_sim[i,j]>0) #we have RACs only where there is contact surface 
        {
          
          n[i,j]<- rpois(1,lambda_sim[j]*a_sim[i]*chi_sim[i,j])
          
        }
        
      }
    }
    
    lambdasRandom[t,]<-Random_estimator(shoe=NULL,subset=NULL,contact_mat=chi_sim,n_ij=n)$random
    lambdasCML[t,]<-CML_estimator(shoe=NULL,subset=NULL,contact_mat=chi_sim,n_ij=n)$CML
    lambdasNaive[t,]<-Naive_estimator(shoe=NULL,subset=NULL,contact_mat=chi_sim,n_ij=n)
    
  }
  
  # calculating bias and MSE
  meanRandom<-.colMeans(lambdasRandom, num_sim, numberOfAreasInShoe, na.rm = FALSE)
  meanNaive<-.colMeans(lambdasNaive, num_sim, numberOfAreasInShoe, na.rm = FALSE)
  meanCML<-.colMeans(lambdasCML, num_sim, numberOfAreasInShoe, na.rm = FALSE)
  
  VarRandom<- apply(lambdasRandom,2,var)
  VarCML<-apply(lambdasCML,2,var)
  VarNaive<-apply(lambdasNaive,2,var)
  
    mj<-colSums(chi_sim!=0)
  chi_sim_NA<-chi_sim
  chi_sim_NA[chi_sim==0] <- NA
  VarNaive_th<-(lambda_sim^2*var_a/mj)+(lambda_sim/(mj^2))* colSums(1/chi_sim_NA,na.rm = TRUE)
    
  
  biasRandom<-meanRandom-lambda_sim
  biasNaive<-meanNaive-lambda_sim
  biasCML<-meanCML-lambda_sim
  
  MSERandom<-(biasRandom)^2+VarRandom
  MSENaive<-(biasNaive)^2+VarNaive
  MSECML<-(biasCML)^2+VarCML
  
  meanMSERandom<-mean(MSERandom)
  meanMSENaive<-mean(MSENaive)
  meanMSECML<-mean(MSECML)
  
  meanbiasRandom<-mean(biasRandom)
  meanbiasNaive<-mean(biasNaive)
  meanbiasCML<-mean(biasCML)
  
  time<-as.vector(proc.time() - ptm)[3]
  return(list(lambda_sim=lambda_sim,meanRandom=meanRandom,meanNaive=meanNaive ,meanCML=meanCML,biasRandom=biasRandom,biasNaive=biasNaive,biasCML=biasCML,VarRandom=VarRandom,VarNaive=VarNaive,VarCML=VarCML,MSERandom=MSERandom,MSENaive=MSENaive,MSECML=MSECML,lambdasRandom=lambdasRandom,lambdasNaive=lambdasNaive,lambdasCML=lambdasCML, a_sim=a_sim , num_sim=num_sim, numberOfShoes=numberOfShoes ,numberOfAreasInShoe=numberOfAreasInShoe,meanMSERandom=meanMSERandom, meanMSENaive=meanMSENaive,meanMSECML=meanMSECML,meanbiasRandom=meanbiasRandom,meanbiasNaive=meanbiasNaive,meanbiasCML=meanbiasCML,VarNaive_th=VarNaive_th,time=time))
}
output1<-estimateLambda()

pdf("Figure7_bias_JASA.pdf", height=6, width=6)
bias_comparison<-rbind(output1$biasCML/output1$lambda_sim,output1$biasRandom/output1$lambda_sim,output1$biasNaive/output1$lambda_sim)
bias_comparison<-as.matrix(bias_comparison)
bias_comparison1<-rbind(bias_comparison[1,order(output1$lambda_sim)],bias_comparison[2,order(output1$lambda_sim)],bias_comparison[3,order(output1$lambda_sim)])
rownames(bias_comparison1) <- c( "CML", "Random","Naive")
colours <- c("blue", "green","yellow")
barplot(as.matrix(bias_comparison1), main="Relative bias", ylab = "Bias",xlab="Lambda values", cex.lab = 1.2, cex.main = 1.4, beside=TRUE, col=colours,names.arg=round(sort(output1$lambda_sim),1),cex.axis=0.8, cex.names=0.7)
legend("topright",  c( "CML", "Random","Naive"), cex=0.7, bty="n", fill=colours)
dev.off()

pdf("Figure7_MSE_JASA.pdf", height=6, width=6)
MSE_comparison<-rbind(output1$MSECML/output1$VarNaive_th,output1$MSERandom/output1$VarNaive_th,output1$MSENaive/output1$VarNaive_th)
MSE_comparison<-as.matrix(MSE_comparison)
MSE_comparison1<-rbind(MSE_comparison[1,order(output1$lambda_sim)],MSE_comparison[2,order(output1$lambda_sim)],MSE_comparison[3,order(output1$lambda_sim)])
rownames(MSE_comparison1) <- c( "CML", "Random","Naive")
print(MSE_comparison1)
colours <- c( "blue", "green","yellow")
out<-sort(output1$lambda_sim)
barplot(as.matrix(MSE_comparison1), main="MSE ratio", ylab = "MSE",xlab="Lambda values", cex.lab = 1.2, cex.main = 1.4, beside=TRUE, col=colours,names.arg=round(out,1),cex.axis=0.8, cex.names=0.7,ylim=c(0,1.35))
legend("topright",  c( "CML", "Random","Naive"), cex=0.7, bty="n", fill=colours)
abline(h =1 , untf = FALSE,col=1)
dev.off()


