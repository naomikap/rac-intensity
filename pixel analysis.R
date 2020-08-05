#This is the code for estimation using maximum resolution (pixel analysis)
#The code  here produces Figures 2,3 and the global test of all coefficients equal to zero using the random effects model conducted in 4.5.1
#In addition, the code produces Figures 1-4 in the Web appendix 
################################################################################################
#Two data sets are used here:
#1. locations_data.CSV: A data set of RAC locations. 
    #The first three columns are used. 
    #The first column indicates the shoe number, 
    #the second indicates the x axis of the RAC location 
    #the third indicates the Y axis of the RAC location. 
#2. contacts_data.txt: A data set of the contact surface 
    #This is a pixel data where 1 indicates there is a contact surface and 0 otherwise
    #There are 307 columns in each shoe and 395 is the number of rows
    #The number of shoes is 387 but 386 is the number of shoes with RACs - shoe 127 has no RACS
    
################################################################################################
setwd("C:\\Users\\kapna\\Dropbox\\naomi-micha\\shoe_data\\Codes for JASA")
###################
set.seed(312)

#install.packages("splines")
#install.packages("lme4")
#install.packages("rgl")
#install.packages("fields")
#install.packages("survival")
#install.packages("smoothie")
#install.packages("car")
#install.packages("ggplot2") 
#install.packages("reshape2")

library(splines)
library(lme4)
library(rgl)
library(fields)
library(survival)
library(smoothie)
library(car)
library(ggplot2) 
library(reshape2)
##################

col_shoe<-307 #307 is the number of columns in each shoe
row_shoe<-395 #395 is the number of rows in each shoe
num_shoe<-387  #387 is the number of shoes but 386 is the number of shoes with RACs - shoe 127 has no RACS
rel_col_shoe<-150  #out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
rel_row_shoe<-300  #out of the 395 rows only 300 are relevant (contain non zero pixels in some shoes)
rel_x_cord<-0.25 #using coordinates as in the locations_data.CSV file the relevant x coordinates are between -.25 and 0.25
rel_Y_cord<-0.5 #the relevant Y coordinates are between -0.5 and 0.5



#The following two functions converts the x and Y coordinate of the location of a RAC to the X and Y pixels

##################################################################################################################
# aspix_x converts the x coordinate to the x pixel
# INPUT:
# ======
# x - the x coordinate
# col_shoe - the number of columns in each shoe
# rel_col_shoe -the number of relevant columns 
    #out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
# rel_x_cord - the relevant coordintes 
    #(using coordinates as in the locations_data.CSV file. The relevant x coordinates are between -.25 and 0.25)
##################################################################################################################
aspix_x <-function(x,col_shoe=307,rel_col_shoe=150,rel_x_cord=0.25)
{
  not_rel_col<-ceiling((col_shoe - rel_col_shoe)/2)
  delx <- (2*rel_x_cord)/rel_col_shoe 
  pix_x <- col_shoe-(floor((x+rel_x_cord)/delx)+not_rel_col) #The plus rel_x_cord is because it is --rel_x_cord (the x starts from -rel_x_cord) 
  return(pix_x)    
}

##################################################################################################################
# aspix_y converts the Y coordinate to the Y pixel
# INPUT:
# ======
# y - the y coordinate
# row_shoe - the number of rows in each shoe
# rel_row_shoe -the number of relevant rows 
    #out of the 395 rows only 300 are relevant (contain non zero pixels in some shoes)
# rel_Y_cord - the relevant coordintes 
    #(using coordinates as in the locations_data.CSV file. the relevant Y coordinates are between -0.5 and 0.5)
##################################################################################################################

aspix_y<-function(y,row_shoe=395,rel_row_shoe=300,rel_Y_cord=0.5)
{
  not_rel_row<-ceiling((row_shoe-rel_row_shoe)/2)
  dely <- (2*rel_Y_cord)/rel_row_shoe 
  pix_y <- row_shoe-(floor((y+rel_Y_cord)/dely)+not_rel_row) # The plus rel_Y_cord is because it is --rel_Y_cord (the y starts from -0.5)
  return(pix_y) 
}  


#############################
#organizing the contacts_data 
#############################

#We are importing the contacts_data as character and creating a list of contact shoe matrices 
d <- readChar("contacts_data.txt",nchars=(col_shoe*row_shoe+2)*num_shoe) 
data <- list() 
for(i in 1:num_shoe)
{
  data[[i]] <- matrix(as.numeric(unlist(strsplit(substr(d, 1+(col_shoe*row_shoe+2)*(i-1), (col_shoe*row_shoe+2)*i-2), split="")) ),row_shoe,col_shoe,byrow=1)
}

#Shoe 9 should be mirrored as all other shoes
shoe9rev <- data[[9]] #(compare image(data[[8]]) and image(data[[9]]))
data[[9]] <- data[[9]][,ncol(data[[9]]):1]

#########cleaning the data set###########################################
#There are identifying stamps the police put on each shoeprint
#These are not part of the shoe's contact surface and thus are omitted  
#The first stage in cleaning the stamps was to try to separate them from the actual contact surface 
#We try to find the lower bound of the cumulative contact surfce to separate the stamps from the actual contact surface
#we found that if we look only at contact surface that appeared in more than 8 shoes it provided a relatively good separation 
allcont <- data[[1]]
for(i in 2:num_shoe)
{
  allcont <- allcont+data[[i]] #this is the contact of all shoes  
}
 
allcont <- (allcont>=8)*1 #here we see pixels that appear in more than 8 shoes

#Removing the stamps
#finding the lower bound of the contact surface
h_width<-floor(row_shoe/2) #this is half the width of the shoe
lb<- rep(NA,h_width) 
j<-1
while(allcont[h_width,j]==0) j<-j+1
lb[1] <- j-1 
for(i in 2:h_width) {
  j<- lb[i-1]
  if(allcont[h_width-i+1,j]==0) {
    while((allcont[h_width-i+1,j]==0)&&j<rel_row_shoe) j <- j+1
    lb[i] <- j-1
  }else{
    while((allcont[h_width-i+1,j]==1)&&j>0) j <- j-1
    lb[i] <- j
  }
}

for(i in 1:h_width) allcont[h_width-i+1,1:lb[i]] <- 0 #removing the lower stamp
#the upper bound of the contact surface
ub<- rep(NA,h_width) 

j<-col_shoe
while(allcont[h_width,j]==0) j<-j-1
ub[1] <- j+1 
for(i in 2:h_width) {
  j<- ub[i-1]
  if(allcont[h_width-i+1,j]==0) {
    while((allcont[h_width-i+1,j]==0)&&j>0) j <- j-1
    ub[i] <- j+1
  }else{
    while((allcont[h_width-i+1,j]==1)&&j<rel_row_shoe) j <- j+1
    ub[i] <- j
  }
}

for(i in 1:h_width) allcont[h_width-i+1,ub[i]:col_shoe] <- 0 #removing the upper stamp

for(i in 1:num_shoe) {
  data[[i]] <- data[[i]]*allcont
}

###################Working with the locations data##############
acciden<-read.csv("locations_data.CSV",header=TRUE)
acci <- list()
delx <- 2*rel_x_cord/rel_col_shoe 
dely <- 2*rel_Y_cord/rel_row_shoe 

for (i in (c(1:126,128:num_shoe)) )#shoe 127 doesn't have RACs
  {
  acci[[i]] <- matrix(0,row_shoe,col_shoe) 
  locations <- cbind(acciden$x[acciden$shoe==i],acciden$y[acciden$shoe==i]) # the coordinates of the RAC
  for(j in 1:nrow(locations)) {
    xpix <- aspix_x(locations[j,1])
    ypix<-aspix_y(locations[j,2])
    acci[[i]][ypix,xpix] <- acci[[i]][ypix,xpix]+1 #if there is more than one RAC (accidental) in a pixel we will count it as well
  }
}

###RACs can be observed only on the contact surface, but as we show below, the data has RACs where there is no contact surface
m <- rep(NA,num_shoe)
for(i in (c(1:126,128:num_shoe)))
{
m[i] <- min(data[[i]][acci[[i]]>=1]) # checking to see if there are RACs where there is no contact surface
}
# 0 means that there is at least one RAC that is not on the contact surface

# As noted in Section 4. When RACs are created they may tear the shoe sole such that the location of the RAC appears to be on an area with 
#no contact surface and thus the value of the contact surface is set to 1 in all cases where there are RACs
data_temp <- list() # a "solution", add contact surface where there is a RAC.  
for(i in (c(1:126,128:387))) {
data_temp[[i]] <- data[[i]]
data_temp[[i]][acci[[i]]>=1] <- 1
}

data_pix<-list() 
# each data_pix[[i]] is a matrix with column 1 indicating the shoe, 2 the x, 3 the y, 4 the amount of RACs in that pixel
# we include only data where there is contact surface (after adjusting for the case that if there is a RAC there will be contact surface)
for(i in (c(1:126,128:num_shoe))) 
  {
  xcoor <- t(matrix(rep((-col_shoe/2+1:col_shoe)*rel_Y_cord/rel_col_shoe,row_shoe),col_shoe,row_shoe))
  ycoor <- -matrix(rep((-row_shoe/2+1:row_shoe)*rel_Y_cord/rel_col_shoe,col_shoe),row_shoe,col_shoe)
  shoe<-rep(i,length(data_temp[[i]][data_temp[[i]]==1]))
  data_pix[[i]]<-cbind(shoe,xcoor[data_temp[[i]]==1],ycoor[data_temp[[i]]==1],acci[[i]][data_temp[[i]]==1])# the data is only  where there is contact surface
  }

data_pix_use<-numeric()
for (i in (c(1:126,128:num_shoe)))
{
  data_pix_use<-rbind(data_pix_use,data_pix[[i]]) 
}

#As noted in Section 4 of the article, the number of RACS is set to 1 in 38 cases where there are 2 RACs. 
#Appearance of two RACs in the same pixel may be due to the way the data were pre-processed and the location was defined.
n_Acc<-data_pix_use[,4]
#data_pix_use[n_Acc==2,] -> These are the 38 pixels with 2 RACs
n_Acc[n_Acc>=1] <-1 # more than one RAC in a shoe is considered as 1
x<- data_pix_use[,2]
y<- data_pix_use[,3]
shoe<-as.factor(data_pix_use[,1]) #it should be noted that as factor changes the numbering
                                  #since shoe 127 doesnt exist, as factor makes the numbering of shoes 128 to 387 to decrease by 1. (shoe 128 is now 127 etc.)
mydata <- data.frame(cbind(n_Acc, x, y,shoe))  #This is the data that will be used                             
for(j in 1:nrow(locations)) {
     xpix <-aspix_x(locations[j,1])   
     ypix <- aspix_y(locations[j,2])  
     acci[[i]][ypix,xpix] <- acci[[i]][ypix,xpix]+1 #if there is more than one RAC in a pixel we will count it as well
   }
sumacci <- acci[[1]]
for(i in c(2:126,128:387))
{
  sumacci <- sumacci+acci[[i]]
}

sumcont <- data[[1]]
for(i in c(2:126,128:387)) 
{
  sumcont <- sumcont+data[[i]]
  
}
###############creating case control data##########################################################################
# As noted in Section 4.4, estimating the intensity function at a high resolution is computationally challenging 
    #and thus case-control sub-sampling techniques are used
#The calculations were based on within-cluster case-control sub-sampling, 
#which includes all cases (pixels with RACs, nij = 1) and 20 random controls (pixels without RACs, nij = 0) from each shoe
        dataCC <- numeric()
        for(i in 1:length(unique(shoe)))
          {
            case <- mydata[mydata$shoe==i&mydata$n_Acc>0,]
            control <- mydata[mydata$shoe==i&mydata$n_Acc==0,]
            control <- control[sample(nrow(control),size=20,replace=FALSE),]
            dataCC <- rbind(dataCC,case,control)
          }
##################################################################################################################

##################################################################################################################
# The naive smooth estimator used on the basis of the entire data, not the case control
#A uniform kernel is used 
  #where each entry of the smoothed matrix is calculated as the average of its 21^2 neighbor entries in the original matrix.
# INPUT:
# ======
#cumRAC is the comulative matrix of RAC locations of all shoes
#cumContact is the comulative matrix of all contact surfaces of all shoes
#areaShoe is the area of the of the shoes which defines the contour of all shoes
  #In our case is all pixels that appear in more than 8 shoes
##################################################################################################################

Naive<-function(cumRAC=sumacci,cumContact=sumcont,areaShoe=allcont) 
{
  
  Naivemat<- cumRAC/cumContact
  Naivemat[areaShoe==0] <- NA
  est <- kernel2dsmooth(Naivemat, kernel.type="boxcar",n=21)
  est[areaShoe==0] <- NA
  return(est)
}

naive_smooth<-Naive()
image.plot(naive_smooth,axes=FALSE)

#The random effects and the CML estimates were calculated using a product of natural cubic splines
#Three knots for the X-axis and five knots for the Y-axis were used and their positions were set according to equal quantiles. 
#These numbers of knots enabled flexibility and still avoided computational problems

##################################################################################################################
# The random effects estimator 
# INPUT:
# ======
#nknotsx the number of x knots using a product of natural cubic splines
#nknotsy the number of y knots using a product of natural cubic splines
#dat is the data used for estimation, we are using here the case control data
##################################################################################################################
####The random effects estimator
Random<-function(nknotsx=3,nknotsy=5,dat=dataCC)
{
  knotsx <- as.numeric(quantile(dat$x,1:nknotsx/(1+nknotsx)))
  knotsy <-as.numeric(quantile(dat$y,1:nknotsy/(1+nknotsy)))
  shoe<-dat$shoe
  est<- glmer(dat$n_Acc ~ ns(dat$x,knots=knotsx):ns(dat$y,knots=knotsy)+(1 | shoe) , data= dat , family=binomial(link="logit"),control = glmerControl(optimizer = "bobyqa"))
  return(est)
}
rand<-Random()

#plot of the random effects estimator
nknotsx <- 3
nknotsy <- 5
knotsx <- as.numeric(quantile(dataCC$x,1:nknotsx/(1+nknotsx)))
knotsy <-as.numeric(quantile(dataCC$y,1:nknotsy/(1+nknotsy)))
basx <- ns(dataCC$x,knots=knotsx)
basy <- ns(dataCC$y,knots=knotsy)
xy <- expand.grid(xcoor[1,],ycoor[,1]) 
newdesignmat <- rep(1,length(xy[,1]))
for(i in 1:length(predict(basy,1))) {
  for(j in 1:length(predict(basx,1))) {
    newdesignmat <-  cbind(newdesignmat,predict(basx, xy[,1])[,j]*predict(basy, xy[,2])[,i])
  }
}
pred.case_control <- newdesignmat%*%fixef(rand)+log(0.005) #log(0.005) is the offset
pred.case_control[t(allcont)==0] <- NA #areas out of the contour (less than 8 shoes has contact surface in these pixels) are given NA
prob.pred <- exp(matrix(pred.case_control ,row_shoe,col_shoe,byrow=1))/(1+exp(matrix(pred.case_control ,row_shoe,col_shoe,byrow=1)))
intens <- -log(1-prob.pred) #turning it to intensity
image.plot(intens,axes=FALSE)
m<-mean(pred.case_control,na.rm=TRUE)  #for use in the CML

##################################################################################################################
# The CML estimator 
# INPUT:
# ======
#nknotsx the number of x knots using a product of natural cubic splines
#nknotsy the number of y knots using a product of natural cubic splines
#dat is the data used for estimation, we are using here the case control data
##################################################################################################################
CML<-function(nknotsx=3,nknotsy=5,dat=dataCC)
{
  knotsx <- as.numeric(quantile(dat$x,1:nknotsx/(1+nknotsx)))
  knotsy <-as.numeric(quantile(dat$y,1:nknotsy/(1+nknotsy)))
  shoe<-dat$shoe
  est<- clogit(dat$n_Acc~  ns(dat$x,knots=knotsx):ns(dat$y,knots=knotsy)+strata(shoe) , data=dat)
  return(est)
}
cml<-CML()
#plot cml
newdesignmat1 <- rep(1,length(xy[,1]))
for(i in 1:length(predict(basy,1))) {
  for(j in 1:length(predict(basx,1))) {
    newdesignmat1 <-  cbind(newdesignmat1,predict(basx, xy[,1])[,j]*predict(basy, xy[,2])[,i])
  }
}

pred.cml.bin.case_control <- newdesignmat1%*%c(0,coefficients(cml)) # the intercept cant be estimated since it cancels 
pred.cml.bin.case_control[t(allcont)==0] <- NA #areas out of the contour (less than 8 shoes has contact surface in these pixels) are given NA
m_1<-mean(pred.cml.bin.case_control,na.rm=TRUE)
pred.cml.bin.case_control<-pred.cml.bin.case_control-m_1+m #making the means of randon and cml to be equal
prob.pred_cml <- exp(matrix(pred.cml.bin.case_control ,row_shoe,col_shoe,byrow=1))/(1+exp(matrix(pred.cml.bin.case_control ,row_shoe,col_shoe,byrow=1)))
intens.pred_cml <- -log(1-prob.pred_cml)
image.plot(intens.pred_cml,axes=FALSE) #notice that these probabilities depend on the intercept which is not included since it cancels. 

# Figure 2: the 3 estimators intensities on the same scale 
sub <- 70
cols <- sub:(col_shoe-sub)
#we multiply CML and random so they will be on the same scale 
com_3_est<-cbind(naive_smooth[,cols],exp(-0.9915/2)*intens[,cols],exp(-0.9915/2)*intens.pred_cml[,cols]) #0.9915 is sigma^2 of the random effect. e^(sigma^2/2) is the expectation of a log linear variable lognormal(0,sigma^2). This is the expectation of the random.
image.plot(t(com_3_est[nrow(com_3_est):1,]),axes=FALSE,xlab='Naive,Random,CML')
pdf(file ="pixel_inten_JASA.pdf", height=6, width=6)
image.plot(t(com_3_est[nrow(com_3_est):1,]),axes=FALSE,xlab='Naive                     Random                     CML')
dev.off()

## hypothesis testing## (Section 4.5.1)
co <- fixef(rand)
vc <- vcov(rand)
matr <- diag(length(co))[-1,]
testing<-linearHypothesis(rand,hypothesis.matrix=matr,rhs=rep(0,length(co)-1),test=c("Chisq", "F"),vcov.=vc,coef.=co)
testing$`Pr(>Chisq)` #pvalue is approximately zero

#confidence intervals - Figure 3
##################################################################################################################
# CML confidence interval
# INPUT:
# ======
#dat - the data used for estimation, we are using here the case control data
# col_shoe - the number of columns in each shoe
# row_shoe - the number of rows in each shoe
# rel_Y_cord - the relevant coordintes 
  #(using coordinates as in the locations_data.CSV file. the relevant Y coordinates are between -0.5 and 0.5)
# rel_col_shoe -the number of relevant columns 
  #out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
#nknotsx the number of x knots using a product of natural cubic splines
#nknotsy the number of y knots using a product of natural cubic splines
##################################################################################################################
CI_cml <- function(dat=dataCC,col_shoe=307,row_shoe=395,rel_Y_cord=0.5,rel_col_shoe=150,nknotsx=3,nknotsy=5) 
{
  xcoor <- t(matrix(rep((-col_shoe/2+1:col_shoe)*rel_Y_cord/rel_col_shoe,row_shoe),col_shoe,row_shoe))
  ycoor <- -matrix(rep((-row_shoe/2+1:row_shoe)*rel_Y_cord/rel_col_shoe,col_shoe),row_shoe,col_shoe)
  cml_bin_fit<-CML(dat=dat,nknotsx=nknotsx,nknotsy=nknotsy)
  rand<-Random(dat=dat,nknotsx=nknotsx,nknotsy=nknotsy)
 
  knotsx <- as.numeric(quantile(dat$x,1:nknotsx/(1+nknotsx)))
  knotsy <-as.numeric(quantile(dat$y,1:nknotsy/(1+nknotsy)))
  basx <- ns(dat$x,knots=knotsx)
  basy <- ns(dat$y,knots=knotsy)
  xy <- expand.grid(xcoor[1,],ycoor[,1]) 
  newdesignmat <- rep(1,length(xy[,1]))
  for(i in 1:length(predict(basy,1))) {
    for(j in 1:length(predict(basx,1))) {
      newdesignmat <-  cbind(newdesignmat,predict(basx, xy[,1])[,j]*predict(basy, xy[,2])[,i])
    }
  }
  pred.case_control_r <- newdesignmat%*%fixef(rand)+log(0.005) #log(0.005) is the offset
  pred.case_control_r[t(allcont)==0] <- NA #areas out of the contour (less than 8 shoes has contact surface in these pixels) are given NA
  prob.pred <- exp(matrix(pred.case_control_r ,row_shoe,col_shoe,byrow=1))/(1+exp(matrix(pred.case_control_r ,row_shoe,col_shoe,byrow=1)))
  m<-mean(pred.case_control_r,na.rm=TRUE) 

  cov_mat<-matrix(as.numeric(vcov(cml_bin_fit)),length(coefficients(cml_bin_fit)))  
  avg <- newdesignmat%*%c(0,coefficients(cml_bin_fit))
  avg[t(allcont)==0] <- NA #areas out of the contour (less than 8 shoes has contact surface in these pixels) are given NA
  m_1 <- mean(avg,,na.rm=TRUE)
  avg<- avg-m_1+m
  
  newdesignmat<-newdesignmat[,-1]
  std <- as.matrix(sqrt(rowSums((newdesignmat%*%cov_mat)*newdesignmat)))
  ex <- exp(avg%*%c(1,1,1)+1.96*std%*%c(-1,0,1))
  pr <- ex/(1+ex)
  inte <- -log(1-pr)
  for(i in 1:3) inte[as.vector(t(allcont)==0),i] <- NA
  ret <- list(matrix(inte[,1],row_shoe,col_shoe,byrow=1),matrix(inte[,2],row_shoe,col_shoe,byrow=1),matrix(inte[,3],row_shoe,col_shoe,byrow=1))
  names(ret) <- c("Low","Mid","High")
  return(ret)
}

CI <- CI_cml()

##################################################################################################################
# Random confidence interval
# INPUT:
# ======
#dat - the data used for estimation, we are using here the case control data
# col_shoe - the number of columns in each shoe
# row_shoe - the number of rows in each shoe
# rel_Y_cord - the relevant coordintes 
  #(using coordinates as in the locations_data.CSV file. the relevant Y coordinates are between -0.5 and 0.5)
# rel_col_shoe -the number of relevant columns 
  #out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
#nknotsx the number of x knots using a product of natural cubic splines
#nknotsy the number of y knots using a product of natural cubic splines
##################################################################################################################
CI_random <- function(dat=dataCC,col_shoe=307,row_shoe=395,rel_Y_cord=0.5,rel_col_shoe=150,nknotsx=3,nknotsy=5) 
{
  xcoor <- t(matrix(rep((-col_shoe/2+1:col_shoe)*rel_Y_cord/rel_col_shoe,row_shoe),col_shoe,row_shoe))
  ycoor <- -matrix(rep((-row_shoe/2+1:row_shoe)*rel_Y_cord/rel_col_shoe,col_shoe),row_shoe,col_shoe)
  
  fit.rand<-Random(dat=dat,nknotsx=nknotsx,nknotsy=nknotsy)
  cov_mat<-matrix(as.numeric(vcov(fit.rand)),length(fixef(fit.rand))) 
  knotsx <- as.numeric(quantile(dat$x,1:nknotsx/(1+nknotsx)))
  knotsy <-as.numeric(quantile(dat$y,1:nknotsy/(1+nknotsy)))
  basx <- ns(dat$x,knots=knotsx)
  basy <- ns(dat$y,knots=knotsy)
  xy <- expand.grid(xcoor[1,],ycoor[,1]) 
  newdesignmat <- rep(1,length(xy[,1]))
  
  for(i in 1:length(predict(basy,1))) 
  {
    for(j in 1:length(predict(basx,1))) 
    {
      newdesignmat <-  cbind(newdesignmat,predict(basx, xy[,1])[,j]*predict(basy, xy[,2])[,i])
    }
  }
  
  avg <- newdesignmat%*%fixef(fit.rand)+log(0.005)
  std <- as.matrix(sqrt(rowSums((newdesignmat%*%cov_mat)*newdesignmat)))
  ex <- exp(avg%*%c(1,1,1)+1.96*std%*%c(-1,0,1))
  pr <- ex/(1+ex)
  inte <- -log(1-pr)
  for(i in 1:3) inte[as.vector(t(allcont)==0),i] <- NA
  ret <- list(matrix(inte[,1],row_shoe,col_shoe,byrow=1),matrix(inte[,2],row_shoe,col_shoe,byrow=1),matrix(inte[,3],row_shoe,col_shoe,byrow=1))
  names(ret) <- c("Low","Mid","High")
  return(ret)
}

CI_r <- CI_random()

#As shown in Figure 3 the confidence interval is calculated for 3 cut points
sho <- CI$Mid*exp(-0.9915/2) #As done in the presentation of the estimators we multiply so they will be on the same scale
sho[110,] <- 0
sho[190,] <- 0
sho[250,] <- 0
#The first part of figure 3 
image.plot(sho,axes=FALSE)

#intervals for the cut points:
cut1<-110
cut2<-190
cut3<-250


##################################################################################################################
# CI_cut creates an image of the confidence interval of the random effects estimator and the CML estimator 
  #for a specific row of the shoe
# INPUT:
# ======
#cut - the specific row of the shoe on the basis of which the CI would be calculated
#ran - the range of the colums on the basis of which the CI would be calculated
##################################################################################################################
range<-c(1:col_shoe) 
CI_cut<-function(cut=110,ran=range)
{
  CI_plot_cml<-data.frame(range,CI$Low[cut,]*exp(-0.9915/2),CI$Mid[cut,]*exp(-0.9915/2),CI$High[cut,]*exp(-0.9915/2),rep("CML",length(range)))  #As done in the presentation of the estimators we multiply so they will be on the same scale
  CI_plot_rand<-data.frame(range,CI_r$Low[cut,]*exp(-0.9915/2),CI_r$Mid[cut,]*exp(-0.9915/2),CI_r$High[cut,]*exp(-0.9915/2),rep("Random",length(range)))
  names(CI_plot_cml)<-c("x","Low_CI","Mid_CI","High_CI","Estimator")
  names(CI_plot_rand)<-c("x","Low_CI","Mid_CI","High_CI","Estimator")
  CI_plot<-rbind(CI_plot_cml,CI_plot_rand)
  CI_cut<-ggplot(CI_plot) + geom_line(aes(x=x,y=Low_CI,colour=Estimator)) +
  geom_line(aes(x=x,y=High_CI,colour=Estimator)) + geom_ribbon(aes(x=x,ymin=Low_CI,ymax=High_CI,fill=Estimator),alpha=0.5)  +
  geom_line(aes(x=x,y=Mid_CI,colour=Estimator),size=1) +
  theme_bw() + theme(plot.title = element_text(color="black", size=14, face="bold")) +
  scale_fill_brewer(palette="Set1") + scale_color_brewer(palette="Set1") +
  ylab("Probability") + labs(colour="Estimator",fill="Estimator") +
  coord_cartesian(xlim = c(100, 220),ylim=c(0.001,0.0082)) +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  return(CI_cut)
}
CI_cut1<-CI_cut(cut=110,ran=range)
CI_cut2<-CI_cut(cut=190,ran=range)
CI_cut3<-CI_cut(cut=250,ran=range)



# Producing Figures 1-4 in the Web appendix
#Figure 1 in the Web appendix
qplot(as.vector(table(acciden$shoe)), geom="histogram",bins=60,alpha=I(.5),col=I("black"))+xlab("RACs")+ylab("count")
pdf(file ="hist_racs_JASA.pdf", height=6, width=6)
qplot(as.vector(table(acciden$shoe)), geom="histogram",bins=60,alpha=I(.5),col=I("black"))+xlab("RACs")+ylab("count")
dev.off()
min_num_RACs<-min(table(acciden$shoe))
max_num_RACs<-max(table(acciden$shoe))
mean_num_RACs<-mean(table(acciden$shoe))

#Figure 2 in the Web appendix
mat_shoe_acc<- as.matrix(table(n_Acc,shoe)) #matrix of shoes, number of pixels with zero and number of pixels with 1
cont_pix<- as.vector(mat_shoe_acc[1,]+mat_shoe_acc[2,]) # The number of pixels with contact surface is the sum of pixels with zero and with one (pixels with no contact surface are not part of the data)
qplot(cont_pix, geom="histogram",bins=20,alpha=I(.5),col=I("black"))+xlab("Contact surface (number of pixels)")+ylab("count")
pdf(file ="hist_npix_JASA.pdf", height=6, width=6)
qplot(cont_pix, geom="histogram",bins=20,alpha=I(.5),col=I("black"))+xlab("Contact surface (number of pixels)")+ylab("count")
dev.off()
min_num_pix<-min(cont_pix)
max_num_pix<-max(cont_pix)
mean_num_pix<-mean(cont_pix)
#Figure 3 in the Web appendix
tmp <- sumcont
tmp[allcont==0]<-NA
image.plot(t(tmp[nrow(tmp):1,]),axes=FALSE,xlab = 'Comulative Contact Surface') #image of cumulative contact surface
pdf(file ="cum_contact_JASA.pdf", height=6, width=6)
image.plot(t(tmp[nrow(tmp):1,]),axes=FALSE,xlab = 'Comulative Contact Surface') #image of cumulative contact surface
dev.off()
#Figure 4 in the Web appendix
dat<-data.frame(cbind(cont_pix,mat_shoe_acc[2,]))
p1 <- ggplot( dat,aes(x = cont_pix, y = mat_shoe_acc[2,]))
p2<- p1 + geom_point(color="dark grey")+geom_smooth(method = "lm", se = FALSE,color=" black") + labs(x="Contact surface (number of pixels)", y = "number of RACs") 
m <- lm(dat$V2 ~ dat$cont_pix)
a <- signif(coef(m)[1], digits = 4)
b <- signif(coef(m)[2], digits = 2)
textlab <- paste("y = ",b,"x + ",a, sep="")
R_spea<-cor(dat$cont_pix, y = dat$V2,  method =  "spearman")
p3<- p2 + geom_text(aes(x = 15400, y = 210, label = textlab), color=" black", size=6, parse = FALSE)
p3   +annotate("text", x = 16000, y = 190, label = "spearman's r = 0.1161",color="black", size=6,fontface =2)
         
pdf(file ="scat_contact_rac_JASA.pdf", height=6, width=6)
p3<- p2 + geom_text(aes(x = 15400, y = 210, label = textlab), color=" black", size=6, parse = FALSE)
p3   +annotate("text", x = 16000, y = 190, label = "spearman's r = 0.1161",color="black", size=6,fontface =2)
dev.off()