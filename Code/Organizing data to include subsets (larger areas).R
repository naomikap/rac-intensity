#In this code we adjust the following files to include 14 subsets and 36 subsets

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
#There are 307 columns in each shoe and 395 is the number of rows

setwd("C:\\Users\\kapna\\Dropbox\\naomi-micha\\shoe_data\\Codes for JASA")


col_shoe<-307 #307 is the number of columns in each shoe
row_shoe<-395 #395 is the number of rows in each shoe
num_shoe<-387  #387 is the number of shoes but 386 is the number of shoes with RACs - shoe 127 has no RACS
rel_col_shoe<-150  #out of the 307 columns only 150 are relevant (contain non zero pixels in some shoes)
rel_row_shoe<-300  #out of the 395 rows only 300 are relevant (contain non zero pixels in some shoes)
rel_x_cord<-0.25 #using coordinates as in the locations_data.CSV file the relevant x coordinates are between -.25 and 0.25
rel_Y_cord<-0.5 #the relevant Y coordinates are between -0.5 and 0.5

data_RAC<-read.csv("locations_data.CSV",header=TRUE) 

############################organizing the two data sets such that the RACs are divided to 14 subsets. #############################3
#Adding to the location_data, a variable (sub_area) that indicates the subset the RAC belongs to. 
#The 14 subsets were obtained according to expert knowledge as seen in Figure 4
for (i in 1:length(data_RAC$x))
{
  if (data_RAC$x[i]< 0)
  {
    if ((data_RAC$y[i]< -0.35))
    {
      data_RAC$sub_area[i]<- 1
    }
    if ((data_RAC$y[i]>= -0.35) && (data_RAC$y[i]< -0.15))
    {
      data_RAC$sub_area[i]<- 2
    }
    if ((data_RAC$y[i]>= -0.15) && (data_RAC$y[i]< 0))
    {
      data_RAC$sub_area[i]<- 3
    }
    if ((data_RAC$y[i]>= 0) && (data_RAC$y[i]< 0.15))
    {
      data_RAC$sub_area[i]<- 4
    }
    if ((data_RAC$y[i]>= 0.15))
    {
      data_RAC$sub_area[i]<- 5
    }   
  } 
  else 
  {
    if ((data_RAC$y[i]< -0.35))
    {
      data_RAC$sub_area[i]<- 10
    }
    if ((data_RAC$y[i]>= -0.35) && (data_RAC$y[i]< -0.15))
    {
      data_RAC$sub_area[i]<- 9
    }
    if ((data_RAC$y[i]>= -0.15) && (data_RAC$y[i]< 0))
    {
      data_RAC$sub_area[i]<- 8
    }
    if ((data_RAC$y[i]>= 0) && (data_RAC$y[i]< 0.15))
    {
      data_RAC$sub_area[i]<- 7
    }
    if ((data_RAC$y[i]>= 0.15))
    {
      data_RAC$sub_area[i]<- 6
    }   
    
  }
}

for (i in 1:length(data_RAC$x))
{
  if(data_RAC$sub_area[i]==4)
  {
    if(data_RAC$y[i]>12*(data_RAC$x[i]+0.02)^2-0.07)
    {data_RAC$sub_area[i]<-11}
  }
}

for (i in 1:length(data_RAC$x))
{
  if(data_RAC$sub_area[i]==5)
  {
    if(data_RAC$y[i]< -12*(data_RAC$x[i]+0.01)^2+0.4)
    {data_RAC$sub_area[i]<-12}
  }
}

for (i in 1:length(data_RAC$x))
{
  if(data_RAC$sub_area[i]==6)
  {
    if(data_RAC$y[i]< -20*(data_RAC$x[i] +0.01)^2+0.4)
    {data_RAC$sub_area[i]<-13}
  }
}

for (i in 1:length(data_RAC$x))
{
  if(data_RAC$sub_area[i]==7)
  {
    if(data_RAC$y[i] > 16*(data_RAC$x[i]+0.02)^2-0.07)
    {data_RAC$sub_area[i]<-14}
  }
}

############# We create another partition of 36 subsets which will be used later for simulations
###############Dividing the shoe to 36 sub areas##################
for (i in 1:length(data_RAC$x))
{
  if (data_RAC$x[i]< -0.09)
  {
    if (data_RAC$y[i]< -0.35)
    {
      data_RAC$sub_area2[i]<- 1
    }
    if ((data_RAC$y[i]>= -0.35) && (data_RAC$y[i]< -0.25))
    {
      data_RAC$sub_area2[i]<- 5
    }
    if ((data_RAC$y[i]>= -0.25) && (data_RAC$y[i]< -0.15))
    {
      data_RAC$sub_area2[i]<- 9  
    }
    if ((data_RAC$y[i]>= -0.15) && (data_RAC$y[i]< 0))
    {
      data_RAC$sub_area2[i]<- 13
    }
    if ((data_RAC$y[i]>= 0) && (data_RAC$y[i]< 0.15))
    {
      if(data_RAC$y[i]>(12*(data_RAC$x[i]+0.02)^2-0.07))
      {
        data_RAC$sub_area2[i]<- 18 
      }
      else
      {
        data_RAC$sub_area2[i]<- 17
      }
      
    }
    if ((data_RAC$y[i]>= 0.15) && (data_RAC$y[i]< 0.25))
    {
      if(data_RAC$y[i]< -12*(data_RAC$x[i]+0.01)^2+0.4)
      {
        data_RAC$sub_area2[i]<-24
      }
      else
      {
        data_RAC$sub_area2[i]<-23 
      }
    }   
    
    if ((data_RAC$y[i]>= 0.25))
    {
      if(data_RAC$y[i]< -12*(data_RAC$x[i]+0.01)^2+0.4)
      {
        data_RAC$sub_area2[i]<-30
      }
      else
      {
        data_RAC$sub_area2[i]<-29 
      }
    }   
  } 
  
  if ((data_RAC$x[i]>= -0.09) && (data_RAC$x[i]< 0))
  {
    if ((data_RAC$y[i]< -0.35))
    {
      data_RAC$sub_area2[i]<- 2
    }
    if ((data_RAC$y[i]>= -0.35) && (data_RAC$y[i]< -0.25))
    {
      data_RAC$sub_area2[i]<- 6
    }
    if ((data_RAC$y[i]>= -0.25) && (data_RAC$y[i]< -0.15))
    {
      data_RAC$sub_area2[i]<- 10  
    }
    if ((data_RAC$y[i]>= -0.15) && (data_RAC$y[i]< 0))
    {
      data_RAC$sub_area2[i]<- 14
    }
    if ((data_RAC$y[i]>= 0) && (data_RAC$y[i]< 0.15))
    {
      data_RAC$sub_area2[i]<-19
      
    }
    if ((data_RAC$y[i]>= 0.15) && (data_RAC$y[i]< 0.25))
    {
      data_RAC$sub_area2[i]<-25
    }   
    
    if ((data_RAC$y[i]>= 0.25))
    {
      if(data_RAC$y[i]< -12*(data_RAC$x[i]+0.01)^2+0.4)
      {
        data_RAC$sub_area2[i]<-32
      }
      else
      {
        data_RAC$sub_area2[i]<-31 
      }
    }   
  } 
  
  #####
  
  if ((data_RAC$x[i]>= 0) && (data_RAC$x[i]< 0.05))
  {
    if ((data_RAC$y[i]< -0.35))
    {
      data_RAC$sub_area2[i]<- 3
    }
    if ((data_RAC$y[i]>= -0.35) && (data_RAC$y[i]< -0.25))
    {
      data_RAC$sub_area2[i]<- 7
    }
    if ((data_RAC$y[i]>= -0.25) && (data_RAC$y[i]< -0.15))
    {
      data_RAC$sub_area2[i]<- 11  
    }
    if ((data_RAC$y[i]>= -0.15) && (data_RAC$y[i]< 0))
    {
      data_RAC$sub_area2[i]<- 15
    }
    if ((data_RAC$y[i]>= 0) && (data_RAC$y[i]< 0.15))
    {
      data_RAC$sub_area2[i]<-20
      
    }
    if ((data_RAC$y[i]>= 0.15) && (data_RAC$y[i]< 0.25))
    {
      data_RAC$sub_area2[i]<-26
    }   
    
    if ((data_RAC$y[i]>= 0.25))
    {
      if(data_RAC$y[i]< -20*(data_RAC$x[i] +0.01)^2+0.4)
      {
        data_RAC$sub_area2[i]<-34
      }
      else
      {
        data_RAC$sub_area2[i]<-33 
      }
      
    }
  }
  ########
  
  if (data_RAC$x[i]>= 0.05)
  {
    if ((data_RAC$y[i]< -0.35))
    {
      data_RAC$sub_area2[i]<- 4
    }
    if ((data_RAC$y[i]>= -0.35) && (data_RAC$y[i]< -0.25))
    {
      data_RAC$sub_area2[i]<- 8
    }
    if ((data_RAC$y[i]>= -0.25) && (data_RAC$y[i]< -0.15))
    {
      data_RAC$sub_area2[i]<- 12  
    }
    if ((data_RAC$y[i]>= -0.15) && (data_RAC$y[i]< 0))
    {
      data_RAC$sub_area2[i]<- 16
    }
    if ((data_RAC$y[i]>= 0) && (data_RAC$y[i]< 0.15))
    {
      if(data_RAC$y[i] > 16*(data_RAC$x[i]+0.02)^2-0.07)
      {
        data_RAC$sub_area2[i]<- 21  
      }
      else
      {
        data_RAC$sub_area2[i]<- 22  
      }
      
    }
    if ((data_RAC$y[i]>= 0.15) && (data_RAC$y[i]< 0.25))
    {
      if(data_RAC$y[i]< -20*(data_RAC$x[i] +0.01)^2+0.4)
      {
        data_RAC$sub_area2[i]<- 27  
      }
      else
      {
        data_RAC$sub_area2[i]<- 28  
      }
    }   
    
    if ((data_RAC$y[i]>= 0.25))
    {
      if(data_RAC$y[i]< -20*(data_RAC$x[i] +0.01)^2+0.4)
      {
        data_RAC$sub_area2[i]<- 36  
      }
      else
      {
        data_RAC$sub_area2[i]<- 35  
      }   
      
    }
  }
}

write.csv(data_RAC, "locations_data_include_subset.CSV",row.names=FALSE)

################ Working on the contacts_data.txt file and dividing the shoe to 14 subsets obtained according to expert knowledge as seen in Figure 4 ##################
contact_dat <- readChar("contacts_data.txt",nchars=(col_shoe*row_shoe+2)*num_shoe) #importing character
data <- list() 
#creating a list of the matrices of contact surface of the different shoes  
for(i in 1:num_shoe) {
  data[[i]] <- matrix(as.numeric(unlist(strsplit(substr(contact_dat, 1+(col_shoe*row_shoe+2)*(i-1), (col_shoe*row_shoe+2)*i-2), split="")) ),row_shoe,col_shoe,byrow=1)
}
shoe9rev <- data[[9]]
#flipping shoe 9 -Shoe 9 should be mirrored to correspond to all other shoes
data[[9]] <- data[[9]][,ncol(data[[9]]):1]
shoe9 <- data[[9]]

for(i in 1:num_shoe) {
  temp <- matrix(0,row_shoe,col_shoe)
  temp[floor(row_shoe/2-rel_Y_cord*rel_row_shoe):floor(row_shoe/2+rel_Y_cord*rel_row_shoe),floor(col_shoe/2-rel_x_cord*rel_row_shoe):floor(col_shoe/2+rel_x_cord*rel_row_shoe)] <- data[[i]][floor(row_shoe/2-rel_Y_cord*rel_row_shoe):floor(row_shoe/2+rel_Y_cord*rel_row_shoe),floor(col_shoe/2-rel_x_cord*rel_row_shoe):floor(col_shoe/2+rel_x_cord*rel_row_shoe)]
  data[[i]]<-temp
  data[[i]]<-data[[i]][,ncol(data[[i]]):1]
}

subshape <- function(j,int) 
{
  res <- matrix(0,row_shoe,col_shoe)
  sep <- -c(-0.5,-0.35,-0.15,0,0.15,0.5)
  sepPix <- floor(rel_row_shoe*sep+row_shoe/2)
  
  if(j<=5) {
    y <- sepPix[c(j+1,j)]
    x<- c(1,(col_shoe-1)/2)
  }
  else {
    y <- sepPix[c(12-j,11-j)]
    x<- c((col_shoe-1)/2+1,col_shoe)
  }
  
  res[((y[1])+1):(y[2]),(x[1]):(x[2])] <- 1
  
  bound <- function(x,k) {
    ifelse(k==4,return(min(max(floor(-(12*((x-col_shoe/2)/rel_row_shoe+0.02)^2-0.07)*rel_row_shoe+row_shoe/2),1),row_shoe)),
           ifelse(k==5,return(min(max(floor(-(-12*((x-col_shoe/2)/rel_row_shoe+0.01)^2+0.4)*rel_row_shoe+row_shoe/2),1),row_shoe)),
                  ifelse(k==6,return(min(max(floor(-(-20*((x-col_shoe/2)/rel_row_shoe+0.01)^2+0.4)*rel_row_shoe+row_shoe/2),1),row_shoe)),
                         ifelse(k==7,return(min(max(floor(-(16*((x-col_shoe/2)/rel_row_shoe+0.02)^2-0.07)*rel_row_shoe+row_shoe/2),1),row_shoe)),
                                return(1)
                         ))))
  }
  
  if(j==4) {
    if(int==0){
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,4)-1),i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[bound(i,4):row_shoe,i] <- 0
      }
    }
  }
  if(j==5)  {
    if(int ==0){
      for(i in (x[1]):(x[2])) {
        res[bound(i,5):row_shoe,i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,5)-1),i] <- 0
      }
    }
  }
  
  if(j==6){
    if(int ==0){
      for(i in (x[1]):(x[2])) {
        res[bound(i,6):row_shoe,i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,6)-1),i] <- 0
      }
    }
  }
  
  if(j==7) {
    if(int==0){
      for(i in (x[1]):(x[2])) {
        res[1:(bound(i,7)-1),i] <- 0
      }
    }
    else {
      for(i in (x[1]):(x[2])) {
        res[bound(i,7):row_shoe,i] <- 0
      }
    }
  }
  
  return(res)
  
}

subshapes <- list()
for(i in 1:10) {
  subshapes[[i]] <- subshape(i,0)
}

for(i in 4:7) {
  subshapes[[i+7]] <- subshape(i,1)
}

subAreaCont <- function(i,j) 
{
  #building the matrix of the sum of contact surface pixels in each region. 
  data[[i]]*subshapes[[j]]
}


chi<-matrix(0,num_shoe,14) #a martix of the contact surface in each of the 14 subsets
#For each shoe and each subset, the number of pixels that have a contact surface 
for (i in 1:num_shoe)
{
  for(j in 1:14)
  {
    chi[i,j]<-sum(subAreaCont(i,j))
  }
}
chi<-chi[-127,] # there is no  RACS in shoe 127 and therefore it is removed
chi<-chi/10000 # chi is divided in 10000 and so now it is be between zero and 1

write.csv(chi, "chi.CSV",row.names=FALSE)


#we now create chi2, a matrix of contact surfaces for a partition of 36 subsets that will also be used and compared to a partition of 14 subsets
part_36 <- list()
xpart <- c(0,floor(-0.09*rel_row_shoe+col_shoe/2),floor(0.05*rel_row_shoe+col_shoe/2),col_shoe+1)
ypart <- c(0,floor(-0.25*rel_row_shoe+row_shoe/2),floor(0.25*rel_row_shoe+row_shoe/2),row_shoe+1)
for(i in 1:3){
  for(j in 1:3) {
    part_36[[3*(3-i)+j]] <- matrix(0,row_shoe,col_shoe)
    part_36[[3*(3-i)+j]][ypart[i]:(ypart[i+1]-1),xpart[j]:(xpart[j+1]-1)] <- 1
  }  
}
subshapes_new <- list()
subshapes_new[[1]] <- part_36[[1]]*subshapes[[1]]
subshapes_new[[2]] <- part_36[[2]]*subshapes[[1]]
subshapes_new[[3]] <- part_36[[2]]*subshapes[[10]]
subshapes_new[[4]] <- part_36[[3]]*subshapes[[10]]
subshapes_new[[5]] <- part_36[[1]]*subshapes[[2]]
subshapes_new[[6]] <- part_36[[2]]*subshapes[[2]]
subshapes_new[[7]] <- part_36[[2]]*subshapes[[9]]
subshapes_new[[8]] <- part_36[[3]]*subshapes[[9]]

subshapes_new[[9]] <- part_36[[4]]*subshapes[[2]]
subshapes_new[[10]] <- part_36[[5]]*subshapes[[2]]
subshapes_new[[11]] <- part_36[[5]]*subshapes[[9]]
subshapes_new[[12]] <- part_36[[6]]*subshapes[[9]]
subshapes_new[[13]] <- part_36[[4]]*subshapes[[3]]
subshapes_new[[14]] <- part_36[[5]]*subshapes[[3]]
subshapes_new[[15]] <- part_36[[5]]*subshapes[[8]]
subshapes_new[[16]] <- part_36[[6]]*subshapes[[8]]

subshapes_new[[17]] <- part_36[[4]]*subshapes[[4]]
subshapes_new[[18]] <- part_36[[4]]*subshapes[[11]]
subshapes_new[[19]] <- part_36[[5]]*subshapes[[11]]
subshapes_new[[20]] <- part_36[[5]]*subshapes[[14]]
subshapes_new[[21]] <- part_36[[6]]*subshapes[[14]]
subshapes_new[[22]] <- part_36[[6]]*subshapes[[7]]

subshapes_new[[23]] <- part_36[[4]]*subshapes[[5]]
subshapes_new[[24]] <- part_36[[4]]*subshapes[[12]]
subshapes_new[[25]] <- part_36[[5]]*subshapes[[12]]
subshapes_new[[26]] <- part_36[[5]]*subshapes[[13]]
subshapes_new[[27]] <- part_36[[6]]*subshapes[[13]]
subshapes_new[[28]] <- part_36[[6]]*subshapes[[6]]

subshapes_new[[29]] <- part_36[[7]]*subshapes[[5]]
subshapes_new[[30]] <- part_36[[7]]*subshapes[[12]]
subshapes_new[[31]] <- part_36[[8]]*subshapes[[5]]
subshapes_new[[32]] <- part_36[[8]]*subshapes[[12]]
subshapes_new[[33]] <- part_36[[8]]*subshapes[[6]]
subshapes_new[[34]] <- part_36[[8]]*subshapes[[13]]
subshapes_new[[35]] <- part_36[[9]]*subshapes[[6]]
subshapes_new[[36]] <- part_36[[9]]*subshapes[[13]]

subAreaCont_36 <- function(i,j) data[[i]]*subshapes_new[[j]]

chi2<-matrix(0,num_shoe,36)
for (i in 1:num_shoe)
{
  for(j in 1:36)
  {
    chi2[i,j]<-sum(subAreaCont_36(i,j))
  }
}
chi2<-chi2/10000 # chi2 is divided in 10000 and so now it is be between zero and 1
chi2<-chi2[-127,] # there is no shoe 127 - this is a shoe with no RACs

write.csv(chi2, "chi2.CSV",row.names=FALSE)


##########The following code is used for adjustment of plot 5.  We determine the contour of a shoe and create a nice plot
n_sub<-ncol(chi)
allcont <- data[[1]]
for(i in 2:num_shoe) allcont <- allcont+data[[i]] #this is the contact of all shoes 
allcont <- (allcont>=18)*1 #We are interseted in the contour of the shoe.  As an approximation to the cumulative contour we are looking at pixels that appear in more than 18 shoes 
cont_use<-list()
for(i in 1:n_sub)
{
  cont_use[[i]]<-subshapes[[i]]*allcont
  write.csv(cont_use[[i]], paste("cont_use",i,".CSV" ,sep = "", collapse = NULL),row.names=FALSE)
  
}
write.csv(allcont, "allcont.CSV",row.names=FALSE)

