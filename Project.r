##### Project Multivariate Analysis #####
############# Group 5 ###################

# Import data #####
# Change Path
data <- read.csv("",header=TRUE)
data$Region <- as.factor(data$Region)
data$Country <- as.factor(data$Country)
data$Class<- as.factor(data$Class)
data$GDP.per.capita <- as.numeric(data$GDP.per.capita)


##### Pre-Process ###########################################################
### Imputation : PMM ####################################################

library("mice")

#Applying Predictive Mean Matching 5 times and choosing one of them
imp <- mice(data[,-c(1:4)], m = 5, method = "pmm",seed=20) # Impute missing values
data_imp <- complete(imp,2)
data_imp <- cbind(data[,c(1:4)],data_imp)

#Displaying the different PMM imputations
stripplot(imp, pch = 20, cex = 1.2)[1:9]
stripplot(imp, pch = 20, cex = 1.2)[10:17]

### Standardizing Imputed Data ##########################
data <- scale(data_imp[,-c(1:4)])
data <- cbind(data_imp[,c(1:4)],data)

### Outlier Analysis ##########################

out <- boxplot.stats(data$Mean.years.of.schooling)$out
out_school <- which(data$Mean.years.of.schooling %in% c(out))

out <- boxplot.stats(data$Carbon.Dioxide.Emissions.per.GDP)$out
out_carbon <- which(data$Carbon.Dioxide.Emissions.per.GDP %in% c(out))

out <- boxplot.stats(data$Change.Florest.Area)$out
out_florest <- which(data$Change.Florest.Area %in% c(out))

out <- boxplot.stats(data$Gender.Inequality.Index)$out
out_gender <- which(data$Gender.Inequality.Index %in% c(out))

out <- boxplot.stats(data$Child.malnutrition)$out
out_child <- which(data$Child.malnutrition %in% c(out))

out <- boxplot.stats(data$Current.health.expenditure)$out
out_health <- which(data$Current.health.expenditure %in% c(out))

out <- boxplot.stats(data$Homicide.Rate)$out
out_homicide <- which(data$Homicide.Rate %in% c(out))

out <- boxplot.stats(data$Net.migration.rate)$out
out_migration <- which(data$Net.migration.rate %in% c(out))

out <- boxplot.stats(data$Unemployment)$out
out_unemployment <- which(data$Unemployment %in% c(out))

out <- boxplot.stats(data$Vulnerable.employment)$out
out_employment <- which(data$Vulnerable.employment %in% c(out))

out <- boxplot.stats(data$Suicide.Rate)$out
out_suicide <- which(data$Suicide.Rate %in% c(out))

out <- boxplot.stats(data$Contraceptive.prevalence)$out
out_contraceptive <- which(data$Contraceptive.prevalence %in% c(out))

out <- boxplot.stats(data$GDP.per.Capita)$out
out_gdp <- which(data$GDP.per.Capita %in% c(out))

out <- boxplot.stats(data$Life.expectancy.at.birth)$out
out_life <- which(data$Life.expectancy.at.birth %in% c(out))

out <- boxplot.stats(data$Median.age)$out
out_median <- which(data$Median.age %in% c(out))

out <- boxplot.stats(data$Refugees.by.country.of.origin)$out
out_ref <- which(data$Refugees.by.country.of.origin %in% c(out))

out <- boxplot.stats(data$Share.parliament.by.women)$out
out_share <- which(data$Share.parliament.by.women %in% c(out))


table(c(out_gdp,out_share,out_ref, out_median, out_life, out_school, out_carbon, out_florest, out_gender, out_child, out_health, out_homicide, out_migration, out_unemployment, out_employment, out_suicide, out_contraceptive))

library("ggplot2")
library("tidyr")
ggplot(gather(data[,-c(1:4)],"x","y"), aes(x, y)) + geom_boxplot() + coord_flip()


### Data Reduction ###################################
# PCA : data_pca #############################################################

#Applying PCA to the normalized data
pca <- prcomp(data[,-c(1:4)], center=TRUE, scale=TRUE)
summary(pca) #Helps determine how many PC's to retain

#Transformation of the data into PCA coordenates
data_pca<-data.frame(predict(pca,data))
summary(data_pca)

#Retaining 7 PC's --> more than 80% variability
data_pca <- data_pca[,1:7] 
data_pca<-cbind(data[,c(1:4)],data_pca)

#Magnitude of data reduction
ncol(data) - ncol(data_pca)

# Corr : data_corr ############################################################

library(corrplot)

#Correlation matrix of standardized data
corrplot(cor(data[,-c(1:4)]),tl.col = "black",tl.srt = 45,col=colorRampPalette(c("#F55118","#FF9933","#FFFFFF","#66CCCC","#298DF7"))(50) )
#Correlation matrix after removing high correlated variables
corrplot(cor(data[,-c(1:4,8,9,14,19)]),tl.col = "black",tl.srt = 45,col=colorRampPalette(c("#F55118","#FF9933","#FFFFFF","#66CCCC","#298DF7"))(50))
data_corr <- data[,-c(8,9,14,19)]

#Magnitude of data reduction
ncol(data) - ncol(data_corr)

# Corr + PCA : data_corr_pca #######################################################

#Applying PCA to the data_corr
corr_pca <- prcomp(data_corr[,-c(1:4)], center=TRUE, scale=TRUE)
summary(corr_pca)
corr_pca

#Transformation of the data_corr into PCA coordenates
data_corr_pca<-data.frame(predict(corr_pca,data_corr))
summary(data_corr_pca)

#Retaining 7 PC's --> more than 80% variability
data_corr_pca <- data_corr_pca[,1:7]
data_corr_pca<-cbind(data[,c(1:4)],data_corr_pca)

#Magnitude of data reduction
ncol(data) - ncol(data_corr_pca)


##### Unsupervised Learning : Cluster Analysis ########################################

library("cluster")
library("clusterCrit")

### Internal Criteria : Number of clusters #########

#Auxiliary fuction for internal criteria, helping determine the number of clusters to use
## dataset : dataset to consider while clustering
## max_cluster_number : maximum best number of clusters possible
## criteria : one of the criteria available in intCriteria(package clusterCrit)
#1 This function returns an order list of values of the criteria for k clusters with k between 2 and max_number_cluster
#2 This function returns the maximum value of the previous list
#3 This function returns the k for which we get the maximum criteria
# The function returns #1,#2,#3 for the following clustering method : single, average and complete linkage, Ward Method, K-Means and K-Medoids  
intCrit <- function(dataset,max_cluster_number,criteria) {
  setClass(Class="Criteria",representation(all_crit_kmeans="vector",all_crit_kmedoids="vector",all_crit_single="vector",all_crit_complete="vector",all_crit_average="vector",all_crit_ward="vector",
                                           best_kmeans="numeric",best_kmedoids="numeric",best_single="numeric",best_complete="numeric",best_average="numeric",best_ward="numeric",
                                           nc_kmeans="numeric",nc_kmedoids="numeric",nc_single="numeric",nc_complete="numeric",nc_average="numeric",nc_ward="numeric"))
  
  all_crit_kmeans <- c()
  all_crit_kmedoids <- c()
  all_crit_single <- c()
  all_crit_complete <- c()
  all_crit_average <- c()
  all_crit_ward <- c()
  
  for(b in 2:max_cluster_number){
    
    data_single<-agnes(dataset, metric = "euclidean",stand = FALSE, method = "single", keep.data = FALSE)
    data_complete<-agnes(dataset, metric = "euclidean",stand = FALSE, method = "complete", keep.data = FALSE)
    data_average<-agnes(dataset, metric = "euclidean",stand = FALSE, method = "average", keep.data = FALSE)
    data_ward<-agnes(dataset, metric = "euclidean",stand = FALSE, method = "ward", keep.data = FALSE)
    
    data_single_cluster<-cutree(data_single,b)
    data_complete_cluster<-cutree(data_complete,b)
    data_average_cluster<-cutree(data_average,b)
    data_ward_cluster<-cutree(data_ward,b)
    
    crit_single <- as.numeric(intCriteria(as.matrix(dataset),data_single_cluster,c(criteria)))
    crit_complete <- as.numeric(intCriteria(as.matrix(dataset),data_complete_cluster,c(criteria)))
    crit_average <- as.numeric(intCriteria(as.matrix(dataset),data_average_cluster,c(criteria)))
    crit_ward <- as.numeric(intCriteria(as.matrix(dataset),data_ward_cluster,c(criteria)))
    
    a1 <-0
    a2 <-0
    for(i in 1:10){
      #
      # K-means
      #
      k1 <- kmeans(dataset,centers = b,nstart = 20)
      crit1 <- as.numeric(intCriteria(as.matrix(dataset),k1$cluster,c(criteria)))
      if(crit1 > a1){
        a1 <- crit1 }
      #
      # K-medoids
      #
      k2 <- pam(dataset,k=b,metric = "euclidean",stand = FALSE)
      crit2 <- as.numeric(intCriteria(as.matrix(dataset),k2$cluster,c(criteria)))
      if(crit2 > a2){
        a2 <- crit2 }
    }
    all_crit_kmeans <- c(all_crit_kmeans,a1)
    all_crit_kmedoids <- c(all_crit_kmedoids,a2)
    all_crit_single <- c(all_crit_single,crit_single)
    all_crit_complete<- c(all_crit_complete,crit_complete)
    all_crit_average <- c(all_crit_average,crit_average)
    all_crit_ward <- c(all_crit_ward,crit_ward)
  }
  return(new("Criteria",
             all_crit_kmeans=all_crit_kmeans,
             all_crit_kmedoids=all_crit_kmedoids,
             all_crit_single=all_crit_single,
             all_crit_complete=all_crit_complete,
             all_crit_average=all_crit_average,
             all_crit_ward=all_crit_ward,
             best_kmeans=max(all_crit_kmeans),
             best_kmedoids=max(all_crit_kmedoids),
             best_single=max(all_crit_single),
             best_complete=max(all_crit_complete),
             best_average=max(all_crit_average),
             best_ward=max(all_crit_ward),
             nc_kmeans=which.max(all_crit_kmeans)+1,
             nc_kmedoids=which.max(all_crit_kmedoids)+1,
             nc_single=which.max(all_crit_single)+1,
             nc_complete=which.max(all_crit_complete)+1,
             nc_average=which.max(all_crit_average)+1,
             nc_ward=which.max(all_crit_ward)+1))
}

crit1 <- intCrit(data_corr_pca[,-c(1:4)],30,"Calinski_Harabasz");crit1
crit2 <- intCrit(data_corr_pca[,-c(1:4)],30,"Dunn");crit2
crit3 <- intCrit(data_corr_pca[,-c(1:4)],30,"Gamma");crit3
crit4 <- intCrit(data_corr_pca[,-c(1:4)],30,"Ratkowsky_Lance");crit4
crit5 <- intCrit(data_corr_pca[,-c(1:4)],30,"Wemmert_Gancarski");crit5

### Clustering ###################################################

# Single, complete and average linkage + Ward Method
# data_corr_pca can be substituted by data, data_pca or data_corr
data_single<-agnes(data_corr_pca[,-c(1:4)], metric = "euclidean",stand = FALSE, method = "single", keep.data = FALSE)
data_complete<-agnes(data_corr_pca[,-c(1:4)], metric = "euclidean",stand = FALSE, method = "complete", keep.data = FALSE)
data_average<-agnes(data_corr_pca[,-c(1:4)], metric = "euclidean",stand = FALSE, method = "average", keep.data = FALSE)
data_ward<-agnes(data_corr_pca[,-c(1:4)], metric = "euclidean",stand = FALSE, method = "ward", keep.data = FALSE)

#Dendogram for Single, complete and average linkage + Ward Method
par(mfrow=c(2,2))
pltree(data_single,main="Single linkage", cex=0.83,xlab="")
pltree(data_complete,main="Complete linkage",cex=0.83,xlab="")
pltree(data_average,main="Average linkage", cex=0.83,xlab="")
pltree(data_ward,main="Ward Method", cex=0.83,xlab="")
par(mfrow=c(1,1))

#Partioning dataset into two clusters (now considering k-means and k-medoids as well)
data_single_cluster<-cutree(data_single,2)
data_complete_cluster<-cutree(data_complete,2)
data_average_cluster<-cutree(data_average,2)
data_ward_cluster<-cutree(data_ward,2)
k1<-kmeans(data_corr_pca[,-c(1:4)],centers=2)
k2<-pam(data_corr_pca[,-c(1:4)],k=2,metric = "euclidean",stand = FALSE)

#Confusion matrix of different partitioning methods, considering the real classes of the original data
table(data_single_cluster,data[,4])
table(data_complete_cluster,data[,4])
table(data_average_cluster,data[,4])
table(data_ward_cluster,data[,4])
table(k1$cluster,data[,4])
table(k2$cluster,data[,4])

### External Criteria ##############################################

#Auxiliary fuction for external criteria
## clusters : classes of the data regarding clustering method
## classification : classes of the data regarding original classes of the data
# This fuction returns a list with several external criteria : (overall accuracy, negative predicted value, positive predicted value, specificity, sensitivity, overall F1-score)
extCrit <- function(clusters,classification){
  switch = FALSE
  conf <- table(clusters,classification)
  if(sum(diag(conf)) < sum(diag(conf[,c(2,1)])) ) {
    conf <- conf[,c(2,1)]
    switch = TRUE
  }
  accuracy <- round( sum(diag(conf)) / sum(conf) , 3)
  precision <- round( diag(conf) / colSums(conf) , 3)
  if(switch) {
    precision <- precision[c(2,1)]
  }
  recall <-  round( diag(conf) / rowSums(conf) , 3)
  F1 <- round( 2 / (precision^(-1)+recall^(-1)) ,3)
  return(c(accuracy,precision,recall,mean(F1)))
}

extCrit(data_single_cluster,data[,4])
extCrit(data_complete_cluster,data[,4])
extCrit(data_average_cluster,data[,4])
extCrit(data_ward_cluster,data[,4])
extCrit(k1$cluster,data[,4])
extCrit(k2$cluster,data[,4])

### Relative : Noise ####
library(mvtnorm)
library(corrplot)
library(psych)
library(cluster)
library(aricode)

magree_ward<-c()
magree_kmean<-c()
magree_kmedoid<-c()

# Step of the variance of the noise
v<-seq(0, 0.2, 0.001)

# Adding normal distributed noise to data_corr_pca
# Calculate agreement between the clusters obtained by data_corr_pca and data_corr_pca + noise, using Ward's Method, K-Means and K-Medoids
for(t in v){
  agree_ward<-c()
  agree_kmean<-c()
  agree_kmedoid<-c()
  for(i in 1:20){
    error<-rmvnorm(154,mean=c(1:7)*0,sigma=diag(c(1:7)*t))
    data2<-data_corr_pca[,-c(1:4)]+error
    data2<-cbind(data_corr_pca[,c(1:4)], data2)
    data_ward<-agnes(data_corr_pca[,-c(1:4)],metric="euclidean",stand = FALSE,method="ward",keep.data = FALSE)
    data2_ward<-agnes(data2[,-c(1:4)],metric="euclidean",stand = FALSE,method="ward",keep.data = FALSE)
    kmean<-kmeans(data_corr_pca[,-c(1:4)],centers=2)
    kmean_2<-kmeans(data2[,-c(1:4)],centers=2)
    kmedoid<-pam(data_corr_pca[,-c(1:4)],k=2,metric = "euclidean",stand = FALSE)
    kmedoid_2<-pam(data2[,-c(1:4)],k=2,metric = "euclidean",stand = FALSE)
    
    data_ward_cluster<-cutree(data_ward,2)
    data2_ward_cluster<-cutree(data2_ward,2)
    
    aux_ward<-sum(diag(table(data_ward_cluster,data2_ward_cluster)))/nrow(data_corr_pca)
    if (aux_ward < 0.5) {
      aux_ward<- 1-aux_ward
    }
    agree_ward<-append(agree_ward,aux_ward)
    
    aux_kmean<-sum(diag(table(kmean$cluster, kmean_2$cluster)))/nrow(data_corr_pca )
    if (aux_kmean< 0.5) {
      aux_kmean<- 1-aux_kmean
    }
    agree_kmean<-append(agree_kmean,aux_kmean)
    
    aux_kmedoid<-sum(diag(table(kmedoid$cluster, kmedoid_2$cluster)))/nrow(data_corr_pca )
    if (aux_kmedoid< 0.5) {
      aux_kmedoid<- 1-aux_kmedoid
    }
    agree_kmedoid<-append(agree_kmedoid,aux_kmedoid)
    
    
    
  }
  
  magree_ward<- append(magree_ward, (mean(agree_ward)))
  magree_kmean<- append(magree_kmean, (mean(agree_kmean)))
  magree_kmedoid<- append(magree_kmedoid, (mean(agree_kmedoid)))
  
}

#Plot of the noise considering the different methods
plot(v, magree_ward,col="#FF9933", type="l", ylim=c(0.6,1), xlab="Variance of Noise", ylab="Agreement")
lines(v, magree_kmean, col="#66CCCC", type="l")
lines(v, magree_kmedoid, col="#298DF7", type="l")
legend( 0, 0.7, c("Ward's Method","K-Means","K-Medoids"), lwd=c(2,2), col=c("#FF9933","#66CCCC","#298DF7"))

##### Supervised Learning #####
### DISCRIMINATOR ANALYSIS with original classes:####
# Linear Discriminator####
library(caret)
library(MASS)
library(comprehenr)

#function that trains the linear discriminator, with leave one out cross validation
fit  <- lda(x = data_corr_pca[,-c(1,2,3,4)],grouping = data_corr_pca[,c(4)],CV=TRUE)

#calculating the criteria, accuracy, recalls, precisions and F1-score
tt<-table(fit$class,data_corr_pca[,c(4)]);tt
acc<-sum(diag(tt))/sum(tt);acc
recalls <- to_vec(for(i in 1:2) if(TRUE) diag(tt)[i]/sum(tt[,i]));recalls
precisions  <- to_vec(for(i in 1:2) if(TRUE) diag(tt)[i]/sum(tt[i,]));precisions
F1_measure  <-to_vec(for(i in 1:2) if(TRUE) 2*recalls[i]*precisions[i]/(recalls[i]+precisions[i])  );F1_measure
F1 <- mean(c(F1_measure));F1


# Quadratic Discriminator####

#function that trains the quadratic discriminator, with leave one out cross validation
fit2  <- qda(x = data_corr_pca[,-c(1,2,3,4)],grouping = data_corr_pca[,c(4)],CV=TRUE)

#calculating the criteria, accuracy, recalls, precisions and F1-score
tt2<-table(fit2$class,data_corr_pca[,c(4)]);tt2
acc2<-sum(diag(tt2))/sum(tt2);acc2
recalls2 <- to_vec(for(i in 1:2) if(TRUE) diag(tt2)[i]/sum(tt2[,i]));recalls2
precisions2 <- to_vec(for(i in 1:2) if(TRUE) diag(tt2)[i]/sum(tt2[i,]));precisions2
F1_measure2 <-to_vec(for(i in 1:2) if(TRUE) 2*recalls2[i]*precisions2[i]/(recalls2[i]+precisions2[i])  );F1_measure2
F1_2<- mean(c(F1_measure2));F1_2

# Decision Trees#######
library(rpart)
library(rpart.plot)
library(caTools)

#we will add the measure of every simulation and then average in the end
#this code initializes the variables that will keep the added measures.
tt_tree_aux <- table(c(0,1), c(0,1)) - table(c(0,1), c(0,1))
acc_tree_aux <- 0
recalls_tree_aux <- c(0,0)
precisions_tree_aux  <- c(0,0)
F1_measure_tree_aux  <- c(0,0)
F1_tree_aux <- 0

#this cycle will train a tree 20 times and then average the criteria
for (i in 1:20) {
  
  #split the data into training set and testing set
  spl = sample.split(data_corr_pca, SplitRatio = 0.9)
  train = subset(data_corr_pca, spl==TRUE)
  test = subset(data_corr_pca, spl==FALSE)
  
  #training the tree using the training set
  fit_tree <- rpart(train[,c(4)]~., data = train[,-c(1,2,3,4)], method = 'class');fit_tree

  #predictions on the testing set
  predictions_tree <- predict(fit_tree, newdata = test[,-c(1,2,3,4)], type = 'class')
  
  #calculating the criteria, accuracy, recalls, precisions and F1-score
  tt_tree <- table(predictions_tree, test[,c(4)])
  acc_tree <- sum(diag(tt_tree))/sum(tt_tree)
  recalls_tree <- to_vec(for(i in 1:2) if(TRUE) diag(tt_tree)[i]/sum(tt_tree[,i]))
  precisions_tree  <- to_vec(for(i in 1:2) if(TRUE) diag(tt_tree)[i]/sum(tt_tree[i,]))
  F1_measure_tree  <- to_vec(for(i in 1:2) if(TRUE) 2*recalls_tree[i]*precisions_tree[i]/(recalls_tree[i]+precisions_tree[i])  )
  F1_tree <- mean(c(F1_measure_tree))
  
  
  #This is to add the measure with previous simulations, in the end the average will be made
  tt_tree_aux <- tt_tree_aux + tt_tree
  acc_tree_aux <- acc_tree_aux + acc_tree
  recalls_tree_aux <- recalls_tree_aux + recalls_tree
  precisions_tree_aux  <- precisions_tree_aux + precisions_tree
  F1_measure_tree_aux  <- F1_measure_tree_aux + F1_measure_tree
  F1_tree_aux <- F1_tree_aux + F1_tree
}
#this plots the last calculated tree
rpart.plot(fit_tree,box.palette = colorRampPalette(c("#F55118","#FF9933","#FFFFFF","#66CCCC","#298DF7"))(50))

#finaly, the criteria are averaged out
tt_tree_final <- tt_tree_aux/20;tt_tree_final
acc_tree_final <- acc_tree_aux/20;acc_tree_final
recalls_tree_final <- recalls_tree_aux/20;recalls_tree_final
precisions_tree_final  <- precisions_tree_aux/20;precisions_tree_final
F1_measure_tree_final  <- F1_measure_tree_aux/20;F1_measure_tree_final
F1_tree_final <- F1_tree_aux/20;F1_tree_final

# Min-max normalize the data for the neural network analyses####

#this auxiliar function makes the maxinum equal to 1, the minimum equal to 0 and changes
#every other value proportionaly
minmax <- function(x){(x-min(x))/(max(x)-min(x))}

#here we make the min-max normalization of the data so the neural network provides good results
data_corr_pca2 <- data_corr_pca
data_corr_pca2$PC1 <- minmax(data_corr_pca2$PC1)
data_corr_pca2$PC2 <- minmax(data_corr_pca2$PC2)
data_corr_pca2$PC3 <- minmax(data_corr_pca2$PC3)
data_corr_pca2$PC4 <- minmax(data_corr_pca2$PC4)
data_corr_pca2$PC5 <- minmax(data_corr_pca2$PC5)
data_corr_pca2$PC6 <- minmax(data_corr_pca2$PC6)
data_corr_pca2$PC7 <- minmax(data_corr_pca2$PC7)


# Neural Network##########################################################################################################
library(neuralnet)
#the same idea as we did on the decision trees, add the measures of every simulation and average in the end
tt_nn_aux <- table(c(0,1), c(0,1)) - table(c(0,1), c(0,1))
acc_nn_aux <- 0
recalls_nn_aux <- c(0,0)
precisions_nn_aux  <- c(0,0)
F1_measure_nn_aux  <- c(0,0)
F1_nn_aux <- 0


for (i in 1:20) {
  
  #split the data into raining and testing
  spl = sample.split(data_corr_pca2, SplitRatio = 0.9)
  train = subset(data_corr_pca2, spl==TRUE)
  test = subset(data_corr_pca2, spl==FALSE)
  
  #train the neural network usint the training set, the hidden variable can be changed to change the number
  #of layers ofthe network and the number of neurons in each layer.
  fit_nn <- neuralnet(train[,c(4)]~., data=train[,-c(1,2,3,4)], hidden=c(5,2,4), linear.output=FALSE, threshold=0.01)
  
  #make the predictions on the testing set
  predictions_nn <- predict(fit_nn,  test[,-c(1,2,3,4)])#;predictions_nn
  
  #the neural network activates every output neuron with a value from 0 to 1
  #the neuron with th higher value is the chosen classification
  predictions_nn2 <- to_vec(for(i in 1:nrow(predictions_nn)) if(TRUE) which.max(predictions_nn[i,]) )#;predictions_nn2
  
  #calculating the criteria, accuracy, recalls, precisions and F1-score
  tt_nn <- table(predictions_nn2, test[,c(4)])
  acc_nn <- sum(diag(tt_nn))/sum(tt_nn)
  recalls_nn <- to_vec(for(i in 1:2) if(TRUE) diag(tt_nn)[i]/sum(tt_nn[,i]))
  precisions_nn  <- to_vec(for(i in 1:2) if(TRUE) diag(tt_nn)[i]/sum(tt_nn[i,]))
  F1_measure_nn  <- to_vec(for(i in 1:2) if(TRUE) 2*recalls_nn[i]*precisions_nn[i]/(recalls_nn[i]+precisions_nn[i])  )
  F1_nn <- mean(c(F1_measure_nn))
  
  
  #This is to add the measure with previous simulations, in the end the average will be made
  tt_nn_aux <- tt_nn_aux + tt_nn
  acc_nn_aux <- acc_nn_aux + acc_nn
  recalls_nn_aux <- recalls_nn_aux + recalls_nn
  precisions_nn_aux  <- precisions_nn_aux + precisions_nn
  F1_measure_nn_aux  <- F1_measure_nn_aux + F1_measure_nn
  F1_nn_aux <- F1_nn_aux + F1_nn
  
}
#plot the final network calculated
plot(fit_nn)

#finaly, the criteria are averaged out
tt_nn_final <- tt_nn_aux/20;tt_nn_final
acc_nn_final <- acc_nn_aux/20;acc_nn_final
recalls_nn_final <- recalls_nn_aux/20;recalls_nn_final
precisions_nn_final  <- precisions_nn_aux/20;precisions_nn_final
F1_measure_nn_final  <- F1_measure_nn_aux/20;F1_measure_nn_final
F1_nn_final <- F1_nn_aux/20;F1_nn_final

###DISCRIMINATOR ANALYSIS with partitions made by the K-means method:####
# Linear Discriminator####
library(caret)
library(MASS)
library(comprehenr)

fit  <- lda(x = data_corr_pca[,-c(1,2,3,4)],grouping = k1$cluster,CV=TRUE)

tt<-table(fit$class,k1$cluster);tt
acc<-sum(diag(tt))/sum(tt);acc

recalls <- to_vec(for(i in 1:2) if(TRUE) diag(tt)[i]/sum(tt[,i]));recalls
precisions  <- to_vec(for(i in 1:2) if(TRUE) diag(tt)[i]/sum(tt[i,]));precisions
F1_measure  <-to_vec(for(i in 1:2) if(TRUE) 2*recalls[i]*precisions[i]/(recalls[i]+precisions[i])  );F1_measure
F1 <- mean(c(F1_measure));F1


# Quadratic Discriminator####


fit2  <- qda(x = data_corr_pca[,-c(1,2,3,4)],grouping = k1$cluster,CV=TRUE)

tt2<-table(fit2$class,k1$cluster);tt2

acc2<-sum(diag(tt2))/sum(tt2);acc2

recalls2 <- to_vec(for(i in 1:2) if(TRUE) diag(tt2)[i]/sum(tt2[,i]));recalls2

precisions2 <- to_vec(for(i in 1:2) if(TRUE) diag(tt2)[i]/sum(tt2[i,]));precisions2

F1_measure2 <-to_vec(for(i in 1:2) if(TRUE) 2*recalls2[i]*precisions2[i]/(recalls2[i]+precisions2[i])  );F1_measure2

F1_2<- mean(c(F1_measure2));F1_2



# Decision Trees#######
library(rpart)
library(rpart.plot)
library(caTools)


tt_tree_aux <- table(c(0,1), c(0,1)) - table(c(0,1), c(0,1))
acc_tree_aux <- 0
recalls_tree_aux <- c(0,0)
precisions_tree_aux  <- c(0,0)
F1_measure_tree_aux  <- c(0,0)
F1_tree_aux <- 0


for (i in 1:20) {
  
  
  spl = sample.split(data_corr_pca, SplitRatio = 0.9)
  train = subset(data_corr_pca, spl==TRUE)
  train_clusters = subset(k1$cluster, spl==TRUE)
  test = subset(data_corr_pca, spl==FALSE)
  test_clusters = subset(k1$cluster, spl==FALSE)

  fit_tree <- rpart(train_clusters~., data = train[,-c(1,2,3,4)], method = 'class')
  
  predictions_tree <- predict(fit_tree, newdata = test[,-c(1,2,3,4)], type = 'class')
  
  
  tt_tree <- table(predictions_tree, test_clusters)
  acc_tree <- sum(diag(tt_tree))/sum(tt_tree)
  recalls_tree <- to_vec(for(i in 1:2) if(TRUE) diag(tt_tree)[i]/sum(tt_tree[,i]))
  precisions_tree  <- to_vec(for(i in 1:2) if(TRUE) diag(tt_tree)[i]/sum(tt_tree[i,]))
  F1_measure_tree  <- to_vec(for(i in 1:2) if(TRUE) 2*recalls_tree[i]*precisions_tree[i]/(recalls_tree[i]+precisions_tree[i])  )
  F1_tree <- mean(c(F1_measure_tree))
  
  
  
  tt_tree_aux <- tt_tree_aux + tt_tree
  acc_tree_aux <- acc_tree_aux + acc_tree
  recalls_tree_aux <- recalls_tree_aux + recalls_tree
  precisions_tree_aux  <- precisions_tree_aux + precisions_tree
  F1_measure_tree_aux  <- F1_measure_tree_aux + F1_measure_tree
  F1_tree_aux <- F1_tree_aux + F1_tree
  
}

rpart.plot(fit_tree,box.palette = colorRampPalette(c("#F55118","#FF9933","#FFFFFF","#66CCCC","#298DF7"))(50))


tt_tree_final <- tt_tree_aux/20;tt_tree_final
acc_tree_final <- acc_tree_aux/20;acc_tree_final
recalls_tree_final <- recalls_tree_aux/20;recalls_tree_final
precisions_tree_final  <- precisions_tree_aux/20;precisions_tree_final
F1_measure_tree_final  <- F1_measure_tree_aux/20;F1_measure_tree_final
F1_tree_final <- F1_tree_aux/20;F1_tree_final

# Min-max normalize the data for the neural network analyses####


minmax <- function(x){(x-min(x))/(max(x)-min(x))}

data_corr_pca2 <- data_corr_pca


data_corr_pca2$PC1 <- minmax(data_corr_pca2$PC1)
data_corr_pca2$PC2 <- minmax(data_corr_pca2$PC2)
data_corr_pca2$PC3 <- minmax(data_corr_pca2$PC3)
data_corr_pca2$PC4 <- minmax(data_corr_pca2$PC4)
data_corr_pca2$PC5 <- minmax(data_corr_pca2$PC5)
data_corr_pca2$PC6 <- minmax(data_corr_pca2$PC6)
data_corr_pca2$PC7 <- minmax(data_corr_pca2$PC7)


# Neural Network##########################################################################################################
library(neuralnet)

tt_nn_aux <- table(c(0,1), c(0,1)) - table(c(0,1), c(0,1))
acc_nn_aux <- 0
recalls_nn_aux <- c(0,0)
precisions_nn_aux  <- c(0,0)
F1_measure_nn_aux  <- c(0,0)
F1_nn_aux <- 0

for (i in 1:20) {
  
  
  spl = sample.split(data_corr_pca2, SplitRatio = 0.9)
  train = subset(data_corr_pca2, spl==TRUE)
  train_clusters = factor(subset(k1$cluster, spl==TRUE))
  test = subset(data_corr_pca2, spl==FALSE)
  test_clusters = factor(subset(k1$cluster, spl==FALSE))
  
  
  fit_nn <- neuralnet(train_clusters~., data=train[,-c(1,2,3,4)], hidden=c(5,2,4), linear.output=FALSE, threshold=0.01)
  
  predictions_nn <- predict(fit_nn,  test[,-c(1,2,3,4)])
  
  predictions_nn2 <- to_vec(for(i in 1:nrow(predictions_nn)) if(TRUE) which.max(predictions_nn[i,]) )
  
  
  tt_nn <- table(predictions_nn2, test_clusters)
  acc_nn <- sum(diag(tt_nn))/sum(tt_nn)
  recalls_nn <- to_vec(for(i in 1:2) if(TRUE) diag(tt_nn)[i]/sum(tt_nn[,i]))
  precisions_nn  <- to_vec(for(i in 1:2) if(TRUE) diag(tt_nn)[i]/sum(tt_nn[i,]))
  F1_measure_nn  <- to_vec(for(i in 1:2) if(TRUE) 2*recalls_nn[i]*precisions_nn[i]/(recalls_nn[i]+precisions_nn[i])  )
  F1_nn <- mean(c(F1_measure_nn))
  
  
  
  tt_nn_aux <- tt_nn_aux + tt_nn
  acc_nn_aux <- acc_nn_aux + acc_nn
  recalls_nn_aux <- recalls_nn_aux + recalls_nn
  precisions_nn_aux  <- precisions_nn_aux + precisions_nn
  F1_measure_nn_aux  <- F1_measure_nn_aux + F1_measure_nn
  F1_nn_aux <- F1_nn_aux + F1_nn
  
}

plot(fit_nn)

tt_nn_final <- tt_nn_aux/20;tt_nn_final
acc_nn_final <- acc_nn_aux/20;acc_nn_final
recalls_nn_final <- recalls_nn_aux/20;recalls_nn_final
precisions_nn_final  <- precisions_nn_aux/20;precisions_nn_final
F1_measure_nn_final  <- F1_measure_nn_aux/20;F1_measure_nn_final
F1_nn_final <- F1_nn_aux/20;F1_nn_final