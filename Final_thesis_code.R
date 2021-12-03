##Load libraries - use install.packages() if these libraries are not installed.
library(mice)
library("heplots")
library(effectsize)
library(umapr)
library(tidyverse)
library(gplots)
library("circlize")
library(parallel)
library(clustertend)
library("psych")
library("ggplot2")
library('Rtsne')
library(dplyr)
library("gplots")
library("factoextra")
library(clValid)
library(NbClust)
library("ComplexHeatmap")
library("fpc")
library("dbscan")
library("rgl")
library(GGally) 
library('kohonen') 
library('Rtsne')
library('rgl')
library('Seurat')
library("kohonen")
library('Rtsne')
library('rgl')
library('Seurat')
library("pivottabler")
library(mclust)
library(hrbrthemes)
library(viridis)
library(plotly)
library(scatterplot3d)
library(plyr)
library(data.table)
library(vtree)

###Load the dataset
###In this section the data for clustering is loaded 
memory.limit() #due to problems with vector sizes this function was used to check the memory limit
memory.limit(100000) #the memory limit is adjusted to 1300000 to enable all functions to run
save(object = list,file="data")
load.Rdata(filename = "data", object =list )
gc() #the garbage collection function can be used to make memory use more efficient
dataDir  <- "../Thesis_D/"
data1 <- read.csv2(paste0(dataDir, "Experimental_finaltestdata.csv"), stringsAsFactors = TRUE)#the data is loaded in the EU format with strings as factors.
summary(data1) #Check the data is loaded, observable and loads in the correct format.

##--Data Cleaning-------------------------------------------------------------##
###Small changes are made to the data set to clean up any problems found during an observation of the summary of the data
###From observing the summary of the data some of the factors had incorrect levels and this was corrected
levels(data1$Result)[levels(data1$Result)=="win"] <- "Win"
levels(data1$Starter_Sub)[levels(data1$Starter_Sub)=="sub"] <- "Sub"
levels(data1$Starter_Sub)[levels(data1$Starter_Sub)=="Sub "] <- "Sub"
levels(data1$Half)[levels(data1$Half)==" First"] <- "First"
levels(data1$Half)[levels(data1$Half)==" First "] <- "First"
levels(data1$Half)[levels(data1$Half)==" First  "] <- "First"
levels(data1$Half)[levels(data1$Half)==" First"] <- "First"
levels(data1$Half)[levels(data1$Half)=="first"] <- "First"
levels(data1$Half)[levels(data1$Half)=="Seond"] <- "Second"
levels(data1$Half)[levels(data1$Half)==" Second"] <- "Second"
levels(data1$Half)[levels(data1$Half)==" Second "] <- "Second"
summary(data1)

##--Exploratory data analysis-------------------------------------------------------------##
### Exploring the dataset, distribution of variables and correlations 
```{r task 1}
summary(data1)
describe(data1) #this function from the psych package is useful for a quick scan of the numerical variables in the dataset

### Vtree propvides a digram showing the count and/or percent of splits in the data by the assigned variables
vtree(model_mem_final, c("Team", "Position"), 
      fillcolor = c(Team = "#42B540FF", Position = "#0099B4FF"),
      horiz = FALSE,
      showpct = FALSE)

#a visualization of the the distribution of the Total Distance Variable 
hist(data1$TD, breaks = 50, xlab = "Total Distance",xlim = c(0,250), main = "Distribution of Total Distance Varaible")

# a visualisation of the number of observations for each team by position, the colourcodes are from the lancet palette
ggplot(data1, aes(x = Team, fill = Position))+
  geom_bar(position = "dodge")+
  scale_fill_manual("Position", values = c("AM"= "#00468BFF","CB"= "#ED0000FF","DM"= "#42B540FF","FB"= "#0099B4FF","ST"=  "#925E9FFF","W"= "#FDAF91FF"))+
  geom_text(aes(label = ..count..),stat = "count", vjust = -0.2, size = 3, position = position_dodge(.9))

###Feature Extraction
###My dataset contained variables that are fully nested in other variables and need to be separated.  An example is that sprint distance (SD) high speed distance and running distance are all part of total distance.

data1$otherdist <- data1$TD - data1$RD #I create a feature for all distance covered at a speed slower than running distance <14.5km/h
data1$rundist <-data1$RD - data1$HSD #I create a feature for all distance covered between 14.5km/h and 19.8km/h
data1$hsdist <- data1$HSD - data1$SPR #I create a feature for all distance covered between 19.8km/h and 25.2km/h
data1$OHSR <- data1$HSR-data1$SPR #I create a feature of only high speed runs and remove the sprints

#The summary and describe functions can be used again to check the new features have been created correctly
summary(data1) #I use the summary function to check the new features have been created correctly
describe(data1[c(4:7,9,20:23)])

features <- data1[c(4:7,9,20:23)] #Create a smaller data frame of the variables to be used as features
model_data <- scale(features)#Features are scaled into a matrix to be used in the model training
head(model_data) #Check the data is scaled correctly

ggpairs(features) #ggpairs allows us to observe the pairwise correlations between our features

hopkins(model_data, nrow(model_data)-1 ) #the hopkins function allows us to assess the clustering tendency of our data set.  The maximum size is number of rows-1

#Visual Assessment of cluster Tendency is another method that enables us to view the clustering tendency of the dataset.  This function is computationally expensive.  Therefore splitting the dataset might be required.  Red colour corresponds to small distance and blue colour indicates big distance between observation.

set.seed(123) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 50% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(model_data), size = floor(.50*nrow(model_data)), replace = F)
VAT <- model_data[sample, ]
#finally we can create a plot of the VAT
fviz_dist(dist(VAT), show_labels = FALSE) + labs(title = "Visual Assessment of clustering Tendency")

#potential clustering can also be viewed using a heatmap.  This colors the scaled columns based on their observations by row
heatmap(model_data, scale = "row")

###-Model Selection-------------------------------------------------------------##
###Selecting the ideal candidate models
###Here you can use the clValid function from the clValid package 
```{r task 1}
gc()
time <-proc.time()  #start a timer
model <- c("kmeans") # there are more options that can be specified
k <- 2:8 #specifiy the number of clusters you want to test
scores<- clValid(model_data, nClust = k, clMethods = model, metric = "euclidean", validation = "internal", maxitems=nrow(model_data))
summary(scores)
#note that the metric and validation can also be changed
proc.time() - time #stop the timer

##--Hyperparameter Tuning-------------------------------------------------------------##
### Selecting the ideal candidate models
### Here you can use the clValid function from the clValid package 
### we can determine the ideal number of clusters using the elbow, average silouhette and gap statistic
elbow2 <- fviz_nbclust(model_data, kmeans, iter.max = 200, method = "wss") +
  labs(subtitle = "Elbow method")
elbow2

# we can determine the ideal number of clusters using the elbow, average silouhette and gap statistic
silhouette <- fviz_nbclust(model_data, kmeans, iter.max = 200, method = "silhouette")+
  labs(subtitle = "Silhoutte method")
silhouette

# we can determine the ideal number of clusters using the elbow, average silouhette and gap statistic

set.seed(123)
gap <- fviz_nbclust(model_data, kmeans, nstart =25, iter.max = 500, method = "gap_stat", nboot = 50) + labs(subtitle = "Gap statistic method")
gap

# we can determine the ideal number of clusters through prediction strength by assesseing the prediction strength in terms of a cutoff
set.seed(123)
pred_strength <- prediction.strength(model_data, Gmin=2, Gmax=10, M=50,
                                     clustermethod=kmeansCBI,
                                     classification="centroid", centroidname = NULL,
                                     cutoff=0.8,nnk=1,
                                     distances=inherits(model_data,"dist"),count=FALSE)

# we can plot a graph to usualise these changes in clusters
plot(seq(2,10),pred_strength$mean.pred[-1],main="Prediction strength method", col=c('steelblue'), type="l", xlab="Number of clusters", ylab="Prediction strength")
text(seq(2,10),pred_strength$mean.pred[-1],labels = round(pred_strength$mean.pred[-1],2), cex = 0.6, pos = 1, col = "black" )

#run the kmeans model with k set to 6 on the dataset 
set.seed(123)
km6 <- eclust(model_data, 
              "kmeans", 
              k = 6, 
              nstart = 50, #attempts n initial configurations of random centroids and reports on the best one. 
              graph = FALSE )
print(km6)

#view the kmeans model cluster means
km6mean <- as.data.frame(aggregate(data1, by=list(km6$cluster), mean))
km6mean

#we can visualise the clusters using the fviz_cluster function
fviz_cluster(list(data = model_test, cluster = km6$cluster), ellipse.type = "norm", geom = "point", stand = FALSE, palette = "lancet", ggtheme = theme_classic())

# we can visualise silhouettes by cluster using the fviz_silhouette function
fviz_silhouette(km6, palette = "lancet", ggtheme = theme_classic(), )

###The second model is the density based model 
##For a baseline I am running the HDBSCAN
set.seed(123)
hdb <- hdbscan(model_data, minPts = 18) # min points = number of dimensions *2
hdb

###The third model 
###We inspect the k-NN distance plot for k = minPts - 1 = 17
set.seed(123)
dbscan::kNNdistplot(model_data, k = 17)
abline(h = 3, col = "steelblue", lty =2) # the h can specified as a plot to help you visualise where the steep vertical increase occurs

###When we know the ideal eps value we can fit the DBSCAN model
set.seed(123)
dbscan <- fpc::dbscan(model_data, eps = 3, MinPts = 18)
dbscan

###The fourth model 
###Agnes using euclidean distance with k set to 6
agnes_euc<- eclust(model_data, "agnes", k = 6, hc_metric = "euclidean")
aggregate(data1, by=list(agnes_euc$cluster), mean)

mc3 <-fviz_mclust(mc, "uncertainty", palette = "lancet")

set.seed(123)
hkmeans <- hkmeans(model_data, 6)
mc <- Mclust(model_data)
summary(mc)
#BIC values for number of clusters
mc_bic <-fviz_mclust(mc, "BIC", palette = "lancet")

#Classification plot that includes uncertainty
mc_unc <-fviz_mclust(mc, "uncertainty", palette = "lancet")

###The fifth model 
###hkmeans with k set to 6
set.seed(123)
hkmeans <- hkmeans(model_data, 6)


#compute the disimilarity matrix
res.dist <- dist(model_data, method = "manhattan")
as.matrix(res.dist)[1:9, 1:9]
```

create the agnes manhattan model
hc <- hclust(res.dist, method = "ward.D2")
summary(hc)


#we can use an iterative process to build up our final dataset with all models cluster assignments .  This process was selecetd because it enables us to easily go back and fix errors.
gc()
model_mem<- cbind(data1, cluster_k6 = km6$cluster)
model_mem2<- cbind(model_mem, cluster_hc = hclustdf$cluster)
model_mem3<- cbind(model_mem2, cluster_ag = agnes_euc$cluster)
model_mem4<- cbind(model_mem3, cluster_hk = hkmeans$cluster)
model_mem5<- cbind(model_mem4, cluster_mc = mc$classification)
model_mem_final <- cbind(model_mem5, cluster_db = dbscan$cluster)
head(model_mem_final[c(24:27)]) #check all the clusters are in the dataframe
qpvt(model_mem_final, "cluster_db", "cluster_ag", "n()") # this generates a pivot table that can compare the cluster assignments

#this code creates the final multivariate match periods by looping through the cluster assignmnets and determining which MMP to allocate an observation to
model_mem_final$finalclusters <- 0
model_mem_final$cluster_db
model_mem_final$cluster_db <- as.numeric(model_mem_final$cluster_db)
model_mem_final$cluster_hk <- as.numeric(model_mem_final$cluster_hk)
for (i in 1:dim(model_mem_final)[1]) {
  if (model_mem_final$cluster_db[i] == 0 & model_mem_final$cluster_hk[i] == 1) {
    model_mem_final$finalclusters[i] <- as.numeric(1)
  } else if (model_mem_final$cluster_db[i] == 0 & model_mem_final$cluster_hk[i] == 2) {
    model_mem_final$finalclusters[i] <- as.numeric(1)
  } else if (model_mem_final$cluster_db[i] == 0 & model_mem_final$cluster_hk[i] == 3) {
    model_mem_final$finalclusters[i] <- as.numeric(1)
  } else if (model_mem_final$cluster_db[i] == 0 & model_mem_final$cluster_hk[i] == 4) {
    model_mem_final$finalclusters[i] <- as.numeric(1)
  } else if (model_mem_final$cluster_db[i] == 0 & model_mem_final$cluster_hk[i] == 5) {
    model_mem_final$finalclusters[i] <- as.numeric(1)
  } else if (model_mem_final$cluster_db[i] == 0 & model_mem_final$cluster_hk[i] == 6) {
    model_mem_final$finalclusters[i] <- as.numeric(1)
  } else if (model_mem_final$cluster_db[i] == 3 & model_mem_final$cluster_hk[i] == 1) {
    model_mem_final$finalclusters[i] <- as.numeric(2)
  } else if (model_mem_final$cluster_db[i] == 3 & model_mem_final$cluster_hk[i] == 2) {
    model_mem_final$finalclusters[i] <- as.numeric(2)
  } else if (model_mem_final$cluster_db[i] == 3 & model_mem_final$cluster_hk[i] == 3) {
    model_mem_final$finalclusters[i] <- as.numeric(2)
  } else if (model_mem_final$cluster_db[i] == 3 & model_mem_final$cluster_hk[i] == 4) {
    model_mem_final$finalclusters[i] <- as.numeric(2)
  } else if (model_mem_final$cluster_db[i] == 3 & model_mem_final$cluster_hk[i] == 5) {
    model_mem_final$finalclusters[i] <- as.numeric(2)
  } else if (model_mem_final$cluster_db[i] == 3 & model_mem_final$cluster_hk[i] == 6) {
    model_mem_final$finalclusters[i] <- as.numeric(2)
  } else if (model_mem_final$cluster_db[i] == 1 & model_mem_final$cluster_hk[i] == 4) {
    model_mem_final$finalclusters[i] <- as.numeric(3)
  } else if (model_mem_final$cluster_db[i] == 1 & model_mem_final$cluster_hk[i] == 6) {
    model_mem_final$finalclusters[i] <- as.numeric(4)
  } else if (model_mem_final$cluster_db[i] == 1 & model_mem_final$cluster_hk[i] == 2) {
    model_mem_final$finalclusters[i] <- as.numeric(5)
  } else if (model_mem_final$cluster_db[i] == 1 & model_mem_final$cluster_hk[i] == 3) {
    model_mem_final$finalclusters[i] <- as.numeric(6)
  } else if (model_mem_final$cluster_db[i] == 1 & model_mem_final$cluster_hk[i] == 5) {
    model_mem_final$finalclusters[i] <- as.numeric(7)
  } else if (model_mem_final$cluster_db[i] == 1 & model_mem_final$cluster_hk[i] == 1) {
    model_mem_final$finalclusters[i] <- as.numeric(8)
  } else if (model_mem_final$cluster_db[i] == 2 & model_mem_final$cluster_hk[i] == 4) {
    model_mem_final$finalclusters[i] <- as.numeric(3)
  } else if (model_mem_final$cluster_db[i] == 2 & model_mem_final$cluster_hk[i] == 6) {
    model_mem_final$finalclusters[i] <- as.numeric(4)
  } else if (model_mem_final$cluster_db[i] == 2 & model_mem_final$cluster_hk[i] == 2) {
    model_mem_final$finalclusters[i] <- as.numeric(5)
  } else if (model_mem_final$cluster_db[i] == 2 & model_mem_final$cluster_hk[i] == 3) {
    model_mem_final$finalclusters[i] <- as.numeric(6)
  } else if (model_mem_final$cluster_db[i] == 2 & model_mem_final$cluster_hk[i] == 5) {
    model_mem_final$finalclusters[i] <- as.numeric(7)
  } else if (model_mem_final$cluster_db[i] == 2 & model_mem_final$cluster_hk[i] == 1) {
    model_mem_final$finalclusters[i] <- as.numeric(8)
  }
}
model_mem_final$finalclusters
model_mem_final$finalclusters <- as.factor(model_mem_final$finalclusters)

aggregate(model_mem_final, by=list(model_mem_final$finalclusters), mean)

#create a tsne representation of the dataset
tsne <- Rtsne(model_data, check_duplicates = FALSE, pca = TRUE, perplexity=50, theta=0.5, dims=3)
#display results of t-SNE
tsne_out <- Rtsne(model_data,dims=3) # Run TSNE

colors1 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF","#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF") #lancet palette
colors <- colors1[as.numeric(model_mem_final$finalclusters)] # allocate colours to the MMPs
plot(tsne_out$Y,col= colors)# Plot the result


# this generates a 3d image of the tsne
scat_km <- scatterplot3d(x=tsne_out$Y[,1],y=tsne_out$Y[,2],z=tsne_out$Y[,3], color = colors)
legend("top", legend = levels(as.factor(model_mem_final$finalclusters)),
       col =  c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF","#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF"), 
       pch = c(16, 17, 18), 
       inset = -0.25, xpd = TRUE, horiz = TRUE)# Plot using scatterplot3d package

# this code opens a 3d plot in an external tab and can be rotated
plot3d(x=tsne_out$Y[,1],y=tsne_out$Y[,2],z=tsne_out$Y[,3], col = c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF","#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF")[model_mem_final$finalclusters])
legend("bottom", legend = levels(as.factor(model_mem_final$finalclusters)),
       col =  c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF","#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF"), 
       pch = c(16, 17, 18), 
       inset = -0.25, xpd = TRUE, horiz = TRUE)

# this code was used to annonymise the data
model_mem_final$player <- model_mem_final$Name
levels(model_mem_final$player)<- c("Player 01", "Player 02", "Player 03", "Player 04", "Player 05", "Player 06", "Player 07", "Player 08", "Player 09", "Player 10", "Player 11", "Player 12", "Player 13", "Player 14", "Player 15", "Player 16", "Player 17", "Player 18", "Player 19", "Player 20", "Player 21", "Player 22", "Player 23", "Player 24", "Player 25", "Player 26", "Player 27", "Player 28", "Player 29", "Player 30", "Player 31", "Player 32", "Player 33", "Player 34", "Player 35", "Player 36", "Player 37", "Player 38", "Player 39", "Player 40", "Player 41", "Player 42", "Player 43", "Player 44", "Player 45", "Player 46")

# this code generates the bar plot used to compare the observations of MMPs by position
model_mem_final$Position <- factor(model_mem_final$Position,levels = c("CB", "FB", "W", "ST", "AM", "DM"))
ggplot(data = model_mem_final)+
  geom_bar(
    mapping = aes(x = Position, fill = finalclusters),
    position = "fill")+
  scale_fill_manual("MMP", values = c("1"= "#00468BFF","2"= "#ED0000FF","3"= "#42B540FF","4"= "#0099B4FF","5"=  "#925E9FFF","6"= "#FDAF91FF","7"= "#AD002AFF","8"= "#ADB6B6FF"))

#we can also find the cluster stats for any model using for example the cluster.stats() function
hk_stats <- cluster.stats(dist(model_data),  hkmeans$cluster)

# this code enables us to generate a umap visualisation of the clusters
embedding <- umap(model_data)
head(embedding)

colorumap <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF91FF") #additional colors if required "#AD002AFF", "#ADB6B6FF", "#1B1919FF"
embedding %>% 
  mutate(Clusters = agnes_euc$cluster) %>%
  ggplot(aes(UMAP1, UMAP2, color = Clusters)) + geom_point()+
  scale_color_gradientn(colours = colorumap)

# this code generates this an image used to compare the MMP's with TD observations >80% of the max of TD for each player
model_mem_final$MMP <- model_mem_final$finalclusters
model_graph <- ddply(model_mem_final, "player", subset, TD >= 0.8 * max(TD))
model_graph
ggplot(data= model_graph)+
  geom_point(
    mapping=aes(x = reorder(player, TD, FUN = max), y = TD, color = MMP, size=1)
  )+
  theme_ipsum()+
  scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6","7","8"),
                     values=c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF",  "#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF"))+
  coord_flip()

getwd()
write.csv(model_mem_final,"C:/Users/rowan/Desktop/Thesis/Thesis_R/finalcorrected.csv", row.names = FALSE)

# this is the code to test the variance in the MMPs using MANOVA 
manova_analysis <- manova(cbind(TD, RD, HSD, SD, PL, ACC, DEC, HSR, SPR)~finalclusters, data = model_mem_final)
summary(manova_analysis)
summary.aov(manova_analysis)

#this is the code that can be used to calculate the eta sq or wilks lambda
effectsize::eta_squared(manova_analysis)
test4effect<- lm(cbind(TD, HSD, RD, SD, PL, ACC, DEC, HSR)~Half, data = model_mem_final)
etasq(Anova(test4effect), anova=TRUE)
etasq(test4effect, test="Wilks")

# this code was used to impute the missing data 
youth <- subset(model_mem_final, Team == "Youth")
youth$HR[youth$HR == 0] <- NA
md.pattern(youth)
imputed_Data2 <- mice(youth, m=5, maxit = 50, method = 'pmm', seed = 123)
summary(imputed_Data2)
imputed_Data2
imputed_Data2$imp$HR
completeData <- mice::complete(imputed_Data2,5)
completeData
describe.by(completeData, completeData$finalclusters)
out0 <- lm(HR ~ 1, data = completeData)
out1 <- lm(HR ~ finalclusters, data = completeData)
out2 <- lm(HR ~ finalclusters + player, data = completeData)
out3 <- lm(HR ~ player, data = completeData)

# this is the code used 
cluster1data <- subset(model_mem_final, model_mem_final$finalclusters ==1 )
cluster2data <- subset(model_mem_final, model_mem_final$finalclusters ==2 )
cluster3data <- subset(model_mem_final, model_mem_final$finalclusters ==3 )
cluster4data <- subset(model_mem_final, model_mem_final$finalclusters ==4 )
cluster5data <- subset(model_mem_final, model_mem_final$finalclusters ==5 )
cluster6data <- subset(model_mem_final, model_mem_final$finalclusters ==6 )
cluster7data <- subset(model_mem_final, model_mem_final$finalclusters ==7 )
cluster8data <- subset(model_mem_final, model_mem_final$finalclusters ==8 )
manova.cluster <- manova(cbind(TD, HSD, RD, SD, PL, ACC, DEC, HSR)~Team, data = cluster6data)
summary(manova.cluster)
summary.aov(manova.cluster) # exam the relationship between the dependent variables

#we can finally write the final data,frame to a csv
getwd()
write.csv(model_mem_final,"C:working directory", row.names = FALSE)




