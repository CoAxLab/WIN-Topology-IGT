
library(gridExtra)
library(grid)
library(outliers)
library(psych)

mergedWINData = read.csv(file="mergedWINData.csv")
mergedWINData$X = NULL

binaryVariables <- mergedWINData[,c(2,3,4,6,8,10,16,18,46,47)]

myWINDataAnalysisFrame <- t(do.call(cbind, lapply(binaryVariables, summary)))

myWINDataAnalysisFrame <- myWINDataAnalysisFrame[,-7] # removing last column (NA's count)

# skew(mergedWINData$density_baseline)
# skew(mergedWINData$clustering_coeff_average.binary._baseline)
# skew(mergedWINData$transitivity.binary._baseline)
# skew(mergedWINData$network_characteristic_path_length.binary._baseline)
# skew(mergedWINData$small.worldness.binary._baseline)
# skew(mergedWINData$global_efficiency.binary._baseline)
# skew(mergedWINData$local_efficiency.binary._baseline)
# skew(mergedWINData$assortativity_coefficient.binary._baseline)
# skew(mergedWINData$baseline_p)
# skew(mergedWINData$baseline_q)

skewValues <- c(0.02684528,0.1126648,1.073562, 1.237727,0.01113487,
                -0.3591412,0.07801355,0.1502269,-0.2112031,-0.779505)

myWINDataAnalysisFrame <- cbind(myWINDataAnalysisFrame,Skew = skewValues)

# grubbs.test(mergedWINData$density_baseline)
# grubbs.test(mergedWINData$clustering_coeff_average.binary._baseline)
# grubbs.test(mergedWINData$transitivity.binary._baseline)
# grubbs.test(mergedWINData$network_characteristic_path_length.binary._baseline)
# grubbs.test(mergedWINData$small.worldness.binary._baseline)
# grubbs.test(mergedWINData$global_efficiency.binary._baseline)
# grubbs.test(mergedWINData$local_efficiency.binary._baseline)
# grubbs.test(mergedWINData$assortativity_coefficient.binary._baseline)
# grubbs.test(mergedWINData$baseline_p)
# grubbs.test(mergedWINData$baseline_q)

grubbsTestOutliers <- c(1,1,1,1,1,1,1,1,1,1)

myWINDataAnalysisFrame <- cbind(myWINDataAnalysisFrame, Grubbs = grubbsTestOutliers)

colnames(myWINDataAnalysisFrame) <- c("Minimum", "Quartile 1", "Median", "Mean", 
                                      "Quartile 3", "Maximum", "Skew", "Grubbs")

rownames(myWINDataAnalysisFrame) <- c("Density", "Clustering Coefficient", "Transitivity", 
                                      "Net. Char. Path Length","Small Worldness", 
                                      "Global Efficiency", "Local Efficiency", 
                                      "Assortativity", "P Scores", "Q Scores")



myDataFrame <- data.frame(round(myWINDataAnalysisFrame, digits = 3))

myTheme <- ttheme_default(base_size = 8, base_colour = "black")

grid.table(myDataFrame, theme = myTheme)

grid.draw(textGrob("Distribution of Binary Structural Topology Measures", 
                   gp = gpar(fontsize = 13), x = unit(0.6, "npc"), y = unit(0.71, "npc")))



library(reshape2)
library(ggplot2)

#Selecting Binary Structural Topology Measures
myData <- mergedWINData[1:126, c(2,3,4,6,8,10,16,18)]

#Making a Correlation Matrix
corrMatrix <- round(cor(myData, use="complete.obs"), 2)

#Only keeping the lower half of the matrix
getLowerTri <- function(corrMatrix){
  corrMatrix[upper.tri(corrMatrix)] <- NA
  return(corrMatrix)}
lowerTri <- getLowerTri(corrMatrix)
meltedCorrMatrix <- melt(lowerTri, na.rm = TRUE )

#Plotting the Correlation Matrix and the plot characteristics
ggplot(data = meltedCorrMatrix, p.mat = p.mat, aes(x=Var1, y=Var2, fill=value)) + 
  
  geom_tile(color="white") +
  
  scale_x_discrete(labels = c("Density", "Clustering Coeff", "Transitivity", "Char Path Length", 
                            "Small Worldness", "Global Efficiency", "Local Efficiency", 
                            "Assortativity")) +
  
  scale_y_discrete(labels = c("Density", "Clustering Coeff", "Transitivity", "Char Path Length", 
                            "Small Worldness", "Global Efficiency", "Local Efficiency", 
                            "Assortativity")) +
  
  scale_fill_gradient(low = "gold", high = "red", 
                      limit = c(-1,1), space = "Lab",
                      name = "WIN Binary Variables\nCorrelation Matrix\n\n") + 
  theme_minimal() + 
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1)) +
  
  coord_fixed() +
  
  geom_text(aes(Var1, Var2, label = value), color = "white", size = 4)


#Standardization of Structural Topology Measures
density_baseline_scaled <- scale(mergedWINData$density_baseline)
clustering_coeff_average.binary._baseline_scaled <- scale(mergedWINData$clustering_coeff_average.binary._baseline)
transitivity.binary._baseline_scaled <- scale(mergedWINData$transitivity.binary._baseline)
network_characteristic_path_length.binary._baseline_scaled <- scale(mergedWINData$network_characteristic_path_length.binary._baseline)
small.worldness.binary._baseline_scaled <- scale(mergedWINData$small.worldness.binary._baseline)
global_efficiency.binary._baseline_scaled <- scale(mergedWINData$global_efficiency.binary._baseline)
local_efficiency.binary._baseline_scaled <- scale(mergedWINData$local_efficiency.binary._baseline)
assortativity_coefficient.binary._baseline_scaled <- scale(mergedWINData$assortativity_coefficient.binary._baseline)


#Creating Data Frame of Standardized Variables
standardizedVariables <- data.frame(density=density_baseline_scaled,
                                    clust_binary=clustering_coeff_average.binary._baseline_scaled,
                                    trans_binary=transitivity.binary._baseline_scaled,
                                    net_binary=network_characteristic_path_length.binary._baseline_scaled,
                                    small_binary=small.worldness.binary._baseline_scaled,
                                    global_binary=global_efficiency.binary._baseline_scaled,
                                    local_binary=local_efficiency.binary._baseline_scaled,
                                    assort_binary=assortativity_coefficient.binary._baseline_scaled)



library("factoextra")

myData <- standardizedVariables

myPCA = princomp(na.omit(myData),
               cor = TRUE,
               scores = TRUE)

summary(myPCA)

myValues <- c(0.421659, 0.6716133, 0.8111032, 0.9232412, 0.99102908, 0.997264861, 0.99961437, 1.0000000000)
myComponents <- c(1,2,3,4,5,6,7,8)

myData <- data.frame(myComponents, myValues)

ggplot(myData, aes(x=myComponents, y=myValues)) +
  geom_bar(stat="identity", fill="steelblue") + 
    geom_text(aes(label=round(myValues, digits=5)), vjust=1.6, color="white", size=3.5)+
        theme_minimal() +
            labs(x="Components", y = "Cumulative Proportion") +
                ggtitle("Cumulative Proportion of Variance Explained") +
                    theme(plot.title = element_text(hjust = 0.5)) +
                        geom_hline(yintercept=.95,linetype="dashed", color = "red")

kfcv.sizes = function(n, k=10) {
  # generate sample sizes for k-fold cross validation on a data set of
  # size n
  
  # author: Matthias C. M. Troffaes
  # date: 22 Nov 2010
  # license: GPLv3
  # usage:
  #
  #   kfcv.sizes(n, k=...)
  #

  sizes = c()
  for (i in 1:k) {
    first = 1 + (((i - 1) * n) %/% k)
    last = ((i * n) %/% k)
    sizes = append(sizes, last - first + 1)
  }
  sizes
}

kfcv.testing = function(n, k=10) {
  # generate testing sample indices for k-fold cross validation on a
  # data set of size n
  
  # author: Matthias C. M. Troffaes
  # date: 22 Nov 2010
  # license: GPLv3
  # usage:
  #
  #   kfcv.testing(n, k=...)
  #

  indices = list()
  sizes = kfcv.sizes(n, k=k)
  values = 1:n
  for (i in 1:k) {
    # take a random sample of given size
    s = sample(values, sizes[i])
    # append random sample to list of indices
    indices[[i]] = s
    # remove sample from values
    values = setdiff(values, s)
  }
  indices
}


library(MASS)
library(boot)
library(glmnet)

#Principal Component Calculation
myData <- standardizedVariables[,c(1,2,3,4,5,6,7,8)]
myPCA <- princomp(na.omit(myData), cor = TRUE, scores = TRUE)
myComponents <- myPCA$scores[, 1:5]

#Fixing P Score Length
myP <- mergedWINData$baseline_p[setdiff(1:126,c(118,126))]


#Combining both lists into one data.frame
myComboData <- data.frame(myPCA$scores[, 1:5], myP)

#GLM Calculation
myFit <- glm(myP~myComponents, family= gaussian, data=myComboData)

myPred <- myFit$fitted.values

#Plotting Predicted versus Observed P Scores
plot(myP, myPred, main="Predicted vs Observed")
abline(lm(myPred~myP), col="red")
cor(myP, myPred)


library(MASS)
library(boot)

#Principal Component Calculation
myData <- standardizedVariables[,c(1,2,3,4,5,6,7,8)]
myPCA <- princomp(na.omit(myData), cor = TRUE, scores = TRUE)
myComponents <- myPCA$scores[, 1:5]

#Fixing Q Score Length
myQ <- mergedWINData$baseline_q[setdiff(1:126,c(118,126))]


#Combining both lists into one data.frame
myComboData <- data.frame(myPCA$scores[, 1:5], myQ)

#GLM Calculation
myFit <- glm(myQ~myComponents, family= gaussian, data=myComboData)

#Prediction values from GLM
myPred <- myFit$fitted.values

#Plotting Predicted versus Observed Q Scores
plot(myQ, myPred, main="Predicted vs Observed")
abline(lm(myPred~myQ), col="red")
cor(myQ, myPred)


#Principal Component Calculation

# Turn off warnings
options(warn=-1)

# Define number of folds
n_folds = 124

# Define number of resamples
n_iter = 10

store_P_predictions = matrix(0,124,n_iter)
store_P_coef = matrix(0,n_iter, 6)
iter_P_corr = matrix(0,n_iter,1)

for (iter in 1:n_iter) {
    
    # Make our 10 fold cross validation indices
    indices = kfcv.testing(124, k=n_folds)
    
    PredP = matrix(0,124,1)
    CVPCoefs = matrix(0,n_folds,6)

    for (fold in 1:n_folds) {
        test = unlist(indices[fold])
        train = unlist(indices[setdiff(1:n_folds, fold)])

        cvfit <- glm(myP~myComponents, data=myComboData, subset=train)
        PredP[test] = predict(cvfit, data=myComboData, subset=test)
        CVPCoefs[fold,] = coef(cvfit)

    }
    
    # Storing P score Predictions
    iter_P_corr[iter] = cor(myP, PredP)
    store_P_predictions[,iter] = PredP
    store_P_coef[iter,] = colMeans(CVPCoefs)
}

# Plotting a histogram of all the correlation coefficients made
hist(iter_P_corr)



#Principal Component Calculation

store_Q_predictions = matrix(0,124,n_iter)
store_Q_coef = matrix(0,n_iter, 6)
iter_Q_corr = matrix(0,n_iter,1)

for (iter in 1:n_iter) {
    
    # Make our k fold cross validation indices
    indices = kfcv.testing(124, k=n_folds)
    
    PredQ = matrix(0,124,1)
    CVQCoefs = matrix(0,n_folds,6)

    for (fold in 1:n_folds) {
        test = unlist(indices[fold])
        train = unlist(indices[setdiff(1:n_folds, fold)])

        cvfit <- glm(myQ~myComponents, family=gaussian, data=myComboData, subset=train)
        PredQ[test] = predict(cvfit, data=myComboData, subset=test)
        CVQCoefs[fold,] = coef(cvfit)

    }
    
    # Storing P score Predictions
    iter_Q_corr[iter] = cor(myQ, PredQ)
    store_Q_predictions[,iter] = PredQ
    store_Q_coef[iter,] = colMeans(CVQCoefs)
}

# Plotting a histogram of all the correlation coefficients made
hist(iter_Q_corr)




PModelCoefs = colMeans(store_P_coef)
PModelCoefs

# Get the matrix to transform back to data space
V = myPCA$loadings
V

# The weight (or influence) of each term on P
wP = V %*% c(PModelCoefs[2:6],0,0,0)
wP

QModelCoefs = colMeans(store_Q_coef)

# The weight (or influence) of each term on Q
wQ = V %*% c(QModelCoefs[2:6],0,0,0)
wQ

# Take all of the princial components
allPCs = myPCA$scores

# Learn the GLM for P
#Combining both lists into one data.frame
#myComboData <- data.frame(myPCA$scores, myP)
P.fit = glm(myP~myPCA$scores)
summary(P.fit)


Q.fit = glm(myQ~myPCA$scores)
summary(Q.fit)

# Make the coefficients of the rest of the components 0
Qcoefs = coef(Q.fit)
Qcoefs[1:5]=0;
Qcoefs[7:9]=0;

# The weight (or influence) of each term on Q
wQ = V %*% Qcoefs[2:9]
wQ
