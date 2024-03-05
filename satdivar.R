### This code repository contains an R script used in Karbstein et al. 2020 (https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.6255) to calculate saturation of diversity variables per population/location (intraspecific functional trait diversity, genetic diversity based on microsatellite/SSR data and within-habitat heterogeneity (rarefaction-like analyses).

#set working directory where R script is located, for example
setwd("~/example/test/")


#### saturation of intraspecific functional trait diversity (iDFCV) ####

# create directory "genetic_diversity"
main_directory <- "~/example/test"

source_file <- "functional_traits.csv"
destination_directory <- "functional_trait_diversity"

# check if the destination directory exists, if not create it
if (!file.exists(destination_directory)) {
  dir.create(destination_directory)
}

# copy the file to the destination directory
file.copy(source_file, paste0(destination_directory, "/", source_file))


#set working directory where csv*files are located, for example

setwd("~/example/test/genetic_diversity")

if (!requireNamespace("adegenet", quietly = TRUE))
  install.packages("adegenet", repo="http://cran.rstudio.com/")
library(adegenet)

#### import data ####

dat_1<-read.csv("functional_traits.csv", header=TRUE, sep = ",", dec=".")

# check data
head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", repo="http://cran.rstudio.com/")
library(dplyr)


#### calculate functional trait diversity (iFDCV) between sample size 1-20 ####


# initialize results list
results <- list()

# create list of traits
traits <- unique(colnames(dat_1[, 3:11]))

for (k in 1:100) {  # replicates
  for (i in 1:20) {  # sample size range
    # Convert 'i' to a string with two digits
    i_str <- sprintf("%02d", i)
    
    # sample data
    subset_i <- dat_1 %>%
      group_by(location) %>%
      sample_n(size = i)
    
    new1 <- NULL  # initialize new1 within the loop
    
    for (trait in traits) {
      # calculate coefficient of variation per trait, sample size range, and iteration
      varK_trait <- tapply(subset_i[[trait]], subset_i$location, sd, na.rm = TRUE) / 
        tapply(subset_i[[trait]], subset_i$location, mean, na.rm = TRUE)
      # store results in new1
      new1[[paste0("varK_", trait)]] <- varK_trait
    }
    
    # append new1 to the results list
    results[[paste(i_str, k, "functional_traits_CV_diversity.csv", sep = "_")]] <- new1
    
    # write output
    write.csv(new1, file = paste(i_str, k, "functional_traits_CV_diversity.csv", sep = "_"))
  }
}


## create 20 folders and copy files into it
#give directory name here
directory <- "~/example/test/functional_trait_diversity/"


# check whether directory exists
if (!file.exists(directory)) {
  dir.create(directory)
}


# create folders in directory
for (i in 1:20) {
  folder_name <- sprintf("%02d", i)  # format the number with leading zeros
  dir.create(file.path(directory, folder_name), recursive = TRUE)
}

# create list of all subdirectories
target_subdirectories <- list.dirs(directory, full.names = TRUE, recursive = FALSE)


# copy files to respective folder
files <- list.files(directory, pattern = "^\\d+_.*\\.csv$", full.names = TRUE)


# iterate over data
for (file in files) {
  file_name <- basename(file)
  
  # extract prefix (e.g., "1_", "2_", etc.)
  prefix <- substr(file_name, 1, 2)
  
  # set target directory based on prefix
  target_subdirectory <- file.path(directory, prefix)
  
  # check whether target directory exists, or create one if not 
  if (!file.exists(target_subdirectory)) {
    dir.create(target_subdirectory)
  }
  # copy to target_subdirectory
  file.rename(file, file.path(target_subdirectory, file_name))
  # or:
  # file.copy(file, file.path(target_subdirectory, file_name))
}


# merge data per category (01:20)

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}


for(i in sprintf("%02d", 1:20)) {
  mymergeddata_i <- multmerge(paste("~/example/test/functional_trait_diversity/", i, sep = ""))
  
  res_i <- aggregate(mymergeddata_i[, 2], list(mymergeddata_i$X), mean)
  
  write.csv(res_i, file = paste("merged_", i, ".csv", sep = ""))
  
}


# create results folder and move intermediate files "merged_*.csv into it 

# create folders in directory
dir.create(file.path(directory, "results"), recursive = TRUE)

current_directories = list.dirs(directory, full.names = TRUE, recursive = FALSE)

results_directory <- grep("results", current_directories, value = TRUE)

files <- list.files(directory, pattern = "^merged_[0-9]{2}\\.csv$", full.names = TRUE)

# move files to the destination directory
for (file in files) {
  file.rename(file, file.path("results", basename(file)))
}



results = multmerge("~/example/test/functional_trait_diversity/results/")

# replicate size = 1 can only lead to CV=0 (NA) (in contrast, in GD, replicate size = 1 leads to positive HE values because genotypes can be heterozygous)

tapply(results$x, results$Group.1, length)

temp <- rep(c(2:20, 1), times = 13)

temp <- temp[-c(78,79,119)] # remove Eh=19,20, Gr=20

results$replicate <- temp

results$x[is.na(results$x)] <- 0

write.csv(results, file="intraspecific_trait_diversity_iFDCV.csv")



#### plot saturation curves ####

if (!requireNamespace("Rfast", quietly = TRUE))
  install.packages("Rfast", repo="http://cran.rstudio.com/")
library(Rfast)


### iFDCV ~ number of chosen samples (across all 13 populations)

#plot(results$x ~ results$replicate, type="p", pch=16, cex=1.0, cex.lab=1.5, cex.axis=1.5, xlab="no. of X", ylab="iFDCV", main="total", cex.main=1.5)


# create subets and plots 

# create list of group/pop_ids given in "results"
group_values <- unique(results$Group.1)

# create an empty list to store subsets
dat_list <- list()

# loop over each unique group value and create subsets
for (group_value in group_values) {
  # subset the results based on the current group value
  dat_list[[group_value]] <- subset(results, Group.1 == group_value)
}


# create an empty list to store the plots
plot_list <- list()

# set graphic parameters

par(mar=c(5.1,4.1,4.1,2.1) +  0.4)
par(mfrow=c(4,4))

# loop over each unique group value and create plots
for (group_value in group_values) {
  # extract the subset corresponding to the current group value
  subset_data <- dat_list[[group_value]]
  
  # create a plot for the current subset
  plot_list[[group_value]] <- plot(subset_data$replicate, subset_data$x, type="p", pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="iFDCV", main=group_value, cex.main=1.5)
  
  #### Assuming you have models fitted for each group value
  # Here's a placeholder for the predict function
  # Replace it with your actual model and prediction code
  #predicted_values <- predict(model, newdata = subset_data)
  
  # Draw lines based on predicted values
  #lines(subset_data$x, predicted_values, col = "red")
  
  # Store the plot in the plot_list
  #plot_list[[group_value]] <- recordPlot()
  
  # Close the current plot
  #dev.off()
}




##################

#### saturation of genetic diversity (GD) ####

#set working directory where R script is located, for example
setwd("~/example/test/")

# create directory "genetic_diversity"
main_directory <- "~/example/test/"
source_file <- "t_montanum_genepop.GEN"
destination_directory <- "genetic_diversity"

# check if the destination directory exists, if not create it
if (!file.exists(destination_directory)) {
  dir.create(destination_directory)
}

# copy the file to the destination directory
file.copy(source_file, paste0(destination_directory, "/", source_file))


# set working directory where csv*files are located, for example

setwd("~/example/test/genetic_diversity")

if (!requireNamespace("adegenet", quietly = TRUE))
  install.packages("adegenet", repo="http://cran.rstudio.com/")
library(adegenet)

#### import data ####
dat_1 <- import2genind("t_montanum_genepop.GEN", ncode = 2L, quiet = FALSE)


# check data structure
head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)


if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", repo="http://cran.rstudio.com/")
library(dplyr)


#### calculate genetic diversity between sample size 1-18 ####

new1<-NULL
subset_i<-NULL

for(k in 1:100){ # 100 replicates
  
  for(i in sprintf("%02d", 1:18)) { # sample size 1-18
    
    
    foo <- seppop(dat_1)# seperate pops
    foo
    
    mySamp <- lapply(foo, function(x) x[sample(1:nrow(x$tab), i)])
    mySamp
    
    
    x <- repool(mySamp)# put subsamples back to genind object
    
    new1<-Hs(genind2genpop(x))# calculate expected heterozygosity
    
    write.csv(new1, file=paste(i,k, "genetic_diversity.csv", sep="_"))# table called "18_100_genetic_diversity.csv"
    
    
  }
}

## create 18 folders and copy files therein
# give directory name here
directory <- "~/example/test/genetic_diversity/"


# check whether directory exists
if (!file.exists(directory)) {
  dir.create(directory)
}


# create folders in directory
for (i in 1:18) {
  folder_name <- sprintf("%02d", i)  # format the number with leading zeros
  dir.create(file.path(directory, folder_name), recursive = TRUE)
}

# create list of all subdirectories
target_subdirectories <- list.dirs(directory, full.names = TRUE, recursive = FALSE)


# copy files to respective folder
files <- list.files(directory, pattern = "^\\d+_.*\\.csv$", full.names = TRUE)


# iterate over data
for (file in files) {
  # extract prefix
  file_name <- basename(file)
  
  # extract prefix (e.g., "1_", "2_", etc.)
  prefix <- substr(file_name, 1, 2)
  
  # set target directory based on prefix
  target_subdirectory <- file.path(directory, prefix)
  
  # check whether target directory exists, or create one if not 
  if (!file.exists(target_subdirectory)) {
    dir.create(target_subdirectory)
  }
  # copy to target_subdirectory
  file.rename(file, file.path(target_subdirectory, file_name))
  # or:
  # file.copy(file, file.path(target_subdirectory, file_name))
}


# merge data per category (01:18)

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}


for(i in sprintf("%02d", 1:18)) {
  mymergeddata_i <- multmerge(paste("~/results/genetic_diversity/", i, sep = ""))
  
  res_i <- aggregate(mymergeddata_i[, 2], list(mymergeddata_i$X), mean)
  
  write.csv(res_i, file = paste("merged_", i, ".csv", sep = ""))
  
}


# create results folder and move intermediate files "merged_*.csv into it 

# create folders in directory
dir.create(file.path(directory, "results"), recursive = TRUE)

current_directories = list.dirs(directory, full.names = TRUE, recursive = FALSE)

results_directory <- grep("results", current_directories, value = TRUE)

files <- list.files(directory, pattern = "^merged_[0-9]{2}\\.csv$", full.names = TRUE)

# move files to the destination directory
for (file in files) {
  file.rename(file, file.path("results", basename(file)))
}



results = multmerge("~/example/test/genetic_diversity/results/")

results$replicate <- rep(1:18, times = 13)

write.csv(results, file="results_genetic_diversity_GD.csv")



#### plot saturation curves ####

if (!requireNamespace("Rfast", quietly = TRUE))
  install.packages("Rfast", repo="http://cran.rstudio.com/")
library(Rfast)


### GD ~ number of chosen samples (across all 13 populations)

#plot(results$x ~ results$replicate, type="p", pch=16, cex=1.0, cex.lab=1.5, cex.axis=1.5, xlab="no. of X", ylab="GD", main="total", cex.main=1.5)


# create subsets and plots 

# create list of group/pop_ids given in "results"
group_values <- unique(results$Group.1)

# create an empty list to store subsets
dat_list <- list()

# loop over each unique group value and create subsets
for (group_value in group_values) {
  # Subset the results based on the current group value
  dat_list[[group_value]] <- subset(results, Group.1 == group_value)
}


# create an empty list to store the plots
plot_list <- list()

# set graphic parameters

par(mar=c(5.1,4.1,4.1,2.1) +  0.4)
par(mfrow=c(4,4))

# loop over each unique group value and create plots
for (group_value in group_values) {
  # extract the subset corresponding to the current group value
  subset_data <- dat_list[[group_value]]

  # create a plot for the current subset
  plot_list[[group_value]] <- plot(subset_data$replicate, subset_data$x, type="p", pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main=group_value, cex.main=1.5)
  
  #### Assuming you have models fitted for each group value
  # Here's a placeholder for the predict function
  # Replace it with your actual model and prediction code
  #predicted_values <- predict(model, newdata = subset_data)
  
  # Draw lines based on predicted values
  #lines(subset_data$x, predicted_values, col = "red")
  
  # Store the plot in the plot_list
  #plot_list[[group_value]] <- recordPlot()
  
  # Close the current plot
  #dev.off()
}



##########

#### saturation of abiotic habitat heterogeneity (HD) ####
#set working directory where R script is located, for example
setwd("~/example/test/")

# create directory "genetic_diversity"
main_directory <- "~/example/test/"

source_file <- "environmental_parameters.csv"
destination_directory <- "habitat_heterogeneity"

# check if the destination directory exists, if not create it
if (!file.exists(destination_directory)) {
  dir.create(destination_directory)
}

# copy the file to the destination directory
file.copy(source_file, paste0(destination_directory, "/", source_file))


# set working directory where csv*files are located, for example

setwd("~/example/test/habitat_heterogeneity/")

if (!requireNamespace("adegenet", quietly = TRUE))
  install.packages("adegenet", repo="http://cran.rstudio.com/")
library(adegenet)

#### import data ####

dat_1<-read.csv("environmental_parameters.csv", header=TRUE, sep = ",", dec=".")

# check data
head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", repo="http://cran.rstudio.com/")
library(dplyr)


#### calculate abiotic habitat heterogeneity (HD) between sample size 1-20 ####

# Initialize results list
results <- list()

# create list of traits
factors <- unique(colnames(dat_1[, 4:13]))

for (k in 1:100) {  # replicates
  for (i in 1:4) {  # sample size range
    # convert 'i' to a string with two digits
    i_str <- sprintf("%02d", i)
    
    # sample data
    subset_i <- dat_1 %>%
      group_by(location) %>%
      sample_n(size = i)
    
    new1 <- NULL  # initialize new1 within the loop
    
    for (factor in factors) {
      # calculate coefficient of variation per factor, sample size range, and iteration
      varK_factor <- tapply(subset_i[[factor]], subset_i$location, sd, na.rm = TRUE) / 
        tapply(subset_i[[factor]], subset_i$location, mean, na.rm = TRUE)
      # store results in new1
      new1[[paste0("varK_", factor)]] <- varK_factor
    }
    
    # append new1 to the results list
    results[[paste(i_str, k, "environmental_parameters_CV_diversity.csv", sep = "_")]] <- new1
    
    # write output
    write.csv(new1, file = paste(i_str, k, "environmental_parameters_CV_diversity.csv", sep = "_"))
  }
}


## create 4 folders and copy files into it
# give directory name here
directory <- "~/example/test/habitat_heterogeneity/"

# check whether directory exists
if (!file.exists(directory)) {
  dir.create(directory)
}


# create folders in directory
for (i in 1:4) {
  folder_name <- sprintf("%02d", i)  # Format the number with leading zeros
  dir.create(file.path(directory, folder_name), recursive = TRUE)
}

# create list of all subdirectories
target_subdirectories <- list.dirs(directory, full.names = TRUE, recursive = FALSE)


# copy files to respective folder
files <- list.files(directory, pattern = "^\\d+_.*\\.csv$", full.names = TRUE)


# iterate over data
for (file in files) {
  file_name <- basename(file)
  
  # extract prefix (e.g., "1_", "2_", etc.)
  prefix <- substr(file_name, 1, 2)
  
  # set target directory based on prefix
  target_subdirectory <- file.path(directory, prefix)
  
  # check whether target directory exists, or create one if not 
  if (!file.exists(target_subdirectory)) {
    dir.create(target_subdirectory)
  }
  # copy to target_subdirectory
  file.rename(file, file.path(target_subdirectory, file_name))
  # or:
  # file.copy(file, file.path(target_subdirectory, file_name))
}


# merge data per category (01:05)

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}


for(i in sprintf("%02d", 1:5)) {
  mymergeddata_i <- multmerge(paste("~/results/habitat_heterogeneity/", i, sep = ""))
  
  res_i <- aggregate(mymergeddata_i[, 2], list(mymergeddata_i$X), mean)
  
  write.csv(res_i, file = paste("merged_", i, ".csv", sep = ""))
  
}


# create results folder and move intermediate files "merged_*.csv into it 

# create folders in directory
dir.create(file.path(directory, "results"), recursive = TRUE)

current_directories = list.dirs(directory, full.names = TRUE, recursive = FALSE)

results_directory <- grep("results", current_directories, value = TRUE)

files <- list.files(directory, pattern = "^merged_[0-9]{2}\\.csv$", full.names = TRUE)

# move files to the destination directory
for (file in files) {
  file.rename(file, file.path("results", basename(file)))
}



results = multmerge("~/example/test/habitat_heterogeneity/results/")

# replicate size = 1 can only lead to CV=0 (NA) (in contrast, in GD, replicate size = 1 leads to positive HE values because genotypes can be heterozygous)

#tapply(results$x, results$Group.1, length)

temp <- rep(c(2, 3, 4, 1), times = 13)

results$replicate <- temp

results$x[is.na(results$x)] <- 0

write.csv(results, file="environmental_parameter_diversity_HD.csv")



#### plot saturation curves ####

if (!requireNamespace("Rfast", quietly = TRUE))
  install.packages("Rfast", repo="http://cran.rstudio.com/")
library(Rfast)


### HD ~ number of chosen samples (across all 13 populations)

#plot(results$x ~ results$replicate, type="p", pch=16, cex=1.0, cex.lab=1.5, cex.axis=1.5, xlab="no. of X", ylab="iFDCV", main="total", cex.main=1.5)


# create subsets and plots 

# create list of group/pop_ids given in "results"
group_values <- unique(results$Group.1)

# create an empty list to store subsets
dat_list <- list()

# loop over each unique group value and create subsets
for (group_value in group_values) {
  # Subset the results based on the current group value
  dat_list[[group_value]] <- subset(results, Group.1 == group_value)
}


# create an empty list to store the plots
plot_list <- list()

# set graphic parameters

par(mar=c(5.1,4.1,4.1,2.1) +  0.4)
par(mfrow=c(4,4))

# loop over each unique group value and create plots
for (group_value in group_values) {
  # extract the subset corresponding to the current group value
  subset_data <- dat_list[[group_value]]
  
  # create a plot for the current subset
  plot_list[[group_value]] <- plot(subset_data$replicate, subset_data$x, type="p", pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main=group_value, cex.main=1.5)
  
  #### Assuming you have models fitted for each group value
  # Here's a placeholder for the predict function
  # Replace it with your actual model and prediction code
  #predicted_values <- predict(model, newdata = subset_data)
  
  # Draw lines based on predicted values
  #lines(subset_data$x, predicted_values, col = "red")
  
  # Store the plot in the plot_list
  #plot_list[[group_value]] <- recordPlot()
  
  # Close the current plot
  #dev.off()
}



