####FD (intraspecific functional trait diversity of populations) ####


##### data import#####

#set working directory where csv*files are located


dat_1<-read.csv("functional_traits.csv", header=TRUE, sep = ",", dec=".")

head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)

#install.packages("dplyr")
library(dplyr)


####total####

new1<-NULL
subset_i<-NULL

for(k in 1:100) { # number of replicates
  
  for(i in 1:20) { ##sample size
  
  subset_i<- dat_1 %>% group_by(dat_1$location) %>%   sample_n(size = i)
  
  
  ####calculate coefficients of variation of functional traits per population (location)
  
  
  # RH
  
  new1$varK_RH<- tapply(subset_i$RH, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$RH, subset_i$location, mean, na.rm=TRUE)
  
  
  
  # AGB
  
  new1$varK_AGB<-tapply(subset_i$AGB, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$AGB, subset_i$location, mean, na.rm=TRUE)
  
  
  
  
  # LA
  new1$varK_LA<-tapply(subset_i$LA, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$LA, subset_i$location, mean, na.rm=TRUE)
  
  
  
  # SLA
  new1$varK_SLA<-tapply(subset_i$SLA, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$SLA, subset_i$location, mean, na.rm=TRUE)
  
  
  
  # LDMC
  new1$varK_LDMC<-tapply(subset_i$LDMC, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$LDMC, subset_i$location, mean, na.rm=TRUE)
  
  
  
  
  # Fv/Fm
  new1$varK_FvFm<-tapply(subset_i$FvFm, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$FvFm, subset_i$location, mean, na.rm=TRUE)
  
  
  
  
  # PI
  new1$varK_PI<-tapply(subset_i$PI, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$PI, subset_i$location, mean, na.rm=TRUE)
  
  
  
  
  # SPS
  
  new1$varK_SPS<-tapply(subset_i$SPS, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$SPS, subset_i$location, mean, na.rm=TRUE)
  
  
  
  # SPI
  
  new1$varK_SPI<-tapply(subset_i$SPI, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$SPI, subset_i$location, mean, na.rm=TRUE)
  
  
  write.csv(new1, file=paste(i, k, "functional_traits_CV_diver_sat.csv", sep="_"))# write all combinations of iterations and sample size in the local directory
 
   
###merge all replicates per sample sizes (I recommend to create subfolders for each sample size)
    
}
}
#1
setwd("~/Desktop/test of diversity saturation/FD/1")

temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

mypath="~/Desktop/test of diversity saturation/FD/1"
multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

full_data = multmerge("~/Desktop/test of diversity saturation/FD/1/merged.csv")



mymergeddata_1 = multmerge("~/Desktop/test of diversity saturation/FD/1/")

res_1<-aggregate(mymergeddata_1[, 2:10], list(mymergeddata_1$X), mean)# ignore first column because its non-numeric


mymergeddata_2 = multmerge("~/Desktop/test of diversity saturation/FD/2/")
res_2<-aggregate(mymergeddata_2[, 2:10], list(mymergeddata_2$X), mean, na.rm=TRUE)

mymergeddata_3 = multmerge("~/Desktop/test of diversity saturation/FD/3/")
res_3<-aggregate(mymergeddata_3[, 2:10], list(mymergeddata_3$X), mean, na.rm=TRUE)


mymergeddata_4 = multmerge("~/Desktop/test of diversity saturation/FD/4/")
res_4<-aggregate(mymergeddata_4[, 2:10], list(mymergeddata_4$X), mean, na.rm=TRUE)

mymergeddata_5 = multmerge("~/Desktop/test of diversity saturation/FD/5/")
res_5<-aggregate(mymergeddata_5[, 2:10], list(mymergeddata_5$X), mean, na.rm=TRUE)

mymergeddata_6 = multmerge("~/Desktop/test of diversity saturation/FD/6/")
res_6<-aggregate(mymergeddata_6[, 2:10], list(mymergeddata_6$X), mean, na.rm=TRUE)

mymergeddata_7 = multmerge("~/Desktop/test of diversity saturation/FD/7/")
res_7<-aggregate(mymergeddata_7[, 2:10], list(mymergeddata_7$X), mean, na.rm=TRUE)

mymergeddata_8 = multmerge("~/Desktop/test of diversity saturation/FD/8/")
res_8<-aggregate(mymergeddata_8[, 2:10], list(mymergeddata_8$X), mean, na.rm=TRUE)

mymergeddata_9 = multmerge("~/Desktop/test of diversity saturation/FD/9/")
res_9<-aggregate(mymergeddata_9[, 2:10], list(mymergeddata_9$X), mean, na.rm=TRUE)

mymergeddata_10 = multmerge("~/Desktop/test of diversity saturation/FD/10/")
res_10<-aggregate(mymergeddata_10[, 2:10], list(mymergeddata_10$X), mean, na.rm=TRUE)

mymergeddata_11 = multmerge("~/Desktop/test of diversity saturation/FD/11/")
res_11<-aggregate(mymergeddata_11[, 2:10], list(mymergeddata_11$X), mean, na.rm=TRUE)

mymergeddata_12 = multmerge("~/Desktop/test of diversity saturation/FD/12/")
res_12<-aggregate(mymergeddata_12[, 2:10], list(mymergeddata_12$X), mean, na.rm=TRUE)

mymergeddata_13 = multmerge("~/Desktop/test of diversity saturation/FD/13/")
res_13<-aggregate(mymergeddata_13[, 2:10], list(mymergeddata_13$X), mean, na.rm=TRUE)


mymergeddata_14 = multmerge("~/Desktop/test of diversity saturation/FD/14/")
res_14<-aggregate(mymergeddata_14[, 2:10], list(mymergeddata_14$X), mean, na.rm=TRUE)

mymergeddata_15 = multmerge("~/Desktop/test of diversity saturation/FD/15/")
res_15<-aggregate(mymergeddata_15[, 2:10], list(mymergeddata_15$X), mean, na.rm=TRUE)

mymergeddata_16 = multmerge("~/Desktop/test of diversity saturation/FD/16/")
res_16<-aggregate(mymergeddata_16[, 2:10], list(mymergeddata_16$X), mean, na.rm=TRUE)

mymergeddata_17 = multmerge("~/Desktop/test of diversity saturation/FD/17/")
res_17<-aggregate(mymergeddata_17[, 2:10], list(mymergeddata_17$X), mean, na.rm=TRUE)

mymergeddata_18 = multmerge("~/Desktop/test of diversity saturation/FD/18/")
res_18<-aggregate(mymergeddata_18[, 2:10], list(mymergeddata_18$X), mean, na.rm=TRUE)

mymergeddata_19 = multmerge("~/Desktop/test of diversity saturation/FD/19/")
res_19<-aggregate(mymergeddata_19[, 2:10], list(mymergeddata_19$X), mean, na.rm=TRUE)

mymergeddata_20 = multmerge("~/Desktop/test of diversity saturation/FD/20/")
res_20<-aggregate(mymergeddata_20[, 2:10], list(mymergeddata_20$X), mean, na.rm=TRUE)






write.csv(res_1, file="merged_1.csv")
write.csv(res_2, file="merged_2.csv")
write.csv(res_3, file="merged_3.csv")
write.csv(res_4, file="merged_4.csv")
write.csv(res_5, file="merged_5.csv")
write.csv(res_6, file="merged_6.csv")
write.csv(res_7, file="merged_7.csv")
write.csv(res_8, file="merged_8.csv")
write.csv(res_9, file="merged_9.csv")
write.csv(res_10, file="merged_10.csv")
write.csv(res_11, file="merged_11.csv")
write.csv(res_12, file="merged_12.csv")
write.csv(res_13, file="merged_13.csv")
write.csv(res_14, file="merged_14.csv")
write.csv(res_15, file="merged_15.csv")
write.csv(res_16, file="merged_16.csv")
write.csv(res_17, file="merged_17.csv")
write.csv(res_18, file="merged_18.csv")
write.csv(res_19, file="merged_19.csv")
write.csv(res_20, file="merged_20.csv")


full_data_2 = multmerge("~/Desktop/test of diversity saturation/FD/res")
full_data_2$replicate<-c(2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1,2:20,1, 2:20,1)

write.csv(full_data_2, file="final.csv")



#install.packages("Rfast")
library(Rfast)

full_data_3<-read.csv(file="final.csv")

full_data_3$ifd_cv<-rowMeans(full_data_3[4:12], na.rm=TRUE)
#final$mean <- rowMeans(final, na.rm=TRUE)



plot(full_data_3$ifd_cv ~ full_data_3$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of X", ylab=main=expression("IFDCV"[2]), main="KW", cex.main=1.5)



#subets

dat_Ba<-subset(full_data_3, full_data_3$Group.1=="Ba")
dat_Bo<-subset(full_data_3, full_data_3$Group.1=="Bo")
dat_Di<-subset(full_data_3, full_data_3$Group.1=="Di")
dat_Eh<-subset(full_data_3, full_data_3$Group.1=="Eh")
dat_Er<-subset(full_data_3, full_data_3$Group.1=="Er")
dat_Gr<-subset(full_data_3, full_data_3$Group.1=="Gr")
dat_Ha<-subset(full_data_3, full_data_3$Group.1=="Ha")
dat_If<-subset(full_data_3, full_data_3$Group.1=="If")
dat_KW<-subset(full_data_3, full_data_3$Group.1=="KW")
dat_Ni<-subset(full_data_3, full_data_3$Group.1=="Ni")
dat_Sa<-subset(full_data_3, full_data_3$Group.1=="Sa")
dat_St<-subset(full_data_3, full_data_3$Group.1=="St")
dat_Wo<-subset(full_data_3, full_data_3$Group.1=="Wo")


# plot
par(mar=c(5.1,4.1,4.1,2.1) +  0.4)

par(mfrow=c(3,5))

plot(dat_KW$ifd_cv ~ dat_KW$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="KW", cex.main=1.5)

plot(dat_Bo$ifd_cv ~ dat_Bo$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Bo", cex.main=1.5)





plot(dat_Ha$ifd_cv ~ dat_Ha$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Ha", cex.main=1.5)



plot(dat_Wo$ifd_cv ~ dat_Wo$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Wo", cex.main=1.5)


plot(dat_Ba$ifd_cv ~ dat_Ba$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Ba", cex.main=1.5)


plot(dat_St$ifd_cv ~ dat_St$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="St", cex.main=1.5)


plot(dat_Sa$ifd_cv ~ dat_Sa$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Sa", cex.main=1.5)


plot(dat_If$ifd_cv ~ dat_If$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="If", cex.main=1.5)




plot(dat_Ni$ifd_cv ~ dat_Ni$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Ni", cex.main=1.5)



plot(dat_Di$ifd_cv ~ dat_Di$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Di", cex.main=1.5)





plot(dat_Er$ifd_cv ~ dat_Er$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Er", cex.main=1.5)



plot(dat_Gr$ifd_cv ~ dat_Gr$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Gr", cex.main=1.5)


plot(dat_Eh$ifd_cv ~ dat_Eh$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab=expression("iFD"[CV]), main="Eh", cex.main=1.5)










#### HD (within-habitat heterogeneity of locations)####


##### data import#####

#set working directory where csv*files are located

setwd("~/Desktop/test of diversity saturation/HD")

dat_1<-read.csv("environmental_parameters.csv", header=TRUE, sep = ",", dec=".")

head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)

#install.packages("dplyr")
library(dplyr)


dat_1<-read.csv("environmental_parameters.csv", header=TRUE, sep = ",", dec=".")

head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)


####total####

new1<-NULL
subset_i<-NULL

for(k in 1:100){ #number of replicates
  
  for(i in 1:5) { #sample size per location
    
    subset_i<- dat_1 %>% group_by(dat_1$location) %>%   sample_n(size = i)
    
    
    
    
    
    # altitude
    
    new1$varK_RH<- tapply(subset_i$altitude, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$altitude, subset_i$location, mean, na.rm=TRUE)
    
    
    
    # slope
    
    new1$varK_slope<-tapply(subset_i$slope, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$slope, subset_i$location, mean, na.rm=TRUE)
    
    
    
    
    # soil_depth
    new1$varK_soil_depth<-tapply(subset_i$soil_depth, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$soil_depth, subset_i$location, mean, na.rm=TRUE)
    
    
    
    # lai
    new1$varK_lai<-tapply(subset_i$lai, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$lai, subset_i$location, mean, na.rm=TRUE)
    
    
    
    # exposition
    new1$varK_exposition<-tapply(subset_i$exposition, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$exposition, subset_i$location, mean, na.rm=TRUE)
    
    
    
    
    # CECpot
    new1$varK_CECpot<-tapply(subset_i$CECpot, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$CECpot, subset_i$location, mean, na.rm=TRUE)
    
    
    
    
    # pHH20
    new1$varK_pHH20<-tapply(subset_i$pHH20, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$pHH20, subset_i$location, mean, na.rm=TRUE)
    
    
    
    
    # N
    
    new1$varK_N<-tapply(subset_i$N, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$N, subset_i$location, mean, na.rm=TRUE)
    
    
    
    # P
    
    new1$varK_P<-tapply(subset_i$P, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$P, subset_i$location, mean, na.rm=TRUE)
    
   
    new1$varK_K<-tapply(subset_i$K, subset_i$location, sd, na.rm=TRUE)/tapply(subset_i$K, subset_i$location, mean, na.rm=TRUE)
    
     
    write.csv(new1, file=paste(i, k, "habitat_heterogeneity_CV_diver_sat.csv", sep="_"))
    
###merge all replicates per sample sizes (I recommend to create subfolders for each sample size)
    
  }
}
#1
#setwd("~/Desktop/test of diversity saturation/HD/1")

#temp = list.files(pattern="*.csv")
#for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

#mypath="~/Desktop/test of diversity saturation/FD/1"

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

#full_data = multmerge("~/Desktop/test of diversity saturation/FD/1/merged.csv")



mymergeddata_1 = multmerge("~/Desktop/test of diversity saturation/HD/1/")
res_1<-aggregate(mymergeddata_1[, 2:11], list(mymergeddata_1$X), mean)


mymergeddata_2 = multmerge("~/Desktop/test of diversity saturation/HD/2/")
res_2<-aggregate(mymergeddata_2[, 2:11], list(mymergeddata_2$X), mean, na.rm=TRUE)

mymergeddata_3 = multmerge("~/Desktop/test of diversity saturation/HD/3/")
res_3<-aggregate(mymergeddata_3[, 2:11], list(mymergeddata_3$X), mean, na.rm=TRUE)


mymergeddata_4 = multmerge("~/Desktop/test of diversity saturation/HD/4/")
res_4<-aggregate(mymergeddata_4[, 2:11], list(mymergeddata_4$X), mean, na.rm=TRUE)

mymergeddata_5 = multmerge("~/Desktop/test of diversity saturation/HD/5/")
res_5<-aggregate(mymergeddata_5[, 2:11], list(mymergeddata_5$X), mean, na.rm=TRUE)

write.csv(res_1, file="merged_1.csv")
write.csv(res_2, file="merged_2.csv")
write.csv(res_3, file="merged_3.csv")
write.csv(res_4, file="merged_4.csv")
write.csv(res_5, file="merged_5.csv")


full_data_2 = multmerge("~/Desktop/test of diversity saturation/HD/res")
full_data_2$replicate<-c(2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1,2:5,1, 2:5,1)

write.csv(full_data_2, file="final.csv")



#install.packages("Rfast")
library(Rfast)

full_data_3<-read.csv(file="final.csv")

full_data_3$ifd_cv<-rowMeans(full_data_3[4:13], na.rm=TRUE)


#final$mean <- rowMeans(final, na.rm=TRUE)



plot(full_data_3$ifd_cv ~ full_data_3$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of X", ylab="HD", main="KW", cex.main=1.5)



#subets

dat_Ba<-subset(full_data_3, full_data_3$Group.1=="Ba")
dat_Bo<-subset(full_data_3, full_data_3$Group.1=="Bo")
dat_Di<-subset(full_data_3, full_data_3$Group.1=="Di")
dat_Eh<-subset(full_data_3, full_data_3$Group.1=="Eh")
dat_Er<-subset(full_data_3, full_data_3$Group.1=="Er")
dat_Gr<-subset(full_data_3, full_data_3$Group.1=="Gr")
dat_Ha<-subset(full_data_3, full_data_3$Group.1=="Ha")
dat_If<-subset(full_data_3, full_data_3$Group.1=="If")
dat_KW<-subset(full_data_3, full_data_3$Group.1=="KW")
dat_Ni<-subset(full_data_3, full_data_3$Group.1=="Ni")
dat_Sa<-subset(full_data_3, full_data_3$Group.1=="Sa")
dat_St<-subset(full_data_3, full_data_3$Group.1=="St")
dat_Wo<-subset(full_data_3, full_data_3$Group.1=="Wo")


# plot

par(mar=c(5.1,4.1,4.1,2.1) +  0.4)
par(mfrow=c(3,5))

plot(dat_KW$ifd_cv ~ dat_KW$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="KW", cex.main=1.5)

plot(dat_Bo$ifd_cv ~ dat_Bo$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Bo", cex.main=1.5)





plot(dat_Ha$ifd_cv ~ dat_Ha$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Ha", cex.main=1.5)



plot(dat_Wo$ifd_cv ~ dat_Wo$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Wo", cex.main=1.5)


plot(dat_Ba$ifd_cv ~ dat_Ba$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Ba", cex.main=1.5)


plot(dat_St$ifd_cv ~ dat_St$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="St", cex.main=1.5)


plot(dat_Sa$ifd_cv ~ dat_Sa$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Sa", cex.main=1.5)


plot(dat_If$ifd_cv ~ dat_If$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="If", cex.main=1.5)




plot(dat_Ni$ifd_cv ~ dat_Ni$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Ni", cex.main=1.5)



plot(dat_Di$ifd_cv ~ dat_Di$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Di", cex.main=1.5)





plot(dat_Er$ifd_cv ~ dat_Er$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Er", cex.main=1.5)



plot(dat_Gr$ifd_cv ~ dat_Gr$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Gr", cex.main=1.5)


plot(dat_Eh$ifd_cv ~ dat_Eh$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="HD", main="Eh", cex.main=1.5)




####GD (genetic diversity (or expected heterozygosity) based on SSR data ####


##### data import#####

#set working directory where csv*files are located

#install.packages("adegenet")
library(adegenet)

setwd("~/Desktop/test of diversity saturation/GD")

dat_1 <- import2genind("T_montanum_genepop.GEN", ncode = 2L, quiet = FALSE)



head(dat_1)
tail(dat_1)

str(dat_1)
summary(dat_1)



#install.packages("dplyr")
library(dplyr)


####total####



new1<-NULL
subset_i<-NULL

for(k in 1:100){ #number of replicates
  
  for(i in 1:18) { #sample size
    
    
    foo <- seppop(dat_1)#seperate pops
    foo
    
    mySamp <- lapply(foo, function(x) x[sample(1:nrow(x$tab), i)])
    mySamp
    
    
    x <- repool(mySamp)# put subsamples back to genind object
    
    new1<-Hs(genind2genpop(x))
    
    write.csv(new1, file=paste(i,k, "genetic_diversity.csv", sep="_"))
   
    
     ###merge all replicates per sample sizes (I recommend to create subfolders for each sample        size)
    
    
  }
}


multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

#full_data = multmerge("~/Desktop/test of diversity saturation/FD/1/merged.csv")



mymergeddata_1 = multmerge("~/Desktop/test of diversity saturation/GD/1/")

res_1<-aggregate(mymergeddata_1[, 2], list(mymergeddata_1$X), mean)


mymergeddata_2 = multmerge("~/Desktop/test of diversity saturation/GD/2/")
res_2<-aggregate(mymergeddata_2[, 2], list(mymergeddata_2$X), mean, na.rm=TRUE)

mymergeddata_3 = multmerge("~/Desktop/test of diversity saturation/GD/3/")
res_3<-aggregate(mymergeddata_3[, 2], list(mymergeddata_3$X), mean, na.rm=TRUE)


mymergeddata_4 = multmerge("~/Desktop/test of diversity saturation/GD/4/")
res_4<-aggregate(mymergeddata_4[, 2], list(mymergeddata_4$X), mean, na.rm=TRUE)

mymergeddata_5 = multmerge("~/Desktop/test of diversity saturation/GD/5/")
res_5<-aggregate(mymergeddata_5[, 2], list(mymergeddata_5$X), mean, na.rm=TRUE)

mymergeddata_6 = multmerge("~/Desktop/test of diversity saturation/GD/6/")
res_6<-aggregate(mymergeddata_6[, 2], list(mymergeddata_6$X), mean, na.rm=TRUE)

mymergeddata_7 = multmerge("~/Desktop/test of diversity saturation/GD/7/")
res_7<-aggregate(mymergeddata_7[, 2], list(mymergeddata_7$X), mean, na.rm=TRUE)

mymergeddata_8 = multmerge("~/Desktop/test of diversity saturation/GD/8/")
res_8<-aggregate(mymergeddata_8[, 2], list(mymergeddata_8$X), mean, na.rm=TRUE)

mymergeddata_9 = multmerge("~/Desktop/test of diversity saturation/GD/9/")
res_9<-aggregate(mymergeddata_9[, 2], list(mymergeddata_9$X), mean, na.rm=TRUE)

mymergeddata_10 = multmerge("~/Desktop/test of diversity saturation/GD/10/")
res_10<-aggregate(mymergeddata_10[, 2], list(mymergeddata_10$X), mean, na.rm=TRUE)

mymergeddata_11 = multmerge("~/Desktop/test of diversity saturation/GD/11/")
res_11<-aggregate(mymergeddata_11[, 2], list(mymergeddata_11$X), mean, na.rm=TRUE)

mymergeddata_12 = multmerge("~/Desktop/test of diversity saturation/GD/12/")
res_12<-aggregate(mymergeddata_12[, 2], list(mymergeddata_12$X), mean, na.rm=TRUE)

mymergeddata_13 = multmerge("~/Desktop/test of diversity saturation/GD/13/")
res_13<-aggregate(mymergeddata_13[, 2], list(mymergeddata_13$X), mean, na.rm=TRUE)


mymergeddata_14 = multmerge("~/Desktop/test of diversity saturation/GD/14/")
res_14<-aggregate(mymergeddata_14[, 2], list(mymergeddata_14$X), mean, na.rm=TRUE)

mymergeddata_15 = multmerge("~/Desktop/test of diversity saturation/GD/15/")
res_15<-aggregate(mymergeddata_15[, 2], list(mymergeddata_15$X), mean, na.rm=TRUE)

mymergeddata_16 = multmerge("~/Desktop/test of diversity saturation/GD/16/")
res_16<-aggregate(mymergeddata_16[, 2], list(mymergeddata_16$X), mean, na.rm=TRUE)

mymergeddata_17 = multmerge("~/Desktop/test of diversity saturation/GD/17/")
res_17<-aggregate(mymergeddata_17[, 2], list(mymergeddata_17$X), mean, na.rm=TRUE)

mymergeddata_18 = multmerge("~/Desktop/test of diversity saturation/GD/18/")
res_18<-aggregate(mymergeddata_18[, 2], list(mymergeddata_18$X), mean, na.rm=TRUE)


write.csv(res_1, file="merged_1.csv")
write.csv(res_2, file="merged_2.csv")
write.csv(res_3, file="merged_3.csv")
write.csv(res_4, file="merged_4.csv")
write.csv(res_5, file="merged_5.csv")
write.csv(res_6, file="merged_6.csv")
write.csv(res_7, file="merged_7.csv")
write.csv(res_8, file="merged_8.csv")
write.csv(res_9, file="merged_9.csv")
write.csv(res_10, file="merged_10.csv")
write.csv(res_11, file="merged_11.csv")
write.csv(res_12, file="merged_12.csv")
write.csv(res_13, file="merged_13.csv")
write.csv(res_14, file="merged_14.csv")
write.csv(res_15, file="merged_15.csv")
write.csv(res_16, file="merged_16.csv")
write.csv(res_17, file="merged_17.csv")
write.csv(res_18, file="merged_18.csv")


full_data_2 = multmerge("~/Desktop/test of diversity saturation/GD/res")
full_data_2$replicate<-c(1:18,1:18,1:18,1:18,1:18,1:18,1:18,1:18,1:18,1:18,1:18,1:18,1:18)

write.csv(full_data_2, file="final.csv")



#install.packages("Rfast")
library(Rfast)

full_data_3<-read.csv(file="final.csv")



plot(full_data_3$x ~ full_data_3$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of X", ylab="GD", main="total", cex.main=1.5)



#subets

dat_Ba<-subset(full_data_3, full_data_3$Group.1=="Ba_19")
dat_Bo<-subset(full_data_3, full_data_3$Group.1=="Bo_20")
dat_Di<-subset(full_data_3, full_data_3$Group.1=="Di_20")
dat_Eh<-subset(full_data_3, full_data_3$Group.1=="Eh_20")
dat_Er<-subset(full_data_3, full_data_3$Group.1=="Er_20")
dat_Gr<-subset(full_data_3, full_data_3$Group.1=="Gr_19")
dat_Ha<-subset(full_data_3, full_data_3$Group.1=="Ha_20")
dat_If<-subset(full_data_3, full_data_3$Group.1=="If_18")
dat_KW<-subset(full_data_3, full_data_3$Group.1=="KW_20")
dat_Ni<-subset(full_data_3, full_data_3$Group.1=="Ni_20")
dat_Sa<-subset(full_data_3, full_data_3$Group.1=="Sa_20")
dat_St<-subset(full_data_3, full_data_3$Group.1=="St_19")
dat_Wo<-subset(full_data_3, full_data_3$Group.1=="Wo_20")


# plot
par(mar=c(5.1,4.1,4.1,2.1) +  0.4)
par(mfrow=c(3,5))

plot(dat_KW$x ~ dat_KW$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="KW", cex.main=1.5)

plot(dat_Bo$x ~ dat_Bo$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Bo", cex.main=1.5)





plot(dat_Ha$x ~ dat_Ha$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Ha", cex.main=1.5)



plot(dat_Wo$x ~ dat_Wo$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Wo", cex.main=1.5)


plot(dat_Ba$x ~ dat_Ba$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Ba", cex.main=1.5)


plot(dat_St$x ~ dat_St$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="St", cex.main=1.5)


plot(dat_Sa$x ~ dat_Sa$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Sa", cex.main=1.5)


plot(dat_If$x ~ dat_If$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="If", cex.main=1.5)




plot(dat_Ni$x ~ dat_Ni$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Ni", cex.main=1.5)



plot(dat_Di$x ~ dat_Di$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Di", cex.main=1.5)





plot(dat_Er$x ~ dat_Er$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Er", cex.main=1.5)



plot(dat_Gr$x ~ dat_Gr$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Gr", cex.main=1.5)


plot(dat_Eh$x ~ dat_Eh$replicate, type="p", pch=16, cex=2.5, cex.lab=1.5, cex.axis=1.5, xlab="no. of samples", ylab="GD", main="Eh", cex.main=1.5)








#### I hope that this exemplary R script was useful for you to calculate saturation of diversity variables. Please do not hesitat to contact me!####











