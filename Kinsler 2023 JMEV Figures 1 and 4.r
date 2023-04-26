
#Load data and packages
##############################################################
library(reshape2)
library(ggplot2)

#First load the fitness data from 1 Big Batch:
read.table(paste("/users/ksamerot/Dropbox\ (ASU)/Noise\ Paper/Final/Kinsler_et_al_2020_fitnessdata.csv"), sep = ",", header = TRUE)->Kinsler

# Pull out fitnesses and reshape table
library(reshape2)
M3F<-Kinsler[,grep("fitness", colnames(Kinsler))]
M3F$barcode<-Kinsler$barcode
M3F$gene<-Kinsler$mutation_type
M<-melt(M3F, id.var = c("barcode", "gene"))

# Add the error column back
M3F<-Kinsler[,grep("error", colnames(Kinsler))]
M3F$barcode<-Kinsler$barcode
M3F$gene<-Kinsler$mutation_type
Merror<-melt(M3F, id.var = c("barcode", "gene"))
M$error<-Merror$value
head(M)

#Pull out appropriate columns that refer to the M3 condition
Replicates<-unique(M$variable)[c(11:16, 19:33, 35, 37, 43, 45, 92,93,94)]
M<-subset(M, M$variable %in% Replicates)
head(M)

#Assign identifying columns back to this table:

#Batch
M$Batch<-as.character(lapply(M$variable, function(x) strsplit(as.character(x), split = "\\.")[[1]][1]))
M$Batch[which(M$Batch == "A_fitness")]<-"BB"
M$Batch[which(M$Batch == "B_fitness")]<-"BB"
M$Batch[which(M$Batch == "C_fitness")]<-"BB"
M$Batch[which(M$Batch == "D_fitness")]<-"BB"
M$Batch[which(M$Batch == "t19")]<-"X19"
unique(M$Batch)

#Replicate
M$Rep<-as.character(lapply(M$variable, function(x) strsplit(as.character(x), split = "\\.")[[1]][2]))
M$Rep<-as.character(lapply(M$Rep, function(x) strsplit(as.character(x), split = "_")[[1]][1]))
M$Rep[which(M$Batch == "BB")]<-as.character(lapply(M$variable[which(M$Batch == "BB")], function(x) strsplit(as.character(x), split = "_")[[1]][1]))
M$Rep[which(M$Rep == "1")]<-"A"
M$Rep[which(M$Rep == "2")]<-"B"
M$Rep[which(M$Rep == "3")]<-"C"
head(M)
unique(M$Rep)

M$BioRep<-paste(M$Batch, M$Rep, sep = "_")

#Date
M$Date<-"none"
M$Date[which(M$Batch=="BB")]<-"12/10/17"
M$Date[which(M$Batch=="X19")]<-"5/9/15"
M$Date[which(M$Batch=="X13")]<-"5/1/15"
M$Date[which(M$Batch=="X21")]<-"8/17/15"
M$Date[which(M$Batch=="X18")]<-"5/9/15"
M$Date[which(M$Batch=="X20")]<-"8/17/15"
M$Date[which(M$Batch=="X23")]<-"9/4/15"
M$Date[which(M$Batch=="X6")]<-"12/26/14"
M$Date[which(M$Batch=="X3")]<-"11/24/14"

################################
#The above table is missing the random samples, so I need to add those back
#We actually omitted these random samples from the paper because they did not add much
#We were trying to show that random sampling results in lower noise than the amount of noise we found, which is already onvious
read.csv("/users/ksamerot/Dropbox\ (ASU)/Noise Paper/Final/Fitness_Inferred_071221.csv")->Random_Samples

#eliminate columsn with the word "error"
library(reshape2)
MRandom<-Random_Samples[,grep("fitness", colnames(Random_Samples))]
MRandom$barcode<-Random_Samples$barcode
MRandom<-melt(MRandom, id.var = "barcode")

#add the error column back
temp<-Random_Samples[,grep("error", colnames(Random_Samples))]
temp$barcode<-Random_Samples$barcode
Merror<-melt(temp, id.var = "barcode")
MRandom$error<-Merror$value
head(MRandom)

#BackupMRandom because I realied I still need it
MSpike<-MRandom

#All we need from this is the random samples
MRandom<-MRandom[grep("Random", MRandom$variable),]

#Assign identifying columns back to this table:
MRandom$Batch<-"BB"
MRandom$Rep<-as.character(lapply(MRandom$variable, function(x) strsplit(as.character(x), split = "\\.")[[1]][2]))
MRandom$BioRep<-paste(MRandom$Batch, MRandom$Rep, sep = "_")
MRandom$gene<-M$gene[match(MRandom$barcode, M$barcode)]
MRandom$Date<-"12/10/17"

#Rbind the randoms back to the full table
MRandom<-MRandom[,c(1,8,2,3,4,5,6,7,9)]
MRandom$Type<-"Random"
head(MRandom)
M$Type<-"Original"
M<-rbind(M, MRandom)
head(M)


head(M)

###########################################
#The above table is also missing the spike ins, so lets add those back too:
MSpike<-MSpike[which(MSpike$barcode > 900000),]
MSpike<-MSpike[grep("Original", MSpike$variable),]

#Assign identifying columns back to this table:
MSpike$Batch<-"BB"
MSpike$Rep<-as.character(lapply(MSpike$variable, function(x) strsplit(as.character(x), split = "\\.")[[1]][2]))
MSpike$BioRep<-paste(MSpike$Batch, MSpike$Rep, sep = "_")
MSpike$gene<-"IRA1mis*"
MSpike$gene[which(MSpike$barcode>900099)]<-"IRA1non*"
MSpike$Date<-"12/10/17"
MSpike<-MSpike[,c(1,8,2,3,4,5,6,7,9)]
MSpike$Type<-"Original"
head(MSpike)
M<-rbind(M, MSpike)
head(M)

#Eliminate barcodes that are abberant meaning they have unique fitness values in one of the conditions studied in Kinsler et al 2020 and likely picked up a mutation.
M<-M[-grep("900010", M$barcode),]
M<-M[-grep("134852", M$barcode),]
M<-M[-grep("900110", M$barcode),]
M<-M[-grep("9000110", M$barcode),]
M<-M[-grep("900109", M$barcode),]
M<-M[-grep("900009", M$barcode),]

#Uniform gene names
M$gene[which(M$gene == "IRA1_nonsense")]<-"IRA1non"
M$gene[which(M$gene == "IRA1_missense")]<-"IRA1mis"
M$gene[which(M$gene == "ExpNeutral")]<-"Unmutated"

M_Backup<-M


##############################
#Do the calculations required for figure 2A and 2B

##################################
#Noise between random samples
###############################################
#Only look at 1BB since I only took random samples from there
#Total of 5 random samples + original = 6
#13 GPB2 mutants get 6 samples aech
#9Ira1mis get 6 samples each
#7IRA1non get 6 samples each

#This is for calculating stedev without using the error from GrantFitSeq
#This random sampling was eventually removed from the manuscript nased on advice from one reviewer.
WI_Flask<-subset(M, M$Type == "Random" & (M$gene == "IRA1non" | M$gene ==  "IRA1mis"| M$gene ==  "GPB2"))
WI_Flask<-na.omit(WI_Flask)
WI_Flask<-subset(WI_Flask, WI_Flask$value != "Inf")
head(WI_Flask)

WId<-aggregate(value~gene+Batch+barcode+Rep, WI_Flask, FUN = length)
WId
WIm<-aggregate(value~gene+Batch+barcode+Rep, WI_Flask, FUN = mean)
WIm
WIs<-aggregate(value~gene+Batch+barcode+Rep, WI_Flask, FUN = sd)
WIcv<-WIm
WIcv$value<-abs(WIs$value/WIm$value)
WIcv$mean<-WIm$value
WIcv$sd<-WIs$value
head(WIcv)
WIcv$n<-6
WIcv$CI<-1.96* ((WIcv$sd)/sqrt(WIcv$n))
WIcv$Percent_Error<-(WIcv$CI/WIcv$mean)*100
SIMcv<-WIcv
head(SIMcv)

#This is for calculating $error using the error from GrantFitSeq

WI_Flask$Weighted_Means<-WI_Flask$value/WI_Flask$error
WI_Flask$Inverse_Variances<-1/WI_Flask$error
Sum_of_weighted_means<-aggregate(Weighted_Means~gene+Batch+barcode+Rep, WI_Flask, FUN = sum)
Sum_of_inverse_variances<-aggregate(Inverse_Variances~gene+Batch+barcode+Rep, WI_Flask, FUN = sum)
Sum_of_weighted_means$Inverse_variances<-Sum_of_inverse_variances$Inverse_Variances
Sum_of_weighted_means$IWM<-Sum_of_weighted_means$Weighted_Means/Sum_of_weighted_means$Inverse_variances
Sum_of_weighted_means$IWV<-1/Sum_of_weighted_means$Inverse_variances
Sum_of_weighted_means$CI<-(1.96*(Sum_of_weighted_means$IWV))
Sum_of_weighted_means$Percent_Error<-(Sum_of_weighted_means$CI/Sum_of_weighted_means$IWM)*100
head(Sum_of_weighted_means)
Sum_of_weighted_means_A<-Sum_of_weighted_means

################################
#Noise between spike ins (or within flask) samples:
#Find the variance for the 4 categories that have it within flasks:
#18 flasks were studied for GPB2
#18 flasks were studied for NSI
#4 flasks only for IRA1misspike
#4flasks only for I1NS
WI_Flask<-subset(M, M$Type == "Original" & (M$gene == "IRA1non*" | M$gene ==  "IRA1mis*"| M$gene ==  "GPB2"))
WI_Flask<-na.omit(WI_Flask)
WI_Flask<-subset(WI_Flask, WI_Flask$value != "Inf")
#WI_Flask<-WI_Flask[-grep("SV", WI_Flask$variable),]

head(WI_Flask)

WIm<-aggregate(value~variable+gene, WI_Flask, FUN = mean)
WIm
WIs<-aggregate(value~variable+gene, WI_Flask, FUN = sd)
WIcv<-WIm
WIcv$value<-abs(WIs$value/WIm$value)
WIcv$mean<-WIm$value
WIcv$sd<-WIs$value
head(WIcv)

WIcv$n<-0
WIcv$n[which(WIcv$gene == "GPB2")]<-length(unique(subset(WI_Flask, WI_Flask$gene == "GPB2")$barcode))
WIcv$n[which(WIcv$gene == "IRA1non*")]<-length(unique(subset(WI_Flask, WI_Flask$gene == "IRA1non*")$barcode))
WIcv$n[which(WIcv$gene == "IRA1mis*")]<-length(unique(subset(WI_Flask, WI_Flask$gene == "IRA1mis*")$barcode))
WIcv$CI<-1.96* ((WIcv$sd)/sqrt(WIcv$n))
WIcv$Percent_Error<-(WIcv$CI/WIcv$mean)*100
head(WIcv)


#Calculating 95% CIs for spike in (within flask) samples:
WI_Flask$Weighted_Means<-WI_Flask$value/WI_Flask$error
WI_Flask$Inverse_Variances<-1/WI_Flask$error
Sum_of_weighted_means<-aggregate(Weighted_Means~variable+gene, WI_Flask, FUN = sum)
Sum_of_inverse_variances<-aggregate(Inverse_Variances~variable+gene, WI_Flask, FUN = sum)
Sum_of_weighted_means$Inverse_variances<-Sum_of_inverse_variances$Inverse_Variances
Sum_of_weighted_means$IWM<-Sum_of_weighted_means$Weighted_Means/Sum_of_weighted_means$Inverse_variances
Sum_of_weighted_means$IWV<-1/Sum_of_weighted_means$Inverse_variances
Sum_of_weighted_means$CI<-(1.96*(Sum_of_weighted_means$IWV))
Sum_of_weighted_means$Percent_Error<-(Sum_of_weighted_means$CI/Sum_of_weighted_means$IWM)*100
head(Sum_of_weighted_means)
Sum_of_weighted_means_B<-Sum_of_weighted_means

################################
#Noise between replicates
#That was for within each flask, now lets do across replicates

WI_Flask<-subset(M,  M$Type == "Original" & (M$gene ==  "IRA1non" | M$gene ==  "IRA1mis"| M$gene ==  "GPB2"))
WI_Flask<-na.omit(WI_Flask)
WI_Flask<-subset(WI_Flask, WI_Flask$value != "Inf")
Bm<-aggregate(value~Batch+barcode+gene, WI_Flask, FUN = mean)
Bs<-aggregate(value~Batch+barcode+gene, WI_Flask, FUN = sd)
Bcv<-Bm
Bcv$value<-abs(Bs$value/Bm$value)
Bcv$mean<-Bm$value
Bcv$sd<-Bs$value
Bcv$n<-3
Bcv$n[which(Bcv$Batch == "X1Big")]<-4
Bcv$CI<-1.96* ((Bcv$sd)/sqrt(Bcv$n))
Bcv$Percent_Error<-(Bcv$CI/Bcv$mean)*100
head(Bcv)

#Calculating 95% CIs for replicates samples:
WI_Flask$Weighted_Means<-WI_Flask$value/WI_Flask$error
WI_Flask$Inverse_Variances<-1/WI_Flask$error
Sum_of_weighted_means<-aggregate(Weighted_Means~Batch+barcode+gene, WI_Flask, FUN = sum)
Sum_of_inverse_variances<-aggregate(Inverse_Variances~Batch+barcode+gene, WI_Flask, FUN = sum)
Sum_of_weighted_means$Inverse_variances<-Sum_of_inverse_variances$Inverse_Variances
Sum_of_weighted_means$IWM<-Sum_of_weighted_means$Weighted_Means/Sum_of_weighted_means$Inverse_variances
Sum_of_weighted_means$IWV<-1/Sum_of_weighted_means$Inverse_variances
Sum_of_weighted_means$CI<-(1.96*(Sum_of_weighted_means$IWV))
Sum_of_weighted_means$Percent_Error<-(Sum_of_weighted_means$CI/Sum_of_weighted_means$IWM)*100
head(Sum_of_weighted_means)
Sum_of_weighted_means_C<-Sum_of_weighted_means


################################
#Noise between batches
##################################Now how about between batches
WI_Flask<-subset(M, M$Type == "Original" & (M$gene ==  "IRA1non" | M$gene ==  "IRA1mis"| M$gene ==  "GPB2"))

#This is simpler, just take sd across all reps and batches without nesting
BBm<-aggregate(value~barcode+gene, WI_Flask, FUN = mean)
BBs<-aggregate(value~barcode+gene, WI_Flask, FUN = sd)

#This is trickier, first need average acros batches
#Then take SD across batches
#BBm<-aggregate(value~barcode+gene, Bm, FUN = mean)
#BBs<-aggregate(value~barcode+gene, Bm, FUN = sd)
BBcv<-BBm
BBcv$value<-abs(BBs$value/BBm$value)
BBcv$mean<-BBm$value
BBcv$sd<-BBs$value
BBcv$n<-5
BBcv$CI<-1.96* ((BBcv$sd)/sqrt(BBcv$n))
BBcv$Percent_Error<-(BBcv$CI/BBcv$mean)*100
head(BBcv)

#Calculating 95% CIs for batch samples:
#Could do it simpler way, this is how we did it.
WI_Flask<-subset(M, M$gene ==  "IRA1non" | M$gene ==  "IRA1mis"| M$gene ==  "GPB2")
WI_Flask$Weighted_Means<-WI_Flask$value/WI_Flask$error
WI_Flask$Inverse_Variances<-1/WI_Flask$error
Sum_of_weighted_means_D<-aggregate(Weighted_Means~barcode+gene, WI_Flask, FUN = sum)
Sum_of_inverse_variances<-aggregate(Inverse_Variances~barcode+gene, WI_Flask, FUN = sum)


#This is trickier, first need average acros batches
#Sum_of_weighted_means<-subset(Sum_of_weighted_means, Sum_of_weighted_means$gene ==  "IRA1non" | Sum_of_weighted_means$gene ==  "IRA1mis"| Sum_of_weighted_means$gene ==  "GPB2")
#Sum_of_weighted_means<-subset(Sum_of_weighted_means, Sum_of_weighted_means$gene != "IRA1nonSpike" & Sum_of_weighted_means$gene != "IRA1misSpike")
#Sum_of_weighted_means$Weighted_Means<-Sum_of_weighted_means$IWM/Sum_of_weighted_means$IWV
#Sum_of_weighted_means$Inverse_Variances<-1/Sum_of_weighted_means$IWV
#Sum_of_weighted_means_D<-aggregate(Weighted_Means~barcode+gene, Sum_of_weighted_means, FUN = sum)
#Sum_of_inverse_variances<-aggregate(Inverse_Variances~barcode+gene, Sum_of_weighted_means, FUN = sum)

#Either way, finish it up here
Sum_of_weighted_means_D$Inverse_variances<-Sum_of_inverse_variances$Inverse_Variances
Sum_of_weighted_means_D$IWM<-Sum_of_weighted_means_D$Weighted_Means/Sum_of_weighted_means_D$Inverse_variances
Sum_of_weighted_means_D$IWV<-1/Sum_of_weighted_means_D$Inverse_variances
Sum_of_weighted_means_D$CI<-(1.96*(Sum_of_weighted_means_D$IWV))
Sum_of_weighted_means_D$Percent_Error<-(Sum_of_weighted_means_D$CI/Sum_of_weighted_means_D$IWM)*100
head(Sum_of_weighted_means_D)

#########################################
#Plot Figure 1D
M$BioRep<-factor(M$BioRep, levels = c("BB_A",  "BB_B", "BB_C",  "BB_D",  "X19_A", "X19_B", "X19_C", "X13_A", "X13_B", "X13_C", "X21_A", "X21_B", "X21_C", "X18_A", "X18_B", "X18_C", "X20_A", "X20_B", "X20_C", "X23_A", "X23_B", "X23_C", "X6_A",  "X6_B",  "X6_C",   "X3_A",  "X3_B",  "X3_C"), order = TRUE)

#This calculates the blue diamond
A<-subset(M, M$Type == "Original"  &(M$gene == "GPB2" | M$gene == "IRA1non"))
B<-aggregate(value~gene+Batch+Rep, A, FUN = mean)
C<-dcast(B, Batch+Rep~gene)
C$difference<-C$IRA1non - C$GPB2
C$BioRep<-paste(C$Batch, C$Rep, sep = "_")

#version with grey dots:
#F2A<-ggplot(subset(M, M$Type == "Original" & (M$gene == "GPB2" | M$gene == "IRA1non" | M$gene == "Unmutated")), aes(x = BioRep, y = value))+geom_point(data = subset(M, M$Type == "Original" & M$gene == "IRA1non"), color = "red", alpha = 0.2)+geom_boxplot(fill = NA, outlier.shape = NA, aes(color = gene))+theme_bw()+scale_color_manual(values = c("GPB2" = "green", "IRA1non" = "red", "Unmutated" = "black"))+geom_point(data = subset(M, M$Type == "Original"  & M$gene == "GPB2"), color = "green", alpha = 0.2, position = position_nudge(-0.25))+geom_point(data = subset(M, M$Type == "Original"  & M$gene == "Unmutated"), color = "black", alpha = 0.2, position = position_nudge(0.25))+ylab("Fitness advantage")+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 10.5)+geom_vline(xintercept = 13.5)+geom_vline(xintercept = 16.5)+geom_vline(xintercept = 19.5)+geom_vline(xintercept = 22.5)+geom_vline(xintercept = 25.5)+scale_x_discrete(labels = c("A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C"))+theme(axis.title.x = element_blank())+geom_point(data = C, aes(x = BioRep, y = difference), color = "blue", shape = 18, size = 2)+ geom_text(x=2.5, y=1.55, label="12/10/17")+ geom_text(x=6, y=1.48, label="5/9/15\nno ancestor")+ geom_text(x=9, y=1.55, label="5/1/15")+ geom_text(x=12, y=1.55, label="8/17/15")+ geom_text(x=15, y=1.55, label="5/9/15")+ geom_text(x=18, y=1.48, label="8/17/15\nno preculture")+ geom_text(x=21, y=1.55, label="9/4/15")+ geom_text(x=24, y=1.55, label="12/26/14")+ geom_text(x=27, y=1.55, label="11/24/14")+xlab("Replicates by batch")+ geom_segment(aes(x = 0.8, y = 0.92, xend = 4.2, yend = 0.92), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 4.8, y = 1.1, xend = 7.2, yend = 1.1), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 0.8, y = 0.22, xend = 4.2, yend = 0.22), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 4.8, y = 0.26, xend = 7.2, yend = 0.26), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 7.6, y = 1.24, xend = 8.4, yend = 1.24),  lty = 3, color = "darkgrey")+ geom_segment(aes(x = 9.6, y = 0.98, xend = 10.4, yend = 0.98), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 7.6, y = 0.26, xend = 8.4, yend = 0.26), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 9.6, y = 0.27, xend = 10.4, yend = 0.27), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 11.2, y = 0.65, xend = 11.2, yend = 1.21), lty = 3, color = "darkgrey")+ geom_segment(aes(x = 13.2, y = 0.76, xend = 13.2, yend = 0.99), lty = 3, color = "darkgrey")

#version with red/blue dots:
Figure_1<-ggplot(subset(M, M$Type == "Original" & (M$gene == "GPB2" | M$gene == "IRA1non" | M$gene == "Unmutated")), aes(x = BioRep, y = value))+geom_point(data = subset(M, M$Type == "Original" & M$gene == "IRA1non"), color = "royalblue3", alpha = 0.4, size = 0.5, shape = 16)+geom_boxplot(fill = NA, outlier.shape = NA, aes(color = gene))+theme_bw()+scale_color_manual(values = c("GPB2" = "forestgreen", "IRA1non" = "royalblue3", "Unmutated" = "black"))+geom_point(data = subset(M, M$Type == "Original"  & M$gene == "GPB2"), color = "forestgreen", alpha = 0.4, position = position_nudge(-0.25), size = 0.5, shape = 16)+geom_point(data = subset(M, M$Type == "Original"  & M$gene == "Unmutated"), color = "black", alpha = 0.4, position = position_nudge(0.25), size = 0.5, shape = 16)+ylab("Fitness advantage")+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 10.5)+geom_vline(xintercept = 13.5)+geom_vline(xintercept = 16.5)+geom_vline(xintercept = 19.5)+geom_vline(xintercept = 22.5)+geom_vline(xintercept = 25.5)+scale_x_discrete(labels = c("A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C"))+theme(axis.title.x = element_blank())+geom_point(data = C, aes(x = BioRep, y = difference), color = "mediumpurple1", shape = 18, size = 1.5)+ geom_text(x=2.5, y=1.55, label="12/10/17", size = 2.5)+ geom_text(x=6, y=1.55, label="5/9/15", size = 2.5)+ geom_text(x=6, y=1.45, label="(no ancestor)", size = 2.5)+ geom_text(x=9, y=1.55, label="5/1/15", size = 2.5)+ geom_text(x=12, y=1.55, label="8/17/15", size = 2.5)+ geom_text(x=15, y=1.55, label="5/9/15", size = 2.5)+ geom_text(x=18, y=1.55, label="8/17/15", size = 2.5)+ geom_text(x=18, y=1.45, label="(no preculture)", size = 2.5)+ geom_text(x=21, y=1.55, label="9/4/15", size = 2.5)+ geom_text(x=24, y=1.55, label="12/26/14", size = 2.5)+ geom_text(x=27, y=1.55, label="11/24/14", size = 2.5)+xlab("Replicates by batch")+ geom_segment(aes(x = 0.8, y = 0.92, xend = 4.2, yend = 0.92), lty = 3, size = 0.2, lwd = 1, color = "royalblue3")+ geom_segment(aes(x = 4.8, y = 1.1, xend = 7.2, yend = 1.1), lty = 3, size = 0.2, lwd = 1, color = "royalblue3")+theme(axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), legend.position = "none")

    #+ geom_segment(aes(x = 0.8, y = 0.22, xend = 4.2, yend = 0.22), lty = 3, size = 0.2, lwd = 1, color = "mediumpurple1")+ geom_segment(aes(x = 4.8, y = 0.26, xend = 7.2, yend = 0.26), lty = 3, size = 0.2, lwd = 1, color = "mediumpurple1")+ geom_segment(aes(x = 7.6, y = 1.24, xend = 8.4, yend = 1.24),  lty = 3, size = 0.2, lwd = 1, color = "royalblue3")+ geom_segment(aes(x = 9.6, y = 0.98, xend = 10.4, yend = 0.98), lty = 3, size = 0.2, lwd = 1, color = "royalblue3")+ geom_segment(aes(x = 7.6, y = 0.26, xend = 8.4, yend = 0.26), lty = 3, size = 0.2, lwd = 1, color = "mediumpurple1")+ geom_segment(aes(x = 9.6, y = 0.27, xend = 10.4, yend = 0.27), lty = 3, size = 0.2, lwd = 1, color = "mediumpurple1")


#########################################
#Figure 4C
#Start merging;You want the column order to be gene, value
Merge<-rbind(WIcv[,c(2,5)], Bcv[,c(3,6)], BBcv[,c(2,5)], SIMcv[,c(1,7)])
Merge$Type<-c(rep("Technical\nnoise", nrow(WIcv)), rep("Variation\nbetween\nreplicates", nrow(Bcv)), rep("Variation\nbetween\nbatches", nrow(BBcv)), rep("Sampling\nnoise", nrow(SIMcv)))
Merge$Type<-factor(Merge$Type, c("Sampling\nnoise", "Technical\nnoise","Variation\nbetween\nreplicates","Variation\nbetween\nbatches"), ordered = TRUE )

#This calculates the purple diamond
A<-subset(M, (M$gene == "GPB2" | M$gene == "IRA1non"))
B<-aggregate(value~gene+Batch+Rep, A, FUN = mean)
B<-aggregate(value~gene+Batch+Rep, A, FUN = sd)
C<-dcast(B, Batch+Rep~gene)
C$difference<-C$IRA1non - C$GPB2
C$BioRep<-paste(C$Batch, C$Rep, sep = "_")
D<-aggregate(difference~Batch, C, FUN = mean)
E<-aggregate(difference~Batch, C, FUN = sd)

#Track which are spike ins
#Merge$Spike<-FALSE
#Merge$Spike[grep("Spiked", Merge$gene)]<-TRUE

Merge$Split<-"None"
Merge$Split[which(Merge$Type == "Technical\nnoise" & Merge$gene == "GPB2")]<-"Yes"
Figure_4C<-ggplot(subset(Merge, Merge$gene != "IRA1mis"), aes(x = Type, y = sd))+geom_jitter(aes(color = gene, dodge = Split), size = 0.5, alpha = 0.7, shape = 16)+theme_bw()+theme(axis.title.x = element_blank())+ylab("Std. deviation on fitness")+geom_boxplot(fill = NA, outlier.color = NA, notch = TRUE)+scale_color_manual(values = c("GPB2" = "forestgreen", "IRA1non" = "royalblue3", "IRA1non*" = "red", "IRA1mis*" = "orange"))+theme(axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), legend.position = "none")

#+geom_point(x=4, y=sd(D$difference), color="mediumpurple1", shape=18, size=1.5)

#+geom_jitter(data = E, aes(x = 3, y = difference), color = "purple1", shape = 18, size = 1.5, alpha = 1)
#PLot
#ggplot(subset(Merge, Merge$gene != "IRA1mis"), aes(x = Type, y = sd))+geom_jitter(aes(color = gene), size = 0.5, alpha = 0.7, shape = 16)+theme_bw()+theme(axis.title.x = element_blank())+ylab("Std. deviation on fitness")+geom_boxplot(fill = NA, outlier.color = NA, notch = TRUE)+scale_color_manual(values = c("GPB2" = "forestgreen", "IRA1non" = "royalblue3", "IRA1non*" = "navyblue", "IRA1mis*" = "lightskyblue1"))+geom_point(x=4, y=sd(D$difference), color="purple1", shape=18, size=0.5)+geom_jitter(data = E, aes(x = 3, y = difference), color = "purple1", shape = 18, size = 0.75, alpha = 1)+theme(axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), legend.position = "none")

#Alternate figure 2B using the percent error derived from GrantFitSeq error:
#Merge<-rbind(Sum_of_weighted_means_B[,c(2,8)], Sum_of_weighted_means_C[,c(3,9)], Sum_of_weighted_means_D[,c(2,8)], Sum_of_weighted_means_A[,c(1,10)])
#Merge$Type<-c(rep("Within\neach\nflask", nrow(Sum_of_weighted_means_B)), rep("Between\nreplicates", nrow(Sum_of_weighted_means_C)), rep("Between\nbatches", nrow(Sum_of_weighted_means_D)), rep("Between\nrandom\nsamples", nrow(Sum_of_weighted_means_A)))
#F1A<-ggplot(Merge, aes(x = Type, y = Percent_Error))+geom_jitter(aes(color = Mutant), size = 0.5, alpha = 0.5)+theme_bw()+theme(axis.title.x = element_blank())+ylab("Percent error\non fitness measurement")+geom_boxplot(fill = NA, outlier.color = NA)+scale_color_manual(values = c("GPB2" = "forestgreen", "IRA1non" = "blue", "IRA1mis" = "red"))

##################################################
#Add a panel A that shows a trajectoru

#First load the frequency data:
load("/Users/ksamerot/Dropbox (ASU)/Samerotte Lab/Barcoding/1bigbatch/BC.Rfile")

M3<-BC
A<-subset(M3, M3$Batch == "1Big" & M3$Replicate == "A" & M3$Media == "M3" & DE == "DE1" & PCR == "PCRa")
Remove<-unique(subset(A, A$Freq < 0.00001)$Barcode)
A<-A[-which(A$Barcode %in% Remove),]
head(A)

F1A<-ggplot(A, aes(x = Time, y = Freq))+geom_line(size = 0.2, color = "grey", alpha = 0.2, aes(group = Barcode))+geom_line(size = 0.2, data = subset(A, gene == "neutral"), color = "black", alpha = 1, aes(group = Barcode))+geom_line(size = 0.2, data = subset(A, gene == "IRA1non"), color = "royalblue3", alpha = 1, aes(group = Barcode))+geom_line(size = 0.2, data = subset(A, gene == "GPB2"), color = "forestgreen", alpha = 1, aes(group = Barcode))+scale_y_log10()+theme_bw()+ylab("Barcode Frequency")+xlab("Timepoint")+theme(axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title.x = element_blank(), axis.text.x = element_blank())

Figure_4A<-ggplot(A, aes(x = Time, y = Freq))+geom_line(size = 0.2, color = "grey", alpha = 0.2, aes(group = Barcode))+geom_line(data = subset(A, gene == "neutral"), color = "black", alpha = 1, aes(group = Barcode))+geom_line(size = 0.2, data = subset(A, gene == "IRA1misSpike"), color = "orange", aes(group = Barcode))+geom_line(size = 0.2, data = subset(A, gene == "IRA1nonSpike"), aes(group = Barcode), color = "red")+scale_y_log10()+theme_bw()+ylab("Barcode Frequency")+xlab("Timepoint")+theme(axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title.x = element_blank(), axis.text.x = element_blank())

A<-subset(M3, M3$Batch == "1Big" & (M3$Replicate == "A" | M3$Replicate == "D") & M3$Media == "M3" & DE == "DE1" & PCR == "PCRa")
Remove<-unique(subset(A, A$Freq < 0.00001)$Barcode)
A<-A[-which(A$Barcode %in% Remove),]
head(A)

Figure_4B<- ggplot(A, aes(x = Time, y = Freq))+geom_line(size = 0.2, color = "grey", alpha = 0.2, aes(group = Barcode))+geom_line(size = 0.2, data = subset(A, gene == "neutral"), color = "black", alpha = 1, aes(group = Barcode))+geom_line(size = 0.5, data = subset(A, Barcode == "7836"), color = "royalblue3", alpha = 1, aes(group = Barcode))+scale_y_log10()+theme_bw()+ylab("Barcode Frequency")+xlab("Timepoint")+theme(strip.text = element_blank(), strip.background = element_blank(), axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 7))+facet_grid(Replicate~.)

