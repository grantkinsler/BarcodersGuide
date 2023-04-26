load("/Users/ksamerot/Dropbox (ASU)/Samerotte Lab/Barcoding/1bigbatch/BC.Rfile")
library(reshape2)
library(ggplot2)

BC_cycles<-BC


#Add column labeling PCR/Cycles
BC_cycles$PCR_Cycles <- paste(BC_cycles$PCR, BC_cycles$Cycles, sep="_")

#Make dcast table, splitting by PCR/Cycles
BC_cycles_dcast <- dcast(BC_cycles, ID+Barcode+Time+Replicate+DE~PCR_Cycles, value.var="Freq")

#Track coverage
BC_cycles_coverage <- dcast(BC_cycles, ID+Barcode+Time+Replicate+DE~PCR_Cycles, value.var="Coverage")
####################################################################

#Comparisons of 22 v 27 cycles (deleted comparisons with 0 results/unique IDs)

#Compare PCRa_22 to PCRb_27

#Delete columns besides PCRa_22 and PCRb_27
BC_cycles_A22.B27 <- BC_cycles_dcast[,-c(7,8,10,11)]
BC_cycles_coverage_A22.B27 <- BC_cycles_coverage[,-c(7,8,10,11)]

#Remove any rows that don't have both PCRa_22 and PCRb_27
BC_cycles_A22.B27 <- na.omit(BC_cycles_A22.B27)
BC_cycles_coverage_A22.B27 <- na.omit(BC_cycles_coverage_A22.B27)

#Check how many unique IDs (comparisons) there are
length(unique(BC_cycles_A22.B27$ID)) #Results = 2 unique IDs

#Add Type column to label which comparison was done
BC_cycles_A22.B27$Type <- "A22vB27"

#Add Unique column to label ID+comparison (because PCR/Cycles columns will be renamed)
BC_cycles_A22.B27$Unique <- paste(BC_cycles_A22.B27$ID, BC_cycles_A22.B27$Type, sep="-")

#Rename PCR/Cycles columns so that rbind can be used later
colnames(BC_cycles_A22.B27)[6]<-"PCR1"
colnames(BC_cycles_A22.B27)[7]<-"PCR2"

#Add coverage
BC_cycles_A22.B27$Min_Coverage<-apply(BC_cycles_coverage_A22.B27[,6:7], 1, FUN = min)

head(BC_cycles_A22.B27)
#Repeat for each comparison ##############################################

BC_cycles_A27.B22 <- BC_cycles_dcast[,-c(6,9,10,11)]
BC_cycles_coverage_A27.B22 <- BC_cycles_coverage[,-c(6,9,10,11)]
BC_cycles_A27.B22 <- na.omit(BC_cycles_A27.B22)
BC_cycles_coverage_A27.B22 <- na.omit(BC_cycles_coverage_A27.B22)
length(unique(BC_cycles_A27.B22$ID)) #Results = 15 unique IDs
BC_cycles_A27.B22$Type <- "A27vB22"
BC_cycles_A27.B22$Unique <- paste(BC_cycles_A27.B22$ID, BC_cycles_A27.B22$Type, sep="-")
colnames(BC_cycles_A27.B22)[6]<-"PCR1"
colnames(BC_cycles_A27.B22)[7]<-"PCR2"
BC_cycles_A27.B22$Min_Coverage<-apply(BC_cycles_coverage_A27.B22[,6:7], 1, FUN = min)

BC_cycles_A27.B27 <- BC_cycles_dcast[,-c(6,8,10,11)]
BC_cycles_coverage_A27.B27 <- BC_cycles_coverage[,-c(6,8,10,11)]
BC_cycles_A27.B27 <- na.omit(BC_cycles_A27.B27)
BC_cycles_coverage_A27.B27 <- na.omit(BC_cycles_coverage_A27.B27)
length(unique(BC_cycles_A27.B27$ID)) #Results = 18 unique IDs
BC_cycles_A27.B27$Type <- "A27vB27"
BC_cycles_A27.B27$Unique <- paste(BC_cycles_A27.B27$ID, BC_cycles_A27.B27$Type, sep="-")
colnames(BC_cycles_A27.B27)[6]<-"PCR1"
colnames(BC_cycles_A27.B27)[7]<-"PCR2"
BC_cycles_A27.B27$Min_Coverage<-apply(BC_cycles_coverage_A27.B27[,6:7], 1, FUN = min)

BC_cycles_A27.C22 <- BC_cycles_dcast[,-c(6,8,9,11)]
BC_cycles_coverage_A27.C22 <- BC_cycles_coverage[,-c(6,8,9,11)]
BC_cycles_A27.C22 <- na.omit(BC_cycles_A27.C22)
BC_cycles_coverage_A27.C22 <- na.omit(BC_cycles_coverage_A27.C22)
length(unique(BC_cycles_A27.C22$ID)) #Results = 4 unique IDs
BC_cycles_A27.C22$Type <- "A27vC22"
BC_cycles_A27.C22$Unique <- paste(BC_cycles_A27.C22$ID, BC_cycles_A27.C22$Type, sep="-")
colnames(BC_cycles_A27.C22)[6]<-"PCR1"
colnames(BC_cycles_A27.C22)[7]<-"PCR2"
BC_cycles_A27.C22$Min_Coverage<-apply(BC_cycles_coverage_A27.C22[,6:7], 1, FUN = min)

BC_cycles_A27.C27 <- BC_cycles_dcast[,-c(6,8,9,10)]
BC_cycles_coverage_A27.C27 <- BC_cycles_coverage[,-c(6,8,9,10)]
BC_cycles_A27.C27 <- na.omit(BC_cycles_A27.C27)
BC_cycles_coverage_A27.C27 <- na.omit(BC_cycles_coverage_A27.C27)
length(unique(BC_cycles_A27.C27$ID)) #Results = 1 unique ID
BC_cycles_A27.C27$Type <- "A27vC27"
BC_cycles_A27.C27$Unique <- paste(BC_cycles_A27.C27$ID, BC_cycles_A27.C27$Type, sep="-")
colnames(BC_cycles_A27.C27)[6]<-"PCR1"
colnames(BC_cycles_A27.C27)[7]<-"PCR2"
BC_cycles_A27.C27$Min_Coverage<-apply(BC_cycles_coverage_A27.C27[,6:7], 1, FUN = min)

BC_cycles_B27.C22 <- BC_cycles_dcast[,-c(6,7,8,11)]
BC_cycles_coverage_B27.C22 <- BC_cycles_coverage[,-c(6,7,8,11)]
BC_cycles_B27.C22 <- na.omit(BC_cycles_B27.C22)
BC_cycles_coverage_B27.C22 <- na.omit(BC_cycles_coverage_B27.C22)
length(unique(BC_cycles_B27.C22$ID)) #Results = 4 unique IDs
BC_cycles_B27.C22$Type <- "B27vC22"
BC_cycles_B27.C22$Unique <- paste(BC_cycles_B27.C22$ID, BC_cycles_B27.C22$Type, sep="-")
colnames(BC_cycles_B27.C22)[6]<-"PCR1"
colnames(BC_cycles_B27.C22)[7]<-"PCR2"
BC_cycles_B27.C22$Min_Coverage<-apply(BC_cycles_coverage_B27.C22[,6:7], 1, FUN = min)

BC_cycles_B27.C27 <- BC_cycles_dcast[,-c(6,7,8,10)]
BC_cycles_coverage_B27.C27 <- BC_cycles_coverage[,-c(6,7,8,10)]
BC_cycles_B27.C27 <- na.omit(BC_cycles_B27.C27)
BC_cycles_coverage_B27.C27 <- na.omit(BC_cycles_coverage_B27.C27)
length(unique(BC_cycles_B27.C27$ID)) #Results = 1 unique ID
BC_cycles_B27.C27$Type <- "B22vC27"
BC_cycles_B27.C27$Unique <- paste(BC_cycles_B27.C27$ID, BC_cycles_B27.C27$Type, sep="-")
colnames(BC_cycles_B27.C27)[6]<-"PCR1"
colnames(BC_cycles_B27.C27)[7]<-"PCR2"
BC_cycles_B27.C27$Min_Coverage<-apply(BC_cycles_coverage_B27.C27[,6:7], 1, FUN = min)

#Number comparisons: A22.B27 = 2, A27.B22 = 15, A27.B27 = 18, A27.C22 = 4, A27.C27 = 1, B27.C22 = 4, B27.C27 = 1
#7 types of comparisons

#############################################################################

#Calculate rsquared values for each comparison
h<-function(ID) summary(lm(PCR1~PCR2, data = ID))$r.squared
a<-lapply(split(BC_cycles_A22.B27, BC_cycles_A22.B27$Unique), h)
b<-lapply(split(BC_cycles_A27.B22, BC_cycles_A27.B22$Unique), h)
c<-lapply(split(BC_cycles_A27.B27, BC_cycles_A27.B27$Unique), h)
d<-lapply(split(BC_cycles_A27.C22, BC_cycles_A27.C22$Unique), h)
e<-lapply(split(BC_cycles_A27.C27, BC_cycles_A27.C27$Unique), h)
f<-lapply(split(BC_cycles_B27.C22, BC_cycles_B27.C22$Unique), h)
g<-lapply(split(BC_cycles_B27.C27, BC_cycles_B27.C27$Unique), h)

#Transpose rsquared tables
R2_a<-as.data.frame(t(as.data.frame(a)))
R2_b<-as.data.frame(t(as.data.frame(b)))
R2_c<-as.data.frame(t(as.data.frame(c)))
R2_d<-as.data.frame(t(as.data.frame(d)))
R2_e<-as.data.frame(t(as.data.frame(e)))
R2_f<-as.data.frame(t(as.data.frame(f)))
R2_g<-as.data.frame(t(as.data.frame(g)))

#Add type column labeling which comparison was made (wrote in 22 v 27 order here so final graph only has 22 v 27 and 27 v 27)
R2_a$Type<-"A22 v B27"
R2_b$Type<-"B22 v A27"
R2_c$Type<-"A27 v B27"
R2_d$Type<-"C22 v A27"
R2_e$Type<-"A27 v C27"
R2_f$Type<-"C22 v B27"
R2_g$Type<-"B27 v C27"

#Add min coverage
R2_a$MC<-aggregate(Min_Coverage~Unique, BC_cycles_A22.B27, FUN = mean)$Min_Coverage
R2_b$MC<-aggregate(Min_Coverage~Unique, BC_cycles_A27.B22, FUN = mean)$Min_Coverage
R2_c$MC<-aggregate(Min_Coverage~Unique, BC_cycles_A27.B27, FUN = mean)$Min_Coverage
R2_d$MC<-aggregate(Min_Coverage~Unique, BC_cycles_A27.C22, FUN = mean)$Min_Coverage
R2_e$MC<-aggregate(Min_Coverage~Unique, BC_cycles_A27.C27, FUN = mean)$Min_Coverage
R2_f$MC<-aggregate(Min_Coverage~Unique, BC_cycles_B27.C22, FUN = mean)$Min_Coverage
R2_g$MC<-aggregate(Min_Coverage~Unique, BC_cycles_B27.C27, FUN = mean)$Min_Coverage

#rbind all rsquared tables together
allPCR_cycles <- rbind(R2_a,R2_b,R2_c,R2_d,R2_e,R2_f,R2_g)

#plot of all rsquareds (1 boxplot)
ggplot(allPCR_cycles, aes(x=1, y=V1))+geom_jitter(aes(color=Type))+geom_boxplot(fill="NA", outlier.color = "NA")+theme_bw()+theme(axis.text.x = element_text(angle=90, size=10))

#plot of all comparisons on separate boxplots
ggplot(allPCR_cycles, aes(x=Type, y=V1))+geom_jitter(aes(color=Type))+geom_boxplot(fill="NA", outlier.color = "NA")+theme_bw()+theme(axis.text.x = element_text(angle=90, size=10))

#Seems clear the low coverage ones have worse R2:
ggplot(allPCR_cycles, aes(x = MC, y = V1))+geom_point()+scale_x_log10()
######################################################################################

#rbind frequency tables together
allPCR_cycles_freq <- rbind(BC_cycles_A22.B27,BC_cycles_A27.B22,BC_cycles_A27.B27,BC_cycles_A27.C22,BC_cycles_A27.C27,BC_cycles_B27.C22,BC_cycles_B27.C27)

#frequencies graph with multiple plots (one for each unique ID)
ggplot(allPCR_cycles_freq, aes(x = PCR1, y = PCR2))+geom_point(alpha = 0.1)+geom_abline(slope = 1, intercept = 0)+theme_bw()+xlab("Barcode frequencies with PCR1")+ylab("Barcode frequencies\nin with PCR2")+facet_wrap(~Unique)+theme(axis.text.x=element_text(angle=90))+scale_y_log10()+scale_x_log10()

F3A_lower<-ggplot(subset(allPCR_cycles_freq, Unique == "A4-A27vB22"), aes(x = PCR1, y = PCR2))+geom_point(alpha = 0.1)+geom_abline(slope = 1, intercept = 0)+theme_bw()+xlab("Barcode freq\nPCR rep 2")+ylab("Barcode freq\nPCR rep 1")+theme(axis.text.x=element_text(angle=90, size = 6), axis.text.y=element_text(size = 6), axis.title.y=element_text(size = 8), , axis.title.x=element_text(size = 8))+scale_y_log10()+scale_x_log10()

#frequencies graph (no log scale to show top frequencies better)
ggplot(allPCR_cycles_freq, aes(x = PCR1, y = PCR2))+geom_point()+geom_abline(slope = 1, intercept = 0)+theme_bw()+xlab("Barcode frequencies with PCR1")+ylab("Barcode frequencies\nin with PCR2")+facet_wrap(~Unique)+theme(axis.text.x=element_text(angle=90))

#####################################################################################

#add in a Cycles column to compare 22 v 27 and 27 v 27

#Pull out first cycle type
allPCR_cycles$Cycles<-substr(allPCR_cycles$Type,2,6)

#Pull out second cycle type
allPCR_cycles$Cycles2<-substr(allPCR_cycles$Type,8,9)

#Paste together into new column
allPCR_cycles$Cycles3<-paste(allPCR_cycles$Cycles, allPCR_cycles$Cycles2, sep="")

#Remove 1st and 2nd cycle type columns
allPCR_cycles<-allPCR_cycles[,-c(4,5)]

#Rename new column as "Cycles"
colnames(allPCR_cycles)[4]<-"Cycles"

#Add new "Category" column w/ same values as Unique (used later to label 2 boxplot graph)
allPCR_cycles$Category <- "Cycles"

head(allPCR_cycles)

#plot of 22 v 27 and 27 v 27 rsquared values (2 box plots)
ggplot(allPCR_cycles,aes(x=Cycles,y=V1))+geom_jitter(aes(color=Cycles))+ylab("R squared")+geom_boxplot(notch = FALSE, fill = NA, outlier.shape = NA)+theme(axis.text.x=element_text(angle = 90))



#plot of rsquared values for cycles (1 box plot)
ggplot(allPCR_cycles,aes(x=1,y=V1))+geom_jitter(aes(color=Cycles))+ylab("R squared")+geom_boxplot(notch = FALSE, fill = NA, outlier.shape = NA)+theme(axis.text.x=element_text(angle = 90))

#####################################################################################

#DNA Extractions

#dcast table of BC to split out DNA extractions
BC_DE <- dcast(BC, ID+Barcode+Time+Replicate+PCR~DE, value.var = "Freq")
BC_DE_coverage<-dcast(BC, ID+Barcode+Time+Replicate+PCR~DE, value.var = "Coverage")
#####################################################################################

#Comparisons of DNA extractions using DE1 as baseline (compare all to DE1)

#####################################################################################

#DE1 v DE2

#remove all DE extractions but DE1 and DE2
BC_DE2 <- BC_DE[,-c(8,9,10)]
BC_DE2_coverage <- BC_DE_coverage[,-c(8,9,10)]

#remove any rows without both DE1 and DE2
BC_DE2 <- na.omit(BC_DE2)
BC_DE2_coverage <- na.omit(BC_DE2_coverage)

#Check how many unique IDs there are (tells us how many comparisons can be made between DE2 and DE1)
length(unique(BC_DE2$ID))

#add a column to label which comparison
BC_DE2$Type <- "DE1vDE2"

#each table will have either DE2, DE3, etc. columns --> need to rename DE column because rbind needs same column names
colnames(BC_DE2)[7]<-"DE_other"
BC_DE2$Unique <- paste(BC_DE2$ID, BC_DE2$Type, sep="-")
BC_DE2$Min_Coverage<-apply(BC_DE2_coverage[,6:7], 1, FUN = min)
head(BC_DE2)


#repeat steps above for DE3, DE4, and DE5 (name the tables differently, need to remake BC_DE each time)

#DE1 v DE3
#BC_DE <- dcast(BC, ID+Barcode+Time+Replicate+PCR~DE, value.var = "Freq")
BC_DE3 <- BC_DE[,-c(7,9,10)]
BC_DE3_coverage <- BC_DE_coverage[,-c(7,9,10)]
BC_DE3 <- na.omit (BC_DE3)
BC_DE3_coverage <- na.omit (BC_DE3_coverage)
length(unique(BC_DE3$ID))
BC_DE3$Type <- "DE1vDE3"
colnames(BC_DE3)[7]<-"DE_other"
BC_DE3$Unique <- paste(BC_DE3$ID, BC_DE3$Type, sep="-")
BC_DE3$Min_Coverage<-apply(BC_DE3_coverage[,6:7], 1, FUN = min)

#DE1 v DE4
#BC_DE <- dcast(BC, ID+Barcode+Time+Replicate+PCR~DE, value.var = "Freq")
BC_DE4 <- BC_DE[,-c(7,8,10)]
BC_DE4_coverage <- BC_DE_coverage[,-c(7,8,10)]
BC_DE4 <- na.omit (BC_DE4)
BC_DE4_coverage <- na.omit (BC_DE4_coverage)
length(unique(BC_DE4$ID))
BC_DE4$Type <- "DE1vDE4"
colnames(BC_DE4)[7]<-"DE_other"
BC_DE4$Unique <- paste(BC_DE4$ID, BC_DE4$Type, sep="-")
BC_DE4$Min_Coverage<-apply(BC_DE4_coverage[,6:7], 1, FUN = min)

#DE1 v DE5
#BC_DE <- dcast(BC, ID+Barcode+Time+Replicate+PCR~DE, value.var = "Freq")
BC_DE5 <- BC_DE[,-c(7,8,9)]
BC_DE5_coverage <- BC_DE_coverage[,-c(7,8,9)]
BC_DE5 <- na.omit (BC_DE5)
BC_DE5_coverage <- na.omit (BC_DE5_coverage)
length(unique(BC_DE5$ID))
BC_DE5$Type <- "DE1vDE5"
colnames(BC_DE5)[7]<-"DE_other"
BC_DE5$Unique <- paste(BC_DE5$ID, BC_DE5$Type, sep="-")
BC_DE5$Min_Coverage<-apply(BC_DE5_coverage[,6:7], 1, FUN = min)

#Number of unique IDs = 19 + 8 + 6 + 1 = 34

#####################################################################################

#Calculate rsquared values
h<-function(ID) summary(lm(DE1~DE_other, data = ID))$r.squared
a<-lapply(split(BC_DE2, BC_DE2$Unique), h)
b<-lapply(split(BC_DE3, BC_DE3$Unique), h)
c<-lapply(split(BC_DE4, BC_DE4$Unique), h)
d<-lapply(split(BC_DE5, BC_DE5$Unique), h)

#Transpose data tables
R2_a<-as.data.frame(t(as.data.frame(a)))
R2_b<-as.data.frame(t(as.data.frame(b)))
R2_c<-as.data.frame(t(as.data.frame(c)))
R2_d<-as.data.frame(t(as.data.frame(d)))

#Type column to label which comparison
R2_a$Type<-"DE1 v DE2"
R2_b$Type<-"DE1 v DE3"
R2_c$Type<-"DE1 v DE4"
R2_d$Type<-"DE1 v DE5"

#Add min coverage
R2_a$MC<-aggregate(Min_Coverage~Unique, BC_DE2, FUN = mean)$Min_Coverage
R2_b$MC<-aggregate(Min_Coverage~Unique, BC_DE3, FUN = mean)$Min_Coverage
R2_c$MC<-aggregate(Min_Coverage~Unique, BC_DE4, FUN = mean)$Min_Coverage
R2_d$MC<-aggregate(Min_Coverage~Unique, BC_DE5, FUN = mean)$Min_Coverage

#rbind the rsquared tables together
allDE <- rbind(R2_a,R2_b,R2_c,R2_d)

#plot of rsquared values for each DE comparison (one boxplot per comparison)
ggplot(allDE, aes(x=Type, y=V1))+geom_jitter(aes(color=Type))+geom_boxplot(fill="NA")

#plot of rsquared values for comparisons all one 1 boxplot
ggplot(allDE, aes(x=1, y=V1))+geom_jitter(aes(color=Type))+geom_boxplot(fill="NA", outlier.color = "NA")+theme_bw()+theme(axis.text.x = element_text(angle=90, size=10))

#Coverage
ggplot(allDE, aes(x = MC, y = V1))+geom_point()+scale_x_log10()

#####################################################################################

#rbind the frequency tables together
DE_freq <- rbind(BC_DE2, BC_DE3, BC_DE4, BC_DE5)
DE_freq$Unique <- paste(DE_freq$ID, DE_freq$Type, sep="-")

#frequency graph with multiple plots (one for each unique ID)
ggplot(DE_freq, aes(x = DE1, y = DE_other))+geom_point(alpha = 0.1)+geom_abline(slope = 1, intercept = 0)+theme_bw()+xlab("Barcode frequencies with DE1")+ylab("Barcode frequencies\nin with DE other")+facet_wrap(~Unique)+theme(axis.text.x=element_text(angle=90))+scale_y_log10()+scale_x_log10()

F3B_lower<-ggplot(subset(DE_freq, DE_freq$Unique == "Y4-DE1vDE4"), aes(x = DE1, y = DE_other))+geom_point(alpha = 0.1)+geom_abline(slope = 1, intercept = 0)+theme_bw()+xlab("Barcode freq\nExtraction rep 2")+ylab("Barcode freq\nExtraction rep 1")+theme(axis.text.x=element_text(angle=90, size = 6), axis.text.y=element_text(size = 6), axis.title.y=element_text(size = 8), , axis.title.x=element_text(size = 8))+scale_y_log10()+scale_x_log10()

#frequency graph without log scale (to see top frequencies better)
ggplot(DE_freq, aes(x = DE1, y = DE_other))+geom_point()+geom_abline(slope = 1, intercept = 0)+theme_bw()+xlab("Barcode frequencies with DE1")+ylab("Barcode frequencies\nin with DE other")+facet_wrap(~Unique)+theme(axis.text.x=element_text(angle=90))

#####################################################################################

#add in a DE column

#Pull out first DE Type
allDE$DE<-substr(allDE$Type,1,6)

#Pull out second DE Type
allDE$DE2<-substr(allDE$Type,7,9)

#Paste together to label comparison
allDE$DE3<-paste(allDE$DE, allDE$DE2, sep="")
allDE<-allDE[,-c(4,5)]
colnames(allDE)[4]<-"DE"

#Renaming/labeling columns for rbinding
head(allDE)
head(allPCR_cycles)

#Unique column to label Cycles column
allPCR_cycles$Unique <- paste(allPCR_cycles$Cycles)
allPCR_cycles <- allPCR_cycles[,-c(4)]

#Unique column to label DE column
allDE$Unique <- paste(allDE$DE)
allDE <- allDE[,-c(4)]

#Add Comparison column to label as DE (used to make/label 2 boxplot graph)
allDE$Category <- "DE"

#Rbind rsquared tables together
Comparison <- rbind(allPCR_cycles, allDE)

####################  GRAPHS  ##########################################################

#plot of rsquared values for each comparison (PCR Cycles & DNA extraction) #Use "Unique" column to plot only cycles, not PCRa vs PCRb
ggplot(Comparison,aes(x=Unique,y=V1))+geom_jitter(aes(color=Unique))+ylab("R squared")+xlab("Comparison Type")+geom_boxplot(notch = FALSE, fill = NA, outlier.shape = NA)+theme(axis.text.x=element_text(angle = 90))

#plot of rsquared values (2 boxplots, one for Cycles and one for DE)
ggplot(Comparison,aes(x=Category,y=V1))+geom_jitter(aes(color=Unique))+ylab("R squared")+xlab("Comparison Type")+geom_boxplot(notch = FALSE, fill = NA, outlier.shape = NA)

#plot of rsquared values (1 boxplot)
ggplot(Comparison,aes(x=1,y=V1))+geom_jitter(aes(color=Unique))+ylab("R squared")+geom_boxplot(notch = FALSE, fill = NA, outlier.shape = NA)+theme(axis.text.x=element_text(angle = 90))


Comparison$Coverage<-"Low"
Comparison$Coverage[which(Comparison$MC>10000)]<-"Mid"
Comparison$Coverage[which(Comparison$MC>300000)]<-"High"

Comparison$Unique[which(Comparison$Unique == "DE1 v DE2")]<-"no modification"
Comparison$Unique[which(Comparison$Unique == "27 v 27")]<-"no modification"
Comparison$Unique[which(Comparison$Unique == "22 v 27")]<-"fewer cycles"
Comparison$Unique[which(Comparison$Unique == "DE1 v DE3")]<-"(-) beads"
Comparison$Unique[which(Comparison$Unique == "DE1 v DE4")]<-"(+) phenol"
Comparison$Unique[which(Comparison$Unique == "DE1 v DE5")]<-"no modification"

Comparison$Coverage<-factor(Comparison$Coverage, levels = c("Low", "Mid", "High"), ordered = TRUE)

Comparison$Unique<-factor(Comparison$Unique, levels = c("no modification", "fewer cycles", "(-) beads", "(+) phenol"), ordered = TRUE)

Comparison$Modification<-Comparison$Unique
F3C<-ggplot(Comparison,aes(x=Category,y=V1))+geom_jitter(aes(color=Coverage, shape = Modification), size = 1)+ylab("R squared")+xlab("Comparison Type")+geom_boxplot(notch = FALSE, fill = NA, outlier.shape = NA)+theme_bw()+scale_x_discrete(labels = c("PCR\n replicates", "Extraction\nreplicates"))+theme(axis.title.x = element_blank())+ylab(expression(paste("Reproducibility ", "(", R^2, ")")))+scale_shape_manual(values = c(15,16,17,18))+scale_color_manual(values = c("red2", "black", "royalblue1"))

########################
#Neeed figure 3D
#First go back to the BC table
#Pull out replicate B from 1BB, just DE1 PCRa
B_Compare<-subset(BC, BC$Replicate == "B" & BC$PCR == "PCRa" & (BC$DE == "DE1" | BC$DE == "DE2")& BC$Batch == "1Big" & BC$Media == "M3")

A<-subset(B_Compare, Freq < 0.00001)
Remove<-unique(A$Barcode)
B_Compare<-B_Compare[-which(B_Compare$Barcode %in% Remove),]
F3D_traj<-ggplot(B_Compare, aes(x = Time, y = Freq))+geom_line(alpha = 0.1, color = "grey", aes(group = Barcode))+scale_y_log10()+theme_bw()+facet_grid(.~DE)+geom_line(data = subset(B_Compare, gene == "IRA1non" | gene == "GPB2" | gene == "neutral"), aes(group = Barcode, color = gene), alpha = 0.4)+scale_color_manual(values = c("neutral" = "black", "IRA1non" = "royalblue3", "GPB2" = "forestgreen"))+ylab("Barcode frequency")+theme(axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 5), axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 5), strip.background = element_blank(), strip.text.x = element_blank(), legend.title = element_text(size = 6), legend.text  = element_text(size = 5))+xlab("Timepoint")

#########
#Calculate fitness
B_Compare$Time<-paste("T", B_Compare$Time, sep = "")
B_Compare$Label<-paste(B_Compare$DE, "R1", B_Compare$Time, sep = "-")
Input<-dcast(B_Compare, Barcode~Label, value.var = "Reads")
colnames(Input)[1]<-"barcode"
colnames(Input)
Input$duplicate<-Input[,5]
colnames(Input)[11]<-"DE2-R1-T3"
head(Input)

write.csv(Input, file = "/users/ksamerot/Dropbox (ASU)/FitnessAssay_forKerry/data/Fig3D_Noise.csv", row.names = FALSE)

#load this data
read.csv("/users/ksamerot/Dropbox (ASU)/FitnessAssay_forKerry/data/Noise_F3D_Output.csv")->A
A$gene<-B_Compare$gene[match(A$barcode, B_Compare$Barcode)]
A<-A[,c(1,3,5,10)]
A<-melt(A, id.var = c("barcode", "gene"))

F3D_right<-ggplot(subset(A, gene == "neutral" | gene == "IRA1non" | gene == "GPB2" ), aes(x = variable, y = value, color = gene))+theme_bw()+scale_color_manual(values = c("neutral" = "black", "IRA1non" = "royalblue3", "GPB2" = "forestgreen"))+geom_boxplot(fill = NA, aes(color = gene), outlier.shape = NA)+scale_x_discrete(labels = c("Extraction\nRep 1", "Extraction\nRep 2"))+geom_jitter(data = subset(A, gene == "neutral"), color = "black", size = 0.5, position = position_nudge(0.2))+geom_point(data = subset(A, gene == "GPB2"), color = "forestgreen", size = 0.5, position = position_nudge(-0.2))+geom_point(data = subset(A, gene == "IRA1non"), color = "royalblue3", size = 0.5)+theme(axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 5), axis.title.x = element_blank(), axis.text.x = element_text(size = 5),legend.position = "none")+ylab("Fitness")

#################################################
# Assemble figure 3
pdf("/users/ksamerot/Dropbox (ASU)/Noise Paper/Figure3_06092022_test.pdf", width = 8, height = 6)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
library(gridExtra)
library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 8, heights = unit(c(1,1,1,1,1, 1), c( "lines", "null", "lines", "null", "lines", "null")), widths = unit(c(1,1,1,1,1,1,1,1), c("lines", "null", "lines", "null", "lines", "null", "lines", "null")))))


print(F3C, vp = vplayout(2:4, 6:8))
print(F3B_lower, vp = vplayout(4, 4))
print(F3A_lower, vp = vplayout(4, 2))
print(F3D_traj, vp = vplayout(6, 2:6))
print(F3D_right, vp = vplayout(6, 8))

# Add panel letters
grid.text("A", vp = vplayout(1,1), gp = gpar(fontsize = 12), vjust = 1)
grid.text("B", vp = vplayout(1,3), gp = gpar(fontsize = 12), vjust = 1)
grid.text("C", vp = vplayout(1,5), gp = gpar(fontsize = 12), vjust = 1)
grid.text("D", vp = vplayout(5,1), gp = gpar(fontsize = 12), vjust = 1)

# close pdf
dev.off()
