setwd()

write("Info for article:", "info.csv") #write some information to a file that you will use for the article.
version<-R.Version() #Fist information which R version are you using?
str(version)
write(version$version.string, "info.csv", append=T) #Write it to the file.

#Read in the variables
variables<-read.csv("Ex_In_Variables.csv", row.names=1, header=T) #Read in your file with the variables per samples
View(variables)

#######################################
######START ANALYSIS OF DATA!
write("statistical results:", "statistics.txt")


################ GENERA
### ALPHA Diversity
genera<-read.csv("Ex_In_Genera.csv", row.names=1)
View(genera)
library("vegan") #vegan has a very good tutorial pdf /documentation
richness<-specnumber(genera)  
shannon<-diversity(genera)
evenness<- shannon/log(specnumber(genera))

richnessA<-specnumber(genera[1:8,])
richnessB<-specnumber(genera[9:16,])

par(mfrow=c(3,1)) #split the plot window in 3
plot(richness,xlab="Samples")
plot(shannon, xlab="Samples")
plot(evenness, xlab="Samples")
dev.off()
Alpha<-cbind(richness,shannon,evenness)
write.csv(Alpha, "alpha.csv") #this is your second result

library(car)
aovT<-aov(log(richness) ~ Plant/Date/Disinfection,data=variables)
summary(aovT)
testT<-t.test(log(richnessA),log(richnessB), paried=T)
testT

library(performance)
check_model(aovT)
OUT<-capture.output(summary(aovT))
write(OUT, "statistics.txt", append=T)
outT<-capture.output(testT)
write(outT, "statistics.txt", append=T)

plot(richness~as.factor(variables$Origin), ylab="", xlab="",type="p")
plot(richness)

###BETADIVERSITY:
beta<-vegdist(genera, method="bray", diag=T, upper=T)
betaA<-vegdist(genera[1:8,], method="bray", diag=T, upper=T)
betaB<-vegdist(genera[9:16,], method="bray", diag=T, upper=T)

plot(hclust(beta, method="average"), main="", xlab="", hang=-1, cex=0.7)
title(main="a)", adj = 0,cex.main=1.5, mgp = c(0, 0, 0))

#PERMANOVA analysis: How much of the variance in betadiversity is explained by different variables
ado<-adonis(beta ~ variables$Plant/variables$Date/variables$Disinfection) 
ado
adonisOUT<-capture.output(ado)
write(adonisOUT, "statistics.txt", append=T)

adoT<-adonis(beta ~ variables$Origin, permutations = 255,strata=variables$Pair)
adoT
adoOutT<-capture.output(adoT)
write(adoOutT, "statistics.txt", append=T)

mtl<-mantel(betaA,betaB,method="spearman") #Test whether the two distance matrixes are correlated
mtl

write(mtl$statistic, "statistics.csv", append=T)
write(mtl$signif, "statistics.csv", append=T)


################ ANTIBIOTIC RESISTOME
### ALPHA Diversity
args<-read.csv("Ex_In_ARGs.csv", row.names=1)
View(args)
richness2<-specnumber(args) 
shannon2<-diversity(args)
evenness2<- shannon2/log(specnumber(args))

richnessA2<-specnumber(args[1:8,])
richnessB2<-specnumber(args[9:16,])

par(mfrow=c(3,1)) #split the plot window in 3
plot(richness2,xlab="Samples")
plot(shannon2, xlab="Samples")
plot(evenness2, xlab="Samples")
dev.off()
Alpha2<-cbind(richness2,shannon2,evenness2)
write.csv(Alpha2, "alpha2.csv") 

aovT2<-aov(log(richness2) ~ Plant/Date/Disinfection,data=variables)
summary(aovT2)
check_model(aovT2)
testT2<-t.test(log(richnessA2),log(richnessB2), paried=T)
testT2

OUT<-capture.output(summary(aovT2))
write(OUT, "statistics.txt", append=T)
outT<-capture.output(testT2)
write(outT, "statistics.txt", append=T)

plot(richness2~variables$Origin, ylab="", xlab="")

###BETADIVERSITY:
beta2<-vegdist(args, method="bray", diag=T, upper=T)
betaA2<-vegdist(args[1:8,], method="bray", diag=T, upper=T)
betaB2<-vegdist(args[9:16,], method="bray", diag=T, upper=T)

plot(hclust(beta2, method="average"), main="", xlab="", ylab="", hang=-1, cex=0.7)
title(main="b)", adj = 0,cex.main=1.5, mgp = c(0, 0, 0))

#PERMANOVA analysis: How much of the variance in betadiversity is explained by different variables
ado2<-adonis(beta2~ variables$Plant/variables$Date/variables$Disinfection) 
ado2
adonisOUT<-capture.output(ado2)
write(adonisOUT, "statistics.txt", append=T)
adoT2<-adonis(beta2 ~ variables$Origin, permutations = 255,strata=variables$Pair)
adoT2
adoOutT<-capture.output(adoT2)
write(adoOutT, "statistics.txt", append=T)

mtl2<-mantel(betaA2,betaB2,method="spearman") #Test whether the two distance matrixes are correlated
mtl2

write(mtl2$statistic, "statistics.csv", append=T)
write(mtl2$signif, "statistics.csv", append=T)


################ HIGH-RISK ARGs
### ALPHA Diversity
hr_args<-read.csv("Ex_In_high_risk_ARGs.csv", row.names=1)
View(hr_args)
richness3<-specnumber(hr_args)  
shannon3<-diversity(hr_args)
evenness3<- shannon3/log(specnumber(hr_args))

richnessA3<-specnumber(hr_args[1:8,])
richnessB3<-specnumber(hr_args[9:16,])

par(mfrow=c(3,1)) #split the plot window in 3
plot(richness3,xlab="Samples")
plot(shannon3, xlab="Samples")
plot(evenness3, xlab="Samples")
dev.off()
Alpha3<-cbind(richness3,shannon3,evenness3)
write.csv(Alpha3, "alpha3.csv") 

aovT3<-aov(log(richness3) ~ Plant/Date/Disinfection,data=variables)
summary(aovT3)
check_model(aovT3)
testT3<-t.test(log(richnessA3),log(richnessB3), paried=T)
testT3

OUT<-capture.output(summary(aovT3))
write(OUT, "statistics.txt", append=T)
outT<-capture.output(testT3)
write(outT, "statistics.txt", append=T)

plot(richness3~variables$Origin, ylab="", xlab="")

###BETADIVERSITY:
beta3<-vegdist(hr_args, method="bray", diag=T, upper=T)
betaA3<-vegdist(hr_args[1:8,], method="bray", diag=T, upper=T)
betaB3<-vegdist(hr_args[9:16,], method="bray", diag=T, upper=T)

plot(hclust(beta3, method="average"), main="", xlab="", ylab="", hang=-1, cex=0.7)
title(main="c)", adj = 0,cex.main=1.5, mgp = c(0, 0, 0))

#PERMANOVA analysis: How much of the variance in betadiversity is explained by different variables
ado3<-adonis(beta3~ variables$Plant/variables$Date/variables$Disinfection) 
ado3
adonisOUT<-capture.output(ado3)
write(adonisOUT, "stats2.txt", append=T)
adoT3<-adonis(beta3 ~ variables$Origin, permutations = 255,strata=variables$Pair)
adoT3
adoOutT<-capture.output(adoT3)
write(adoOutT, "statistics.txt", append=T)

adoT3<-adonis(beta3 ~ variables$Origin+variables$Pair)
adoT3
adoOutT<-capture.output(adoT3)
write(adoOutT, "statistics.txt", append=T)

mtl3<-mantel(betaA3,betaB3,method="spearman") #Test whether the two distance matrixes are correlated
mtl3

write(mtl3$statistic, "statistics.csv", append=T)
write(mtl3$signif, "statistics.csv", append=T)

for (i in c(1:34)){
	model_abb <- aov(hr_args[,i] ~ Plant/Date/Disinfection,data=variables)
	print(colnames(hr_args)[i])
	print(summary(model_abb))
	out<-capture.output(cbind(colnames(hr_args)[i]),summary(model_abb))
	write(out,"statistics.txt",append=T)
}

for (i in c(1:34)){
	model_abb <- t.test(hr_args[c(1:8),i],hr_args[c(9:16),i], paried=T)
	print(colnames(hr_args)[i])
	print(model_abb)
	out<-capture.output(cbind(colnames(hr_args)[i]),model_abb)
	write(out,"statistics.txt",append=T)
}

plot(VraOTU5$BACA~variables$Plant)
plot(VraOTU5$BACA~variables$Date)
plot(VraOTU5$APH(3')-I~variables$Date)

plot(VraOTU5$ERMB~variables$Origin)
plot(VraOTU5$TOLC~variables$Origin)
plot(VraOTU5$LNUB~variables$Origin)
plot(VraOTU5$TEM~variables$Origin)
plot(VraOTU5$AADE~variables$Origin)
plot(VraOTU5$PENA~variables$Origin)


################ GRAPHS
colNew1<-c("darkblue","gray25","blue3","gray45","deepskyblue1","gray65","lightblue1","gray85","darkblue","gray25","blue3","gray45","deepskyblue1","gray65","lightblue1","gray85")

p <- ggplot(data=variables,aes(x=factor(Origin,c("Intra","Extra")), y=richness)) + 
      geom_dotplot(fill=colNew1,binaxis = "y",binwidth = 1.0,stackdir = "center",dotsize = 3.25) +
	theme(legend.position="none")+ labs (x=NULL, y="Richness (Microbial Community)")+ 
	stat_summary(fun = mean, fun.min = mean, fun.max = mean,geom = "crossbar", width = 0.5, col="black")

p2 <- ggplot(variables, aes(x=factor(Origin,c("Intra","Extra")), y=richness2))+
      geom_dotplot(fill=colNew1,
	binaxis = "y",binwidth = 1.0,stackdir = "center",dotsize = 7.75)+scale_fill_manual(values=colNew1)+
	theme(legend.position="none")+labs (x=NULL, y="Richness (Whole Antibiotic Resistome)")+
	stat_summary(fun = mean, fun.min = mean, fun.max = mean,geom = "crossbar", width = 0.5, col="black")

p3 <- ggplot(variables, aes(x=factor(Origin,c("Intra","Extra")), y=richness3)) +
      geom_dotplot(fill=colNew1,
	binaxis = "y",binwidth = 1.0,stackdir = "center",dotsize = 1.25) + scale_fill_manual(values=colNew1)+
	theme(legend.position="none")+ labs (x=NULL, y="Richness (High-risk ARGs)")+
	stat_summary(fun = mean, fun.min = mean, fun.max = mean,geom = "crossbar", width = 0.5, col="black")

library("cowplot")
plot_grid(p, p2,p3, labels=c("A","B","C"),ncol=3)

###
NMDS<-metaMDS(genera, distance="bray", k=3)
plot(NMDS)
colNew2<-c("darkblue","deepskyblue1")
data.scores <- as.data.frame(scores(NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Origin <- factor(variables$Origin,c("Intra","Extra"))  #  add the grp variable created earlier
data.scores$Disinfection <- factor(variables$Disinfection, c("Pre", "Post"))
head(data.scores)  #look at the data
d<-ggplot() + geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=Origin,fill=Origin, shape=Disinfection),size=4) + 
	scale_fill_manual(values=colNew2) +	scale_color_manual(values=colNew2)+ coord_equal() + theme_bw()+ labs (title="Microbial Community")+
	scale_shape_manual(values=c(24,22))+ylim(-1.0,0.75)+xlim(-1.25,2.5)

NMDS2<-metaMDS(args, distance="bray", k=3)
plot(NMDS2)
data.scores2 <- as.data.frame(scores(NMDS2))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores2$Origin <- factor(variables$Origin,c("Intra","Extra"))  #  add the grp variable created earlier
data.scores2$Disinfection <- factor(variables$Disinfection, c("Pre", "Post"))
head(data.scores2)  #look at the data
d2<-ggplot() + geom_point(data=data.scores2,aes(x=NMDS1,y=NMDS2,color=Origin,fill=Origin, shape=Disinfection),size=4) + 
	scale_fill_manual(values=colNew2) +	scale_color_manual(values=colNew2)+ coord_equal() + theme_bw()+ labs (title="Whole Antibiotic Resistome")+
	scale_shape_manual(values=c(24,22))+ylim(-0.85,0.65)+xlim(-1.25,2.25)

NMDS3<-metaMDS(hr_args, distance="bray", k=3)
plot(NMDS3)
data.scores3 <- as.data.frame(scores(NMDS3))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores3$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores3$Origin <- factor(variables$Origin,c("Intra","Extra"))  #  add the grp variable created earlier
data.scores3$Disinfection <- factor(variables$Disinfection, c("Pre", "Post"))
head(data.scores3)  #look at the data
d3<-ggplot() + geom_point(data=data.scores3,aes(x=NMDS1,y=NMDS2,color=Origin,fill=Origin, shape=Disinfection),size=4) + 
	scale_fill_manual(values=colNew2) +	scale_color_manual(values=colNew2)+ coord_equal() + theme_bw()+ labs (title="High-risk ARGs")+
	scale_shape_manual(values=c(24,22))+ylim(-0.7,0.7)+xlim(-1.0,2.0)

plot_grid(d,d2,d3,labels=c("A","B","C"),nrow=3)

##############################
library("reshape2")
library("RColorBrewer")

n <- 34
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector,n)

tARGs<-read.csv("Ex_In_high_risk_ARGs2.csv", row.names=1, header=T)
View(tARGs)
tClass<-read.csv("Ex_In_ARG_class.csv", row.names=1, header=T)
View(tClass)
allST1<-as.data.frame(tARGs)
datm1<-melt(cbind(allST1, ind=row.names(allST1), id.vars=c('ind')))
b<-ggplot(datm1,aes(x = variable, y = value,fill = ind)) + geom_bar(stat="identity") +scale_fill_manual(values= col)+ guides(fill=guide_legend(ncol=3))+theme(text = element_text(size=15), 
	axis.text.x = element_text(angle=90, vjust=0), legend.position="right", legend.text=element_text(size=10),legend.title=element_text(size=12)) + scale_x_discrete(name="") + 
	scale_y_continuous(name="TPM normalised abundance",breaks= c(0,75000, 150000,225000,300000,375000),limits=c(0,375000)) + theme(axis.title.y=element_text(colour="#999999"))+labs(fill="High-risk ARGs",tag ="B")
b
allST2<-as.data.frame(tClass)
datm2<-melt(cbind(allST2, ind=row.names(allST2), id.vars=c('ind')))
a<-ggplot(datm2,aes(x = variable, y = value,fill = ind)) + geom_bar(stat="identity") +scale_fill_manual(values= col)+ guides(fill=guide_legend(ncol=2))+theme(text = element_text(size=15), 
	axis.text.x = element_text(angle=90, vjust=0), legend.position="right", legend.text=element_text(size=10),legend.title=element_text(size=12)) + scale_x_discrete(name="") + 
	scale_y_continuous(name="TPM normalised abundance") + theme(axis.title.y=element_text(colour="#999999"))+labs(fill="Resistance Class",tag = "A")
a
a1<-plot_grid(a)
b1<-plot_grid(b)
plot_grid(a, b, nrow=2)