#The Role of microbes in the Silwood Park plant-pollinator network
#Sofia Riccomagno, CID: 01092789

#Clear environment
rm(list = ls())

#Load libraries
library(dplyr)
library(RCurl)
library(dplyr)
library(stringr)
library(bipartite)
library(viridis)
library(reshape2)
library(ggplot2)

#Load datasets
obs_data <- read.csv("https://raw.githubusercontent.com/sofiariccomagno/plant-pollinator-microbes/main/observations.csv")
microbial_data <- read.csv("https://raw.githubusercontent.com/sofiariccomagno/plant-pollinator-microbes/main/microbial_counts.csv")
morphotype_data <- read.csv("https://raw.githubusercontent.com/sofiariccomagno/plant-pollinator-microbes/main/morphology.csv",
                            row.names = 1, check.names = F)

###############################Hypothesis 1###############################

##Average number of colonies per sample
#Subset to only include bee data
bee_data <- subset(microbial_data, is_bee == 1)

#Assign new column for genus
bee_data$genus<-word(bee_data$species, 1)

#Assign new column for sociality
bee_data$sociality<-ifelse(bee_data$genus=="Apis" | bee_data$genus=="Bombus","Eusocial","Solitary")

#Make summary table
no_colonies_bee <- bee_data %>% group_by(genus) %>% summarise(Samples = n(),
                                                         Samples_with_growth=sum(no_colonies > 0, na.rm = TRUE),
                                                         Total_colonies=sum(no_colonies),
                                                         Mean_types_colonies = mean(no_types_colonies),
                                                         Mean_colonies_per_sample = mean(no_colonies),
                                                         Max_colonies=max(no_colonies))
print(no_colonies_bee)

#Subset by sociality
eusocial <- subset(bee_data, sociality == "Eusocial")

#Subset by genus
apis <- subset(bee_data, genus == "Apis")
bombus <- subset(bee_data, genus == "Bombus")

#Check for normality - p<0.05, data not normally distributed -> Mann-Whitney U test
shapiro.test(eusocial$no_colonies)

#Make new vectors for test
apis_colonies <- apis$no_colonies
bombus_colonies <- bombus$no_colonies

#Two sample Kolmogorov-Smirnov test to check distributions are at least similar shapes
ks.test(apis$no_colonies, bombus$no_colonies)

#Mann-Whitney U test
wilcox.test(apis_colonies, bombus_colonies, paired = F, alternative = "greater")
#W = 300, p-value = 0.1472


##Build association network

#Replace NAs with zeros
morphotype_data[is.na(morphotype_data)] <- 0

#Transpose table and merge morphotype and microbial observations data sets together
morphotype_data<-t(morphotype_data)
h1_df<- merge(x = morphotype_data,y = bee_data, by.x=0, by.y='sample', all.y=T)
h1_df[is.na(h1_df)] <- 0
h1_df <- melt(h1_df, id.vars = c("Row.names", "is_bee","no_colonies",'no_types_colonies','species'))

#Subset only to include data where frequency of association > 0
h1_df <- subset(h1_df, value > 0)

#Prepare matrix for network
in.bees <- h1_df$species 
in.microbes <- h1_df$variable 
df_for_matrix_h1 <- data.frame(in.bees, in.microbes)
bee_matrix<- table(df_for_matrix_h1)
bee_matrix<- bee_matrix[,c(1:12)]


#Prepare colour palette so that the sum of the number of colonies found for a link will be 
#proportionally darker as the number increases
x_h1 <- ncol(bee_matrix) * nrow(bee_matrix)
colour_palette_h1 <- rainbow(x_h1)
positions_controlling_links_h1 <- c(1,2,4,5,6,8,9,12)
intensity_colour_h1 <- c(4,1,1,1,1,2,3,2)
colour_vector_h1 <- viridis_pal(alpha = 0.75,direction = -1, option = 'B')(4)
colours_h1 <- colour_vector_h1[intensity_colour_h1]
colour_palette_h1[positions_controlling_links_h1] <- colours_h1

#Plot association network
plotweb(bee_matrix,method = "cca",
        abuns.type='additional',
        arrow="up.center",
        text.rot=90,
        col.interaction = colour_palette_h1,
        bor.col.interaction = NA, 
        high.lab.dis = 0.05,
        ybig=1.2,
        y.width.high = .06,
        y.lim = c(-0.5,2),
)

#Gather statistics
networklevel(bee_matrix, index = "H2") # H2 0.3668402
specieslevel(bee_matrix, index = "d") #Apis mellifera d' = 0.3115739 #Bombus pascuorum d' = 0.4796249


#Null models

#Set seed to have repeatable results
set.seed(2022)

#Make vector of observed indices
Iobs_bee <- specieslevel(bee_matrix,index = "d")
Iobs_bee <- Iobs_bee$`lower level`[1]
Iobs_apis <- Iobs_bee$d[1]
Iobs_bombus <- Iobs_bee$d[2]

#Run null model of network 10000 times
nulls_bee <- nullmodel(web=bee_matrix, N=10000, method='r2d') 

#Make vector of indices from null models (takes a bit of time!)
Inulls_bee <-sapply(nulls_bee, function(x) specieslevel(x,index = "d"))

#Create empty vectors to be filled
apis_d <- NULL
bombus_d <- NULL

#Populate with vector of indices 
for (i in 1:length(nulls_bee)){
  d_score <- Inulls_bee[,i]$`lower level`
  apis_prov_d <- d_score$d[1]
  apis_d <- append(apis_d, apis_prov_d)
  bombus_prov_d <- d_score$d[2]
  bombus_d <- append(bombus_d, bombus_prov_d)
  i <- i+1
}

#Calculate z-scores
print(zscore_apis<-(Iobs_apis-mean(apis_d))/sd(apis_d)) #-0.0905120
print(zscore_bombus<-(Iobs_bombus-mean(bombus_d))/sd(bombus_d)) #0.52707

#More plots

#Interaction strength matrix plot
visweb(bee_matrix, type="nested", labsize = 0.3)

#Modules plot
bee_matrix <- apply(as.matrix.noquote(bee_matrix),2,as.numeric)
row.names(bee_matrix) <-c("Apis mellifera", "Bombus pascuorum")

modules_bee <- computeModules(bee_matrix)
plotModuleWeb(modules_bee,labsize=0.8) 


###############################Hypothesis 2###############################

##Average number of colonies per sample

#Subset to only include flower data
flower_data <- subset(microbial_data, is_bee == 0)

#Make summary table
no_colonies_flower <- flower_data %>% group_by(control) %>% summarise(Samples = n(),
                                                                 Samples_with_growth=sum(no_colonies > 0, na.rm = TRUE),
                                                                 Total_colonies=sum(no_colonies),
                                                                 Mean_types_colonies = mean(no_types_colonies),
                                                                 Mean_colonies_per_sample = mean(no_colonies),
                                                                 Max_colonies=max(no_colonies))
print(no_colonies_flower)

#Add column for genus
flower_data$genus <- word(flower_data$species, 1)

#Subset to only include species that have samples for both treatments
paired_flower <- subset(flower_data, genus == "Cirsium" | genus == "Glechoma" | genus == "Jacobaea" | genus == "Lotus" | genus == "Pentaglottis" | genus == "Veronica" | genus == "Vicia")

#Summary table for paired samples
no_colonies_paired <- paired_flower %>% group_by(genus, control) %>% summarise(Samples = n(),
                                                                               Samples_with_growth=sum(no_colonies > 0, na.rm = TRUE),
                                                                               Total_colonies=sum(no_colonies),
                                                                               Mean_types_colonies = mean(no_types_colonies),
                                                                               Mean_colonies_per_sample = mean(no_colonies),
                                                                               Max_colonies=max(no_colonies))
print(no_colonies_paired)

#Subset by test and control treatments
test <- subset(no_colonies_paired, control == 0)
control<- subset(no_colonies_paired, control == 1)

#Check for normality - p<0.05, data not normally distributed -> Mann-Whitney U test
shapiro.test(no_colonies_paired$Mean_colonies_per_sample)

#Two sample Kolmogorov-Smirnov test to check distributions are at least similar shapes
ks.test(control$Mean_colonies_per_sample, test$Mean_colonies_per_sample)

#Make new vector for test
control_colonies <- control$Mean_colonies_per_sample
test_colonies <- test$Mean_colonies_per_sample
#Mann-Whitney U test
wilcox.test(control_colonies, test_colonies, paired = T) #V = 10, p-value = 0.5896


##Association networks

#Merge morphotype and microbial observations data sets together
h2_df<- merge(x = morphotype_data,y = flower_data, by.x=0, by.y='sample', all.y=T)
h2_df[is.na(h2_df)] <- 0
h2_df <- melt(h2_df, id.vars = c("Row.names", "is_bee","no_colonies",'no_types_colonies','species'))

#Subset only to include data where frequency of association > 0
h2_df <- subset(h2_df, value > 0)

#Test treatment association network
#Subset to just include test treatment data
test <- subset(h2_df, control == 0)

#Prepare matrix for network
test_in.flowers <- test$species 
test_in.microbes <- test$variable 
test_df <- data.frame(test_in.flowers, test_in.microbes)
test_matrix<- table(test_df)

#Prepare colour palette so that the sum of the number of colonies found for a link will be 
#proportionally darker as the number increases
x <- ncol(test_matrix) * nrow(test_matrix)
colP <- rainbow(x)
controlling_positions <- c(7,8,9,10,11,12,15,18,20,22,24,27,28,29,32,33,35,38,40,44,46)
easy_colours <- c(1,4,8,4,5,7,5,1,1,2,3,1,1,9,2,2,1,1,2,6,1)
color <- viridis_pal(alpha = 0.75,direction = -1, option = 'B')(9)
final_colours <- color[easy_colours]
colP[controlling_positions] <- final_colours

plotweb(test_matrix,method = "cca",
        abuns.type='additional',
        arrow="up.center",
        text.rot=90,
        col.interaction = colP,
        bor.col.interaction = NA, #remove the black border color
        high.lab.dis = 0.05, #top label distance
        ybig=1.2,
        y.width.high = .06,
        y.lim = c(-0.5,2),
)


networklevel(test_matrix, index = "H2") #H2 0.2941544

#Control treatment
#control
control <- subset(h2_df, control == 1)
control_in.flowers <- control$species 
control_in.microbes <- control$variable 
control_df <- data.frame(control_in.flowers, control_in.microbes)
control_matrix<- table(control_df)

x <- ncol(control_matrix) * nrow(control_matrix)
colP <- rainbow(x)
controlling_positions <- c(1,2,4,6,8,9,13,27,30,33,40,42,54,55,57,61,66)
easy_colours <- c(9,6,3,6,3,5,1,1,1,1,3,1,4,2,1,8,7)
color <- viridis_pal(alpha = 0.75,direction = -1, option = 'B')(9)
final_colours <- color[easy_colours]
colP[controlling_positions] <- final_colours

plotweb(control_matrix,method = "cca",
        abuns.type='additional',
        arrow="up.center",
        text.rot=90,
        #col.high=ifelse(colnames(pol.matrix)=="Apis mellifera","red","grey80"),
        col.interaction = colP,
        #col.low = ifelse(colnames(pol.matrix)=="Apis mellifera","red","grey80"),
        bor.col.interaction = NA, #remove the black border color
        high.lab.dis = 0.05, #top label distance
        ybig=1.2,
        y.width.high = .06,
        y.lim = c(-0.5,2),
)


networklevel(control_matrix, index = "H2") #H2 0.3254601 

#Null models
set.seed(2022)
#Control
Iobs_control <- networklevel(control_matrix,index = "H2")
nulls_control <- nullmodel(web=control_matrix, N=10000, method='r2d') 
Inulls_control<-sapply(nulls_control, function(x) networklevel(x,index = "H2"))
zscore_control<-(Iobs_control-mean(Inulls_control))/sd(Inulls_control)
print(zscore_control) #0.7450789

#Test
Iobs_test <- networklevel(test_matrix,index = "H2")
nulls_test <- nullmodel(web=test_matrix, N=10000, method='r2d')
Inulls_test<-sapply(nulls_test, function(x) networklevel(x,index = "H2"))
zscore_test<-(Iobs_test-mean(Inulls_test))/sd(Inulls_test)
print(zscore_test)



###############################Hypothesis 3###############################
#prepare matrix
observations <- obs_data

in.flowers <- observations$flower_species 
in.pollinate <- observations$pollinator_species 



df <- data.frame(in.flowers, in.pollinate) 
pol.matrix<- table(df) 



morph_df <- morphotype_data


inter_df<-microbial_data

df<-merge(x=morph_df,y=inter_df,by.x=0,by.y='sample',all.x=T)


gdf<-df%>% group_by(species) %>% summarise(
  MorphoMorphotypeA=sum(`Morphotype A`),
  MorphoMorphotypeB=sum(`Morphotype B`),
  MorphoMorphotypeC=sum(`Morphotype C`),
  MorphoMorphotypeD=sum(`Morphotype D`),
  MorphoMorphotypeE=sum(`Morphotype E`),
  MorphoMorphotypeF=sum(`Morphotype F`),
  MorphoMorphotypeG=sum(`Morphotype G`),
  MorphoMorphotypeH=sum(`Morphotype H`),
  MorphoMorphotypeI=sum(`Morphotype I`),
  MorphoMorphotypeJ=sum(`Morphotype J`),
  MorphoMorphotypeK=sum(`Morphotype K`),
  MorphoMorphotypeL=sum(`Morphotype L`)) #group to get sum of associations by type

gdf$microbe_Morphotypes=rowSums(gdf > 0)-1 # number Morphotype of microbes associated with

microbe_freq<-select(gdf,c('species','microbe_Morphotypes'))


maxColorValue <- 10
palette <- colorRampPalette(c("grey80","red"))(7)

microbe_freq$colour<-palette[microbe_freq$microbe_MorphoMorphotypes] #set colours

plotweb(pol.matrix,method = "cca", #plot web
        text.rot=90,
        col.high=c("#EE4444","grey80","#DD8888","grey80","grey","grey80","grey80","grey80"),
        col.low =c("#FF0000","grey80","#CCCCCC","#E56666","#EE4444","#DD8888","#E56666","#DD8888","grey80",'grey80'),
        bor.col.interaction = NA, #remove the black border color
        high.lab.dis = 0.05, #top label distance
        ybig=1.2,
        y.width.high = .06,
        y.lim = c(-0.7,3),
        labsize = 1.05
)



#### Calculate network statistics
networklevel(pol.matrix)

sl<-specieslevel(pol.matrix)
sl_bee_pp<-sl[[1]]
sl_plant_pp<-sl[[2]]


##### produce regression for bees
df <- morphotype_data


df_bee <- microbial_data
bee <- subset(df_bee, is_bee == 1) #get bee data

df<-t(df)
new_df <- merge(x=df,y=bee,by.x=0,by.y='sample',all.y=T)
new_df[is.na(new_df)] <- 0 #df with bees and microbial

melt_df <- melt(new_df, id.vars=c("Row.names", "is_bee","no_colonies",'no_types_colonies','species'))
melt_df <- subset(melt_df, value > 0) #reshape df



in.bees <- melt_df$species 
in.microbes <- melt_df$variable 
df <- data.frame(in.bees, in.microbes)
bee_matrix<- table(df)

sl_bee<-specieslevel(bee_matrix)[[2]]%>%select(c(d)) #get d'


df_for_big<-bee %>% group_by(species) %>% summarise(Samples=n())
df_for_big<-merge(df_for_big,sl_bee,by.x="species",by.y=0,all.x=T)
df_for_big<-merge(df_for_big,select(sl_bee_pp,c(d)),by.x="species",by.y=0,all.x=T)
colnames(df_for_big)<-c("species","samples","degree_microbes","degree_pp")

lmbee = lm(degree_microbes~degree_pp, data = df_for_big) #Create the linear regression
summary(lmbee) #Review the results #linear regression

interaction_graph<-ggplot(df_for_big, aes(x=degree_pp, y=degree_microbes))+
  geom_point()+
  xlab("Specialisation d' bee-flower") +
  ylab("Specialisation d' bee-microbe")

interaction_graph <- interaction_graph + geom_smooth(method="lm", col="black")

interaction_graph #plot

##### produce regression for flowers, same as bees
df_flower <- read.csv("microbial_tidy.csv")
flower <- subset(df_flower, is_bee == 0)


df <- read.csv("morphology.csv", row.names = 1, check.names = F)
df<-t(df)
new_df <- merge(x=df,y=flower,by.x=0,by.y='sample',all.y=T)

melt_df <- melt(new_df, id.vars=c("Row.names", "is_bee","no_colonies",'no_types_colonies','species'))
melt_df <- subset(melt_df, value > 0)



in.flowers <- melt_df$species 
in.microbes <- melt_df$variable 
df <- data.frame(in.flowers, in.microbes)
flower_matrix<- table(df)

sl_flower<-specieslevel(flower_matrix)[[2]]%>%select(c(d))


df_for_big_flower<-flower %>% group_by(species) %>% summarise(Samples=n())
df_for_big_flower<-merge(df_for_big_flower,sl_flower,by.x="species",by.y=0,all.x=T)
df_for_big_flower<-merge(df_for_big_flower,select(sl_plant_pp,c(d)),by.x="species",by.y=0,all.x=T)
colnames(df_for_big_flower)<-c("species","samples","degree_microbes","degree_pp")

lmflower = lm(degree_microbes~degree_pp, data = df_for_big_flower) #Create the linear regression
summary(lmflower) #Review the results, linear model


interaction_graph<-ggplot(df_for_big_flower, aes(x=degree_pp, y=degree_microbes))+
  geom_point()+
  xlab("Specialisation d' flower-bee") +
  ylab("Specialisation d' flower-microbe")

interaction_graph <- interaction_graph + geom_smooth(method="lm", col="black")

interaction_graph






































































