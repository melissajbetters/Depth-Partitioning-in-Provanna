###################################
#Title: Provanna Data Analysis
#Author: Melissa J Betters
#Finalized: August 2023
###################################
###################################
#1. Load Libraries
###################################
#If first time running this script, delete the # at the beginning of the following lines and run: 
#install.packages("car")
#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("ade4")
#install.packages("gclus")
#install.packages("ape")
#install.packages("lattice")
#install.packages("GGally")
#install.packages("gsw")
#install.packages("dplyr")
#install.packages("gridExtra")
#install.packages("RColorBrewer")
#install.packages("ggpubr")
#install.packages("ggfortify")
#install.packages("ISLR")
#install.packages("devtools")
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#Load in packages
library(car)
library(vegan)
library(ggplot2)
library(ade4)
library(gclus)
library(ape)
library(lattice)
library(GGally)
library(gsw)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(ISLR)
library(ggfortify)
library(cluster)
library(devtools)
library(pairwiseAdonis)

###################################
#2. Read in Data
###################################

#First, you will have to set the working directory to the folder that holds the datafile
#Load in data
morph_csv <- read.csv("Provanna_measurement_env_DATA.csv", header = TRUE, na.strings = "NA")
#STEP 1: Change columns to correct object types for analyses
morph_csv$Site <- as.factor(morph_csv$Site)
morph_csv$sp <- as.factor(morph_csv$sp)
morph_csv$Substrate <- as.factor(morph_csv$Substrate)
#Check
summary(morph_csv)

#STEP 2: Standardize Shell Texture according to Shell Length
rel_Texture <- (morph_csv$Texture / morph_csv$L_Shell)
as.data.frame(rel_Texture)
#Bind to dataframe
morph_csv <- cbind(morph_csv, rel_Texture)
#Check:
summary(morph_csv)

#STEP 3: Calculate roundness of Aperture
Ap_Round <- (morph_csv$W_Ap / morph_csv$L_Ap)
as.data.frame(Ap_Round)
#Bind to dataframe
morph_csv <- cbind(morph_csv, Ap_Round)
#Check
summary(morph_csv)

#STEP 4: Calculate roundness of Shell
Shell_Round <- (morph_csv$W_Shell / morph_csv$L_Shell)
as.data.frame(Shell_Round)
morph_csv <- cbind(morph_csv, Shell_Round)
summary(morph_csv)

#STEP 5: Calculate practical salinity from Conductivity, temperature, and pressure
Cond_ms_cm <- (morph_csv$Cond * 10)
PSU <- gsw_SP_from_C(Cond_ms_cm, morph_csv$TempC, morph_csv$Depth)
morph_csv <- cbind(morph_csv, PSU)
summary(morph_csv)

#STEP 6: Select relevant columns
colnames(morph_csv)
morph_set <- morph_csv[,c(3,4,5,6,7,9,12,15,16,17,19,20,21,23,24,25,26)]
colnames(morph_set)
#[1] "sp"  "L_Shell" "W_Shell" "L_Ap" "W_Ap" "Site" "Count" "TempC" "O2" "Depth"      
#[11] "Substrate"   "Thickness"   "ORP_sd_100"  "rel_Texture" "Ap_Round"    "Shell_Round" "PSU"        

#STEP 7: Summary of morphology and sites
summary(subset(morph_set, sp == 'laevis'))
summary(subset(morph_set, sp == 'goniata'))

summary(subset(morph_csv, Site == 'Jaco Scar'))
summary(subset(morph_csv, Site == 'Jaco Summit'))
summary(subset(morph_csv, Site == 'Mound 11'))
summary(subset(morph_csv, Site == 'Mound 12'))
summary(subset(morph_csv, Site == 'Quepos Seep'))
summary(subset(morph_csv, Site == 'The Thumb'))


####################################################
#3. Pearson Correlation Coefficients
####################################################
colnames(morph_set)
#Morphology
ggcorr(morph_set[,c(2,3,4,5,12,14,15,16)], label = TRUE, label_round = 3)
#Environment
ggcorr(morph_set[,c(8,9,10,13,17)], label = TRUE, label_round = 3)

################################
#4. Depth Distributions
################################

#Selecting environmental factors
colnames(morph_csv)
all_present <- morph_csv[,c(1,3,9,12,15,16,17,19,21,26)]
colnames(all_present)
#[1] "ID"   "sp" "Site"  "Count" "TempC""O2" "Depth""Substrate" "ORP_sd_100" "PSU" 
all_present <- unique(all_present)
summary(all_present)
#goniata: 18
#laevis: 22
#pacifica: 3

#Duplicate rows based on value in "Count" column:
dup_pres <- data.frame(all_present[rep(seq_len(dim(all_present)[1]), all_present$Count), 2:10, drop = FALSE], row.names=NULL)
summary(dup_pres)
#laevis: 1624
#goniata: 180
#pacifica: 9

#Drop missing values to calculate mean, below
dup_pres.na <- na.omit(dup_pres)
dup_pres.na <- droplevels.data.frame(dup_pres.na, drop = TRUE)
summary(dup_pres.na)  
#laevis: 719
#goniata: 160
#pacifica: 6

#Maximum values for each species:
aggregate(dup_pres[,c(4:6,8,9)], dup_pres[1], max)
#Mean values for each species:
aggregate(dup_pres.na[,c(4:6,8,9)], dup_pres.na[1], mean)

#Summary of values by site:
summary(subset(dup_pres, Site == "Mound 12"))
summary(subset(dup_pres, Site == "Jaco Scar"))

#REORDER FACTOR LEVELS (for graph)
dup_pres$Site <- factor(dup_pres$Site , levels=c("Jaco Scar", "Jaco Summit", "Quepos Seep",  "The Thumb", "Mound 12", "Mound 11"))
dup_pres$sp <- factor(dup_pres$sp, levels = c("laevis", "goniata", "pacifica"))
summary(dup_pres)
#The sites and species should now be listed in the order specified above.

par(mfrow = c(1,1))
ggplot(data = dup_pres, mapping = aes(x = Site, y = Depth)) +
  geom_boxplot(aes(color = sp), size = 1, outlier.size = 1) + 
  labs(y = "Depth (m)") +
  scale_color_brewer(palette = 'Set1') +
  scale_y_reverse() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =16, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 18, angle = -90, color = "black", vjust = 0.2),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))


#####################################################
#5. ANOVAs & Boxplots of species ranges of environmental variables
#####################################################
#Presence weighted by occurrence#
###########################
#Assessment of normality of environmental variables
###########################
#Create dataframe of only laevis
norml <- subset(all_present, sp == "laevis")
norml <- droplevels.data.frame(norml, drop = TRUE)
summary(norml)
#Create dataframe of only goniata
normg <- subset(all_present, sp == "goniata")
normg <- droplevels.data.frame(normg, drop = TRUE)
summary(normg)
#Create dataframe of only pacifica
normp <- subset(all_present, sp == "pacifica")
normp <- droplevels.data.frame(normp, drop = TRUE)
summary(normp)
par(mfrow = c(5,3))
#TempC
qqPlot(norml$TempC)
qqPlot(normg$TempC)
qqPlot(normp$TempC)

#O2
qqPlot(norml$O2)
qqPlot(normg$O2)
qqPlot(normp$O2)

#Depth
qqPlot(norml$Depth)
qqPlot(normg$Depth)
qqPlot(normp$Depth)

#ORP
qqPlot(norml$ORP_sd_100)
qqPlot(normg$ORP_sd_100)
qqPlot(normp$ORP_sd_100)

#PSU
qqPlot(norml$PSU)
qqPlot(normg$PSU)
qqPlot(normp$PSU)

#Specify the anova dataframe:
summary(all_present)
foranov <- subset(all_present)
foranov2 <- subset(all_present, sp != 'pacifica')
foranov2 <- droplevels.data.frame(foranov2, drop = TRUE)
summary(foranov2)
#only laevis & goniata
  
#Depth differences
depth_anov <- lm(Depth ~ sp, data = foranov)
summary(depth_anov)
## Tukey Test
depth_anov.tukey <- TukeyHSD(aov(depth_anov))$sp
# reorder tukey output by p-value
depth_anov.tukey <- depth_anov.tukey[order(depth_anov.tukey[,'p adj']),]
depth_anov.tukey
#diff        lwr       upr        p adj
#laevis-goniata   -800.1144 -856.9856 -743.2432 4.636291e-13
#pacifica-goniata -527.5556 -639.1449 -415.9664 5.624390e-13
#pacifica-laevis   272.5587  162.4282  382.6893 1.296944e-06

ggplot(data = dup_pres, mapping = aes(x = sp, y = Depth)) +
  geom_boxplot(aes(fill = sp), size = 1, outlier.size = 2) +
  labs(y = "Depth (m)") +
  scale_y_reverse() +
  annotate("text", label = "n = 1624", x = 1, y = 200, size = 4, colour = "black") +
  annotate("text", label = "n = 180", x = 2, y = 200, size = 4, colour = "black") +
  annotate("text", label = "n = 9", x = 3, y = 200, size = 4, colour = "black") +
  scale_fill_grey(start = 1, end = .7) +
  theme(
    legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#Temperature differences
temp_anov <- lm(TempC ~ sp, data = foranov)
summary(temp_anov)
## Tukey Test
temp_anov.tukey <- TukeyHSD(aov(temp_anov))$sp
# reorder tukey output by p-value
temp_anov.tukey <- temp_anov.tukey[order(temp_anov.tukey[,'p adj']),]
temp_anov.tukey
#diff       lwr        upr        p adj
#laevis-goniata    2.297893  2.0693158  2.5264695 4.636291e-13
#pacifica-goniata  1.247111  0.7986108  1.6956113 1.175674e-07
#pacifica-laevis  -1.050782 -1.4934189 -0.6081444 2.867926e-06

ggplot(data = dup_pres, mapping = aes(x = sp, y = TempC)) +
  geom_boxplot(aes(fill = sp), size = 1, outlier.size = 2) +
  labs(y = "Temperature (C)") +
  coord_cartesian(ylim = c(0, 6.5)) +
  annotate("text", label = "n = 1624", x = 1, y = 0, size = 4, colour = "black") +
  annotate("text", label = "n = 180", x = 2, y = 0, size = 4, colour = "black") +
  annotate("text", label = "n = 9", x = 3, y = 0, size = 4, colour = "black") +
  scale_fill_grey(start = 1, end = .7) +
  theme(
    legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#ORP differences
orp_anov <- lm(ORP_sd_100 ~ sp, data = foranov2)
summary(orp_anov)
## Tukey Test
orp_anov.tukey <- TukeyHSD(aov(orp_anov))$sp
# reorder tukey output by p-value
orp_anov.tukey <- orp_anov.tukey[order(orp_anov.tukey[,'p adj']),]
orp_anov.tukey
#diff        lwr        upr        p adj
#0.1428420 -0.1556344  0.4413185  0.3365948 

ggplot(data = dup_pres, mapping = aes(x = sp, y = ORP_sd_100)) +
  geom_boxplot(aes(fill = sp), size = 1, outlier.size = 2) +
  labs(y = "Oxidative Reducative Potential") +
  coord_cartesian(ylim = c(-1.5, 2.5)) +
  annotate("text", label = "n = 1624", x = 1, y = 1.8, size = 4, colour = "black") +
  annotate("text", label = "n = 180", x = 2, y = 1.8, size = 4, colour = "black") +
  annotate("text", label = "n = 9", x = 3, y = 1.8, size = 4, colour = "black") +
  scale_fill_grey(start = 1, end = .7) +
  theme(
    legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#Oxygen differences
o2_anov <- lm(O2 ~ sp, data = foranov2)
summary(o2_anov)
## Tukey Test
o2_anov.tukey <- TukeyHSD(aov(o2_anov))$sp
# reorder tukey output by p-value
o2_anov.tukey <- o2_anov.tukey[order(o2_anov.tukey[,'p adj']),]
o2_anov.tukey
#diff        lwr        upr        p adj
#-6.473270e+01 -7.053165e+01 -5.893375e+01  1.002531e-13 

ggplot(data = dup_pres, mapping = aes(x = sp, y = O2)) +
  geom_boxplot(aes(fill = sp), size = 1, outlier.size = 2) +
  labs(y = "Oxygen Concentration (uM/L)") +
  coord_cartesian(ylim = c(0, 125)) +
  annotate("text", label = "n = 1624", x = 1, y = 113, size = 4, colour = "black") +
  annotate("text", label = "n = 180", x = 2, y = 113, size = 4, colour = "black") +
  annotate("text", label = "n = 9", x = 3, y = 113, size = 4, colour = "black") +
  scale_fill_grey(start = 1, end = .7) +
  theme(
    legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#Salinity differences
sal_anov <- lm(PSU ~ sp, data = foranov)
summary(sal_anov)
## Tukey Test
sal_anov.tukey <- TukeyHSD(aov(sal_anov))$sp
# reorder tukey output by p-value
sal_anov.tukey <- sal_anov.tukey[order(sal_anov.tukey[,'p adj']),]
sal_anov.tukey
#diff          lwr         upr        p adj
#laevis-goniata   -0.05755177 -0.073865855 -4.123769e-02 3.862246e-10
#pacifica-goniata -0.03194474 -0.063955273  6.580027e-05 5.056988e-02
#pacifica-laevis   0.02560704 -0.005985042  5.719912e-02 1.321386e-01

ggplot(data = dup_pres, mapping = aes(x = sp, y = PSU)) +
  geom_boxplot(aes(fill = sp), size = 1, outlier.size = 2) +
  labs(y = "Practical Salinity Units") +
  coord_cartesian(ylim = c(34.5, 34.7)) +
  annotate("text", label = "n = 1624", x = 1, y = 34.5, size = 4, colour = "black") +
  annotate("text", label = "n = 180", x = 2, y = 34.5, size = 4, colour = "black") +
  annotate("text", label = "n = 9", x = 3, y = 34.5, size = 4, colour = "black") +
  annotate("text", label = "n = 2", x = 4, y = 34.5, size = 4, colour = "black") +
  scale_fill_grey(start = 1, end = .7) +
  theme(
    legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#######################################################
#6. PCA Plots
#######################################################
####################
#ENVIRONMENTAL VARIABLES
####################
#All species:
colnames(all_present)
#pca.e.sub <- dup_pres[,c(1,4,5,6,8,9)]
pca.e.sub <- all_present[,c(2,5,6,7,9,10)]
pca.na <- na.omit(pca.e.sub)
colnames(pca.na)
#sp    TempC             O2      Depth  ORP_sd_100             PSU      
summary(pca.na)

#Z-scale environmental variables
pcTemp.z <- scale(pca.na$TempC, center = TRUE, scale = TRUE)
pcO2.z <- scale(pca.na$O2, center = TRUE, scale = TRUE)
pcDepth.z <- scale(pca.na$Depth, center = TRUE, scale = TRUE)
pcORP.z <- scale(pca.na$ORP_sd_100, center = TRUE, scale = TRUE)
pcPSU.z <- scale(pca.na$PSU, center = TRUE, scale = TRUE)
#Compile together:
pca.e <- data.frame(pcTemp.z, pcO2.z, pcDepth.z, pcORP.z, pcPSU.z)
summary(pca.e)

#Running PCA:
pca <- rda(pca.e, scaling = 2)
pca
summary(pca)
# Extract PC scores from rda object 
pca.scores <- data.frame(pca.na, summary(pca)$sites)
pca.scores
# pull out variance explained from output for plotting on axes
var_all.expl.e <- pca$CA$eig
var_all.expl.e
pca.plot <- pca.e
pca.plot$sp <- factor(pca.na$sp, levels = c("laevis", "goniata", "pacifica"))
summary(pca.plot)

# PCA plot based on Response (morphological) variables:
pca_toplot2 <- prcomp(pca.e, scale = TRUE)
autoplot(pca_toplot2, data = pca.plot, colour = 'sp', size = 3,
         loadings = TRUE, loadings.size = 3, loadings.colour = 'black', loadings.label = ) +
  scale_color_brewer(palette = 'Set1') +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  stat_conf_ellipse(aes(color = sp), fill = NA, show.legend = TRUE, alpha = 0.3, geom = "polygon", level = 0.99, size = .5) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =16, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))


#perMANOVA of Environmental Variables
summary(pca.e)
env.dist <- vegdist(pca.e, method = "euclidean")
pairwise.adonis2(env.dist ~ sp, data = pca.plot, nperm = 10000)
#sp_vs_sp             Pr(>F)
#laevis_vs_goniata    ***
#laevis_vs_pacifica   NS
#goniata_vs_pacifica  0.05*

#######################################################
#7. ANOVAS - WITHIN SPECIES COMPARISON BY SUBSTRATE
#looking at Shell Roundness, Aperture Roundness, Size, and Relative texture
#######################################################
#Laevis
nowood <- subset(morph_set, sp == 'laevis')
nowood <- subset(nowood, Substrate != "Wood ")
nowood <- droplevels.data.frame(nowood, drop = TRUE)
summary(nowood$Substrate)
#Mussels: 65
#Rock: 5
#Tubeworm: 15

#Goniata
nochips <- subset(morph_set, sp == 'goniata')
nochips <- subset(nochips, Substrate != "Chips")
nochips$Substrate <- droplevels.factor(nochips$Substrate, drop = TRUE)
summary(nochips$Substrate) 
summary(nochips)
#Mussels: 21
#Rock: 8
#Tubeworm: 16

#Assess normality of morphological variables
colnames(nowood)
colnames(nochips)
par(mfrow = c(5,2))
#L_Ap
qqPlot(nowood$L_Ap)
qqPlot(nochips$L_Ap)
#Thickness
qqPlot(nowood$Thickness)
qqPlot(nochips$Thickness)
#Ap_Round
qqPlot(nowood$Ap_Round)
qqPlot(nochips$Ap_Round)
#Shell_Round
qqPlot(nowood$Shell_Round)
qqPlot(nochips$Shell_Round)
#rel_text
qqPlot(nochips$rel_Texture)

##################
#SHELL SIZE
##################
#LAEVIS
L_Ap_anovl <- lm(L_Ap ~ Substrate, data = nowood)
summary(L_Ap_anovl)
## Tukey Test
L_Ap_anovl.tukey <- TukeyHSD(aov(L_Ap_anovl))$Substrate
# reorder tukey output by p-value
L_Ap_anovl.tukey <- L_Ap_anovl.tukey[order(L_Ap_anovl.tukey[,'p adj']),]
L_Ap_anovl.tukey
#No significant differences

#GONIATA
L_Ap_anovg <- lm(L_Ap ~ Substrate, data = nochips)
summary(L_Ap_anovg)
## Tukey Test
L_Ap_anovg.tukey <- TukeyHSD(aov(L_Ap_anovg))$Substrate
# reorder tukey output by p-value
L_Ap_anovg.tukey <- L_Ap_anovg.tukey[order(L_Ap_anovg.tukey[,'p adj']),]
L_Ap_anovg.tukey
#Tubeworm-Mussel *

##################
#THICKNESS
##################
#LAEVIS
Thickness_anovl <- lm(Thickness ~ Substrate, data = nowood)
summary(Thickness_anovl)
## Tukey Test
Thickness_anovl.tukey <- TukeyHSD(aov(Thickness_anovl))$Substrate
# reorder tukey output by p-value
Thickness_anovl.tukey <- Thickness_anovl.tukey[order(Thickness_anovl.tukey[,'p adj']),]
Thickness_anovl.tukey
#No significant differences

#GONIATA
Thickness_anovg <- lm(Thickness ~ Substrate, data = nochips)
summary(Thickness_anovg)
## Tukey Test
Thickness_anovg.tukey <- TukeyHSD(aov(Thickness_anovg))$Substrate
# reorder tukey output by p-value
Thickness_anovg.tukey <- Thickness_anovg.tukey[order(Thickness_anovg.tukey[,'p adj']),]
Thickness_anovg.tukey
#Tubeworm-Mussel  *

##################
#SHELL ROUNDNESS
##################
#LAEVIS
shell_round_anovl <- lm(Shell_Round ~ Substrate, data = nowood)
summary(shell_round_anovl)
## Tukey Test
shell_round_anovl.tukey <- TukeyHSD(aov(shell_round_anovl))$Substrate
# reorder tukey output by p-value
shell_round_anovl.tukey <- shell_round_anovl.tukey[order(shell_round_anovl.tukey[,'p adj']),]
shell_round_anovl.tukey
#Rock-Mussel *

#GONIATA
shell_round_anovg <- lm(Shell_Round ~ Substrate, data = nochips)
summary(shell_round_anovg)
## Tukey Test
shell_round_anovg.tukey <- TukeyHSD(aov(shell_round_anovg))$Substrate
# reorder tukey output by p-value
shell_round_anovg.tukey <- shell_round_anovg.tukey[order(shell_round_anovg.tukey[,'p adj']),]
shell_round_anovg.tukey
#No significant differences

##################
#APERTURE ROUNDNESS
##################
#LAEVIS
Ap_round_anovl <- lm(Ap_Round ~ Substrate, data = nowood)
summary(Ap_round_anovl)
## Tukey Test
Ap_round_anovl.tukey <- TukeyHSD(aov(Ap_round_anovl))$Substrate
# reorder tukey output by p-value
Ap_round_anovl.tukey <- Ap_round_anovl.tukey[order(Ap_round_anovl.tukey[,'p adj']),]
Ap_round_anovl.tukey
#Tubeworm-Rock  *

#GONIATA
Ap_round_anovg <- lm(Ap_Round ~ Substrate, data = nochips)
summary(Ap_round_anovg)
## Tukey Test
Ap_round_anovg.tukey <- TukeyHSD(aov(Ap_round_anovg))$Substrate
# reorder tukey output by p-value
Ap_round_anovg.tukey <- Ap_round_anovg.tukey[order(Ap_round_anovg.tukey[,'p adj']),]
Ap_round_anovg.tukey
#No significant differences

##################
#SHELL TEXTURE
##################
#ONLY GONIATA
rel_Texture_anovg <- lm(rel_Texture ~ Substrate, data = nochips)
summary(rel_Texture_anovg)
## Tukey Test
rel_Texture_anovg.tukey <- TukeyHSD(aov(rel_Texture_anovg))$Substrate
# reorder tukey output by p-value
rel_Texture_anovg.tukey <- rel_Texture_anovg.tukey[order(rel_Texture_anovg.tukey[,'p adj']),]
rel_Texture_anovg.tukey
#No significant differences

##########################################################################
#8. ANOVAS - BETWEEN SPECIES COMPARISON BY SUBSTRATE
##########################################################################
#Automatically omits NA's

morph_all_nochips <- subset(morph_set, Substrate != "Chips" & Substrate != "Wood ")
morph_all_nochips <- subset(morph_all_nochips, sp != "pacifica")
morph_all_nochips <- droplevels.data.frame(morph_all_nochips)
#Reorder levels for graphing
morph_all_nochips$sp <- factor(morph_all_nochips$sp, levels=c("laevis", "goniata"))
summary(morph_all_nochips)
#laevis: 85
#goniata: 45

##################
#SHELL SIZE
##################

L_Ap_anov <- lm(L_Ap ~ Substrate * sp, data = morph_all_nochips)
summary(L_Ap_anov)
## Tukey Test
L_Ap_anov.tukey <- TukeyHSD(aov(L_Ap_anov))
L_Ap_anov.tukey
#Only ones I'm interested in:
#                                      diff        lwr       upr     p adj
#Mussel:goniata-Mussel:laevis      0.12979836 -0.6241331 0.8837299 0.9961392
#Rock:goniata-Rock:laevis          1.21818619 -0.3928040 2.8291763 0.2500835
#Tubeworm:goniata-Tubeworm:laevis  0.42630566 -0.5893031 1.4419144 0.8284195

ggplot(data = morph_all_nochips, mapping = aes(x = Substrate, y = L_Ap)) +
  geom_boxplot(aes(fill = sp), size = 1.25, outlier.size = 2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "Aperture Length (mm)", title = "") +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

##################
#THICKNESS
##################
Thickness_anov <- lm(Thickness ~ Substrate * sp, data = morph_all_nochips)
summary(Thickness_anov)
## Tukey Test
Thickness_anov.tukey <- TukeyHSD(aov(Thickness_anov))
Thickness_anov.tukey
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Only ones I'm interested in:
#                                       diff         lwr       upr     p adj
#Mussel:goniata-Mussel:laevis      2.896825e-02 -0.06700343 0.1249399 0.9519216
#Rock:goniata-Rock:laevis         -7.100000e-02 -0.27571459 0.1337146 0.9155767
#Tubeworm:goniata-Tubeworm:laevis  8.583333e-02 -0.04322390 0.2148906 0.3912929

ggplot(data = morph_all_nochips, mapping = aes(x = Substrate, y = Thickness)) +
  geom_boxplot(aes(fill = sp), size = 1.25, outlier.size = 2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "Shell Thickness (mm)", title = "") +
  coord_cartesian(ylim=c(0,0.75))  +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

##################
#SHELL ROUNDNESS
##################
shell_round_anov <- lm(Shell_Round ~ Substrate * sp, data = morph_all_nochips)
summary(shell_round_anov)
## Tukey Test
shell_round_anov.tukey <- TukeyHSD(aov(shell_round_anov))
shell_round_anov.tukey
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Only ones I'm interested in:
#                                          diff          lwr           upr     p adj
#Mussel:goniata-Mussel:laevis     -0.0280542004 -0.056514261  0.0004058603 0.0557752*
#Rock:goniata-Rock:laevis         -0.0987414802 -0.159657481 -0.0378254789 0.0001026***
#Tubeworm:goniata-Tubeworm:laevis -0.0446880869 -0.083091067 -0.0062851064 0.0126527*

ggplot(data = morph_all_nochips, mapping = aes(x = Substrate, y = Shell_Round)) +
  geom_boxplot(aes(fill = sp), size = 1.25, outlier.size = 2) +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = "", y = "Shell Roundness (W/L)", title = "") +
  coord_cartesian(ylim=c(0.5,0.8))  +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

##################
#APERTURE ROUNDNESS
##################
Ap_round_anov <- lm(Ap_Round ~ Substrate * sp, data = morph_all_nochips)
summary(Ap_round_anov)
## Tukey Test
Ap_round_anov.tukey <- TukeyHSD(aov(Ap_round_anov))
Ap_round_anov.tukey
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Only ones I'm interested in:
#Mussel:goniata-Mussel:laevis      0.034103545  0.002744703 0.06546239 0.0246303*
#Rock:goniata-Rock:laevis          0.042566289 -0.024440841 0.10957342 0.4444171
#Tubeworm:goniata-Tubeworm:laevis  0.010549282 -0.031693699 0.05279226 0.9787421

ggplot(data = morph_all_nochips, mapping = aes(x = Substrate, y = Ap_Round)) +
  geom_boxplot(aes(fill = sp), size = 1.25, outlier.size = 2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "Aperture Roundness (W/L)", title = "") +
  coord_cartesian(ylim=c(0.5,.93))  +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

###############################
#9. RDAs
################################
#PROVANNA LAEVIS
colnames(morph_set)
morphl <- subset(morph_set, sp == 'laevis')
summary(morphl)
#Select all variables you will need for RDA:
colnames(morphl)
morphl_select <- morphl[,c(4,8,9,10,12,13,15,16,17)]
colnames(morphl_select)
#[1] "L_Ap"        "TempC"       "O2"      "Depth"    "Thickness"   "ORP_sd_100" "Ap_Round"    "Shell_Round" "PSU"  
morphl_select.na <- na.omit(morphl_select)
summary(morphl_select.na)
#When missing values are omitted, there is only one value of PSU
#This is removed

#Explanatory variables, removing PSU as it shows no variation
rdal.e <- morphl_select.na[,c(2,3,4,6)]
summary(rdal.e)
nrow(rdal.e)
# TempC  O2 Depth ORP
#59 rows

#Transform Data
lz_TempC <- scale(rdal.e$TempC, center = TRUE, scale = TRUE)
lz_O2 <- scale(rdal.e$O2, center = TRUE, scale = TRUE)
lz_Depth <- scale(rdal.e$Depth, center = TRUE, scale = TRUE)
lz_ORP <- scale(rdal.e$ORP_sd_100, center = TRUE, scale = TRUE)
#Compile together:
rdal.e.z <- data.frame(lz_TempC, lz_O2, lz_Depth, lz_ORP)
summary(rdal.e.z)
nrow(rdal.e.z)
#59 rows

#Response Variables
colnames(morphl_select.na)
rdal.r <- morphl_select.na[,c(1,5,7,8)]
summary(rdal.r)
nrow(rdal.r)
# L_Ap         Thickness       Ap_Round       Shell_Round    
#59 rows

#Transform data
lz_L_Ap <- scale(rdal.r$L_Ap, center = TRUE, scale = TRUE)
lz_Thickness <- scale(rdal.r$Thickness, center = TRUE, scale = TRUE)
lz_Ap_Round <- scale(rdal.r$Ap_Round, center = TRUE, scale = TRUE)
lz_Shell_Round <- scale(rdal.r$Shell_Round, center = TRUE, scale = TRUE)
#Compile together:
rdal.r.z <- data.frame(lz_L_Ap, lz_Thickness, lz_Ap_Round, lz_Shell_Round)
summary(rdal.r.z)

#Run RDA
rdal <- rda(rdal.r.z ~ ., data = rdal.e.z)
summary(rdal)
rdal
#Variance:
#   Constrained: ~16.29% of the variance of Y is explained by X
#   Unconstrained: ~83.71% of the variance of Y is unexplained

#Significance test of whole model
anova.cca(rdal, step = 1000)
#Model: ***

#Significance test of explanatory variables
anova.cca(rdal, step = 1000, by = "term")
#Oxygen** 
#Depth***

#Significance test of axes
anova.cca(rdal, step = 1000, by = "axis")
#RDA1*    
#RDA2.

#Most of the variance in the response data is explained by factors other than the environment.

#PLOTTING

#Extract % explained by the first 2 axes
perc <- round(100*(summary(rdal)$cont$importance[2, 1:2]), 2)
perc

#Abby's plot - laevis
#Extract scores
sc_sil <- as.data.frame(scores(rdal, display="sites", choices=c(1,2), scaling=2))
sc_spl <- as.data.frame(scores(rdal, display="species", choices=c(1,2), scaling=2))
sc_bpl <- as.data.frame(scores(rdal, display="bp", choices=c(1, 2), scaling=2))
rownames(sc_spl) <- c('AP','TH','AP.R','SH.R')
rownames(sc_bpl) <- c('Temp.','O2','Depth','ORP')

#plot
ggplot()+
  geom_point(data=sc_sil,aes(x=RDA1,y=RDA2), color='red',size=3)+
  geom_text(data=sc_spl,aes(x=RDA1,y=RDA2,label=rownames(sc_spl)),
            color='black', size = 5)+
  geom_segment(data = sc_bpl[c('Temp.'),], aes(x = 0, y = 0, xend = RDA1,
      yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
      color = 'black', size = 1.25) +
  annotate('text',x = sc_bpl['Temp.','RDA1']-.3, 
           y = sc_bpl['Temp.','RDA2']+.2, label = 'Temp.',
           size = 5, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpl[c('O2'),], aes(x = 0, y = 0, xend = RDA1,
        yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
        color = 'black', size = 1.25) +
  annotate('text',x = sc_bpl['O2','RDA1']+.1, 
           y = sc_bpl['O2','RDA2']-.1, label = 'O2**',
           size = 6, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpl[c('Depth'),], aes(x = 0, y = 0, xend = RDA1,
          yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
          color = 'black', size = 1.25) +
  annotate('text',x = sc_bpl['Depth','RDA1']-.1, 
           y = sc_bpl['Depth','RDA2']-.1, label = 'Depth***',
           size = 6, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpl[c('ORP'),], aes(x = 0, y = 0, xend = RDA1,
               yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
               color = 'black', size = 1.25) +
  annotate('text',x = sc_bpl['ORP','RDA1']+.3, 
           y = sc_bpl['ORP','RDA2']+.1, label = 'ORP',
           size = 5, color='black',fontface = 'bold')+
  labs(x='RDA1 (9.67%)*',
       y='RDA2 (5.37%)')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = .5),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(color = "white", size = 1, ),
    panel.grid.minor.y = element_line(color = "white", size = 1.2),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", size = 2, fill = NA))


#PROVANNA GONIATA

colnames(morph_set)
morphg <- subset(morph_set, sp == 'goniata')
summary(morphg)
#Select all variables you will use in the RDA:
morphg_select <- morphg[,c(4,8,9,10,12,13,14,15,16,17)]
colnames(morphg_select)
#[1] "L_Ap"        "TempC"       "O2"    "Depth"      "Thickness"   "ORP_sd_100"  "rel_Texture" "Ap_Round"    "Shell_Round" "PSU"       
morphg_select.na <- na.omit(morphg_select)
summary(morphg_select.na)

#Select explanatory variables
rdag.e <- morphg_select.na[,c(2,3,4,6,10)]
summary(rdag.e)
nrow(rdag.e)
# TempC             O2             Depth          ORP_sd_100             PSU   
#55 rows

#Transform Data
gz_TempC <- scale(rdag.e$TempC, center = TRUE, scale = TRUE)
gz_O2 <- scale(rdag.e$O2, center = TRUE, scale = TRUE)
gz_Depth <- scale(rdag.e$Depth, center = TRUE, scale = TRUE)
gz_ORP <- scale(rdag.e$ORP_sd_100, center = TRUE, scale = TRUE)
gz_PSU <- scale(rdag.e$PSU, center = TRUE, scale = TRUE)
#Compile together:
rdag.e.z <- data.frame(gz_TempC, gz_O2, gz_Depth, gz_ORP, gz_PSU)
summary(rdag.e.z)
nrow(rdag.e.z)
#55 rows

#RESPONSE Variables
colnames(morphg_select.na)
rdag.r <- morphg_select.na[,c(1,5,7,8,9)]
summary(rdag.r)
nrow(rdag.r)
# L_Ap         Thickness      rel_Texture    Ap_Round       Shell_Round    
#55 rows

#Create dataset
gz_L_Ap <- scale(rdag.r$L_Ap, center = TRUE, scale = TRUE)
gz_Thickness <- scale(rdag.r$Thickness, center = TRUE, scale = TRUE)
gz_Texture <- scale(rdag.r$rel_Texture, center = TRUE, scale = TRUE)
gz_Ap_Round <- scale(rdag.r$Ap_Round, center = TRUE, scale = TRUE)
gz_Shell_Round <- scale(rdag.r$Shell_Round, center = TRUE, scale = TRUE)
#Compile together:
rdag.r.z <- data.frame(gz_L_Ap, gz_Thickness, gz_Texture, gz_Ap_Round, gz_Shell_Round)
summary(rdag.r.z)

#Run RDA
rdag <- rda(rdag.r.z ~ ., data = rdag.e.z)
summary(rdag)
#Constrained: ~18% of the variance of Y is explained by X
#Unconstrained: ~81.94% of the variance of Y is unexplained

#Significance test of whole model
anova.cca(rdag, step = 1000)
#Model: **

#Significance test of explanatory variables
anova.cca(rdag, step = 1000, by = "term")
#Oxygen*, ORP*                      

#Significance test of axes
anova.cca(rdag, step = 1000, by = "axis")
#RDA1*            

#Most of the variance in the response data is explained by factors other than the environment.

#Abby's plot - goniata
#Extract % explained by the first 2 axes
perc <- round(100*(summary(rdag)$cont$importance[2, 1:2]), 2)
perc
#Extract scores
sc_sig <- as.data.frame(scores(rdag, display="sites", choices=c(1,2), scaling=2))
sc_spg <- as.data.frame(scores(rdag, display="species", choices=c(1,2), scaling=2))
sc_bpg <- as.data.frame(scores(rdag, display="bp", choices=c(1, 2), scaling=2))
rownames(sc_spg) <- c('AP','TH','TX','AP.R','SH.R')
rownames(sc_bpg) <- c('Temp.','O2','Depth','ORP','PSU')

#plot
ggplot()+
  geom_point(data=sc_sig,aes(x=RDA1,y=RDA2), color='cornflowerblue',size=3)+
  geom_text(data=sc_spg,aes(x=RDA1,y=RDA2,label=rownames(sc_spg)),
            color='black', size = 5)+
  geom_segment(data = sc_bpg[c('Temp.'),], aes(x = 0, y = 0, xend = RDA1,
              yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
              color = 'black', size = 1.25) +
  annotate('text',x = sc_bpg['Temp.','RDA1']-.2, 
           y = sc_bpg['Temp.','RDA2']+.2, label = 'Temp',
           size = 5, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpg[c('O2'),], aes(x = 0, y = 0, xend = RDA1,
              yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
               color = 'black', size = 1.25) +
  annotate('text',x = sc_bpg['O2','RDA1']+.2, 
           y = sc_bpg['O2','RDA2']+.2, label = 'O2*',
           size = 6, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpg[c('ORP'),], aes(x = 0, y = 0, xend = RDA1,
               yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
               color = 'black', size = 1.25) +
  annotate('text',x = sc_bpg['ORP','RDA1']+.2, 
           y = sc_bpg['ORP','RDA2']+.2, label = 'ORP*',
           size = 6, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpg[c('PSU'),], aes(x = 0, y = 0, xend = RDA1,
               yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
               color = 'black', size = 1.25) +
  annotate('text',x = sc_bpg['PSU','RDA1']-.2, 
           y = sc_bpg['PSU','RDA2'], label = 'PSU',
           size = 5, color='black',fontface = 'bold')+
  geom_segment(data = sc_bpg[c('Depth'),], aes(x = 0, y = 0, xend = RDA1,
              yend = RDA2), arrow = arrow(length = unit(1/2, 'picas')),
               color = 'black', size = 1.25) +
  annotate('text',x = sc_bpg['Depth','RDA1']+.1, 
           y = sc_bpg['PSU','RDA2']+.3, label = 'Depth',
           size = 5, color='black',fontface = 'bold')+
  labs(x='RDA1 (9.56%)*',
       y='RDA2 (6.09%)')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = .5),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(color = "white", size = 1, ),
    panel.grid.minor.y = element_line(color = "white", size = 1.2),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", size = 2, fill = NA))

#####################
#Supp. Linear relationships
#####################
##########
#LAEVIS
##########
shr.l.o2 <- ggplot(data = morphl, mapping = aes(x = O2, y = Shell_Round)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Shell Roundness (W/L)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

th.l.o2 <- ggplot(data = morphl, mapping = aes(x = O2, y = Thickness)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Shell Thickness (mm)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

ap.l.o2 <- ggplot(data = morphl, mapping = aes(x = O2, y = L_Ap)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Aperture Length (mm)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

apr.l.o2 <- ggplot(data = morphl, mapping = aes(x = O2, y = Ap_Round)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "O2 (uM/L)", y = "Aperture Roundness (W/L)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

shr.l.depth <- ggplot(data = morphl, mapping = aes(x = Depth, y = Shell_Round)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Shell Roundness (W/L)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

th.l.depth <- ggplot(data = morphl, mapping = aes(x = Depth, y = Thickness)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Shell Thickness (mm)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

ap.l.depth <- ggplot(data = morphl, mapping = aes(x = Depth, y = L_Ap)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Aperture Length (mm)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

apr.l.depth <- ggplot(data = morphl, mapping = aes(x = Depth, y = Ap_Round)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "Depth (m)", y = "Aperture Roundness (W/L)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

grid.arrange(shr.l.o2, shr.l.depth,
             th.l.o2, th.l.depth,
             ap.l.o2, ap.l.depth,
             apr.l.o2, apr.l.depth, 
             nrow = 4, ncol = 2)

##########
#GONIATA
##########
shr.g.o2 <- ggplot(data = morphg, mapping = aes(x = O2, y = Shell_Round)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Shell Roundness (W/L)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

shr.g.orp <- ggplot(data = morphg, mapping = aes(x = ORP_sd_100, y = Shell_Round)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

th.g.o2 <- ggplot(data = morphg, mapping = aes(x = O2, y = Thickness)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Shell Thickness (mm)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))


th.g.orp <- ggplot(data = morphg, mapping = aes(x = ORP_sd_100, y = Thickness)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

ap.g.o2 <- ggplot(data = morphg, mapping = aes(x = O2, y = L_Ap)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Aperture Length (mm)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

ap.g.orp <- ggplot(data = morphg, mapping = aes(x = ORP_sd_100, y = L_Ap)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

apr.g.o2 <- ggplot(data = morphg, mapping = aes(x = O2, y = Ap_Round)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "Aperture Roundness (W/L)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

apr.g.orp <- ggplot(data = morphg, mapping = aes(x = ORP_sd_100, y = Ap_Round)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "", y = "") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

tx.g.o2 <- ggplot(data = morphg, mapping = aes(x = O2, y = rel_Texture)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "O2 (uM/L)", y = "Relative Shell Texture") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

tx.g.orp <- ggplot(data = morphg, mapping = aes(x = ORP_sd_100, y = rel_Texture)) +
  geom_point(color = "cornflowerblue", size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(x = "ORP", y = "") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =10, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

grid.arrange(shr.g.o2, shr.g.orp, 
             th.g.o2, th.g.orp, 
             ap.g.o2, ap.g.orp, 
             apr.g.o2, apr.g.orp,
             tx.g.o2, tx.g.orp,
             nrow = 5, ncol = 2)

#################
#END CODE
#################