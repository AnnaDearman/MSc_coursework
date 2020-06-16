# Assessment 2
# Dataset 4: Land use change in SE Asia

library(dplyr)
library(ggplot2)
install.packages("vegan")
library(vegan)

# Bat faeces was tested for insect mit DNA to identify
# which insects (OTUs) they've been eating in three different
# types of land, categorised by the way humans have
# used the land, to see how similar the insect OTUs
# are in the bats' diets, as a proxy to indicate hunting
# strategy, which may be different in the different land
# use types

# Conduct a multivariate visualisation and analysis 
# to explore the effects of land use type on diet of bats. 

##### EXPLORATORY DATA ANALYSIS #####

# Read in raw data file
raw4 <- read.table("Experiment4_data_set 29 .txt")

str(raw4)
# 150 obs of 26 variables
# 150 individual bat faeces samples
# Treatment is a factor with 3 levels: A, B, C
# All other variables are X[num], int
# 25 insect OTUs

head(raw4)

# Check whether any OTU has no counts at all
# Need to make table without first column
# before I can use colsums
raw4justnums <- raw4[,2:length(raw4[1,])]
raw4cs <- colSums(raw4justnums)
print(raw4cs)
# No, all OTUs have counts associated with them
hist(raw4cs) # approx norm dist for each OTU

# Check whether any samples had no counts at all
raw4rs <- rowSums(raw4justnums)
print(raw4rs)
# No, all samples had OTU counts in them
hist(raw4rs) # norm dist for each sample

hist(as.matrix(raw4justnums)) # Distribution of per
# sample abundances for each OTU. Mostly 0s, negative
# skew. 

##### RUN NMDS #####
NMDS_4_1 <- metaMDS(raw4justnums,distance="bray",k=2)
# Object naming format:
# Method = NMDS
# 4 = Dataset 4
# 1 = First attempt at fitting a model

# Extract stress score
NMDS_4_1$stress
# 0.176, could be better

# Try again with higher k
NMDS_4_2 <- metaMDS(raw4justnums,distance="bray",k=3)

# Extract stress score
NMDS_4_2$stress
# 0.137, try again

# Try again with higher k
NMDS_4_3 <- metaMDS(raw4justnums,distance="bray",k=4)

# Extract stress score
NMDS_4_3$stress
# 0.114, try again

# Try again with higher k
NMDS_4_4 <- metaMDS(raw4justnums,distance="bray",k=5)

# Extract stress score
NMDS_4_4$stress
# 0.097, good

# Run stress plot
stressplot(NMDS_4_4,main="Shepard plot") # Looks good

# Run plot
plot(NMDS_4_4) # Can see three clusters, presumably 
# sites A, B, C?

# Use ordiplot to make axes

# use type="n" to suppress labels
plot1 <- ordiplot(NMDS_4_4,type='n')

# Need to colour data points by site.
# Check raw4 again
raw4$Treatment
# First 50:A, second 50:B, last 50:C

# Plot data points

# First make a selection vector based on
# the fact that we know the samples are listed
# as 50x As, then 50x Bs, then 50x Cs
selectvect <- c(rep(1,50),rep(2,50),rep(3,50))

# Then plot the points
points(NMDS_4_4,
       col=c("#167BAB","#E85B06","black")
       [selectvect])

# Sites form very distinct clusters!

# Add a legend
legend("bottomright",legend=c("Secondary rainforest","Oil palm plantations","Urban ornamental gardens"),
       col=c("#167BAB","#E85B06","black"),pch=16)

# Use 'ordihull" to draw polygons to represent the groups
ordihull(NMDS_4_4,
         groups=as.factor(c(rep("A",50),rep("B",50),rep("C",50))),
         display="sites",
         draw="polygon",
         col=c("#167BAB","#E85B06","black"),
         alpha=50)

# Use adonis function to carry out PERMANOVA

# First create distance matrix
dist.matrix.bats <- vegdist(raw4justnums,method="bray")

# Use adonis function to test for significance of between-group differences
# in composition of insect OTUs
permano.bats <- adonis(dist.matrix.bats~raw4$Treatment,permutations=500)
permano.bats # Print p value etc

# Significant adonis result: p-value = 0.002 for Treatment
# R2 is 0.46 for treatment, 0.54 for residuals
# 46% of variation in dissimilarity matrix was attributed to 
# geographical area

# Check similarity of multivariate spread between groups
betadisp.bats <- betadisper(dist.matrix.bats,raw4$Treatment)

# Are the variances different in each group?
permutest(betadisp.bats,pairwise=TRUE)

# permutest function returns a p-value of 0.001, suggesting 
# that the significant adonis result might reflect the
# variance in the individual groups (non-homogeneity of
# variance between groups)

# Assumptions of PERMANOVA were violated

# SIMPER = similarity percentage
# Use SIMPER to find the most influential OTUs (abundant
# and/or variable)

simp.bats <- simper(raw4justnums,
       group=raw4$Treatment,
       permutations=100)

print(simp.bats)

# Copy and paste printed output and paste into Excel 
# for later combination with relative abundances

# NB 70% of contribution to Bray-Curtis measures
# came from 14 species each

# NB these OTUs are not necessarily SIGNIFICANTLY different
# by virtue of appearing here, this just signifies that they 
# contribute to Bray-Curtis measures of diversity

# Find relative abundance of each species per bat faeces 
# sample

r.abund <- data.frame()

for (i in 1:150){
  row <- raw4justnums[i,]/raw4rs[i]
  r.abund <- rbind(r.abund,row)
  }

# Then take a mean for each OTU
mean.r.abund <- colSums(r.abund)/150
# This shows, on average across each bat faeces sample,
# how abundant each OTU was, relative to all OTUs
# in that bat faeces sample

# Make a vector of OTU names that contribute a combined
# 70% to Bray Curtis measures for each pairwise comparison.
# We know from before that we should take the first 14
# from the ordered lists
A_B_species <- simp.bats$A_B$ord[1:14]
A_C_species <- simp.bats$A_C$ord[1:14]
B_C_species <- simp.bats$B_C$ord[1:14]

# Make boolean vectors from these
A_B_boo <- 1:25 %in% A_B_species
A_C_boo <- 1:25 %in% A_C_species
B_C_boo <- 1:25 %in% B_C_species

# Make vectors of mean relative abundances of most 
# relevant OTUs
A_B_abunds <- mean.r.abund[A_B_boo]
A_C_abunds <- mean.r.abund[A_C_boo]
B_C_abunds <- mean.r.abund[B_C_boo]

# Put elements in correct order
A_B_abunds <- A_B_abunds[rank(simp.bats$A_B$ord[1:14])]
A_C_abunds <- A_C_abunds[rank(simp.bats$A_C$ord[1:14])]
B_C_abunds <- B_C_abunds[rank(simp.bats$B_C$ord[1:14])]

# Copy printed vectors into Excel and format
# Table can be pasted in to report
print(A_B_abunds)
print(A_C_abunds)
print(B_C_abunds)
