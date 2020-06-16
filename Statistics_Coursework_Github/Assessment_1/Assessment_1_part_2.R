# Assessment 1: Probiotics and synbiotics
# Investigation 2

##### EXPLORATORY DATA ANALYSIS #####

# Read in raw data file
raw2 <- read.table("Experiment2_data_set 29 .txt")

# Check structure
str(raw2) # As expected, 120 obs of 3 variables, coded
# as Probiotic: Factor w/ 2 levels "No", "Yes",
# Diet: Factor w/ 2 levels "Highfibre", "Normal"
# and "Shannon.diversity": num. Further checking
# to follow:

# Describe variables in relation to present study,
# and sanity-check raw data:

# RESPONSE VARIABLE, continuous (num):
# Participants' microbial diversity was measured by
# calculating the Shannon diversity index in six stool 
# samples and taking a mean.
# View all data points and check length
raw2$Shannon.diversity # All numbers recorded to 2dp
# One negative value! This shouldn't be possible.

# I will remove the observation with a negative Shannon
# diversity score as I cannot infer a correct value
install.packages("dplyr")
library(dplyr)

raw2clean <- raw2 %>%
  filter(Shannon.diversity>0)

length(raw2clean$Shannon.diversity) # 119 as expected
raw2clean$Shannon.diversity # All numbers recorded to 2dp
range(raw2clean$Shannon.diversity) # No negative values
# Range 0.93-7.79

# Have a quick look at the distribution of Shannon
# diversity means:
hist(raw2clean$Shannon.diversity,breaks=50) # Almost a 
# positively skewed norm dist but an additiona, medium 
# sized peak at around 6-7

# EXPLANATORY VARIABLE 1, factor, 2 levels:
# Participants were either given a probiotic or not.
# This is indicated in the raw data file as factor
# "Probiotic" which has two levels "No" and "Yes".
# View all data points, levels and check length
raw2clean$Probiotic # Two levels, looks good
length(raw2clean$Probiotic) # 119 as expected after
# removal. Check whether 59 and 60 obs for each level:
plot(raw2clean$Probiotic) # Yes, 59 and 60 

# EXPLANATORY VARIABLE 2, factor, 2 levels:
# Participants were either on a high fibre diet or not
# This is indicated in the raw data file as factor
# "Diet" which has two levels "Highfibre" and "Normal".
# View all data points, levels and check length
raw2clean$Diet # Two levels, looks good
length(raw2clean$Diet) # 119 as expected
# Check whether 59 and 60 obs for each level:
plot(raw2clean$Diet) # Yes, 59 and 60

# Compare response variable to each explanatory variable:
boxplot(raw2clean$Shannon.diversity~raw2clean$Probiotic)
# Two slightly skewed normal distributions. 
# Negative skew for ProbioticNo, positive skew for
# ProbioticYes. There was a higher Shannon diversity for 
# participants who took the probiotic

boxplot(raw2clean$Shannon.diversity~raw2clean$Diet)
# Two normal distributions. There was a higher Shannon
# diversity for participants with the high fibre diet

# Check that ProbioticYes and ProbioticNo were evenly
# distributed between HighfibreYes and HighfibreNO

str(raw2clean %>%
      filter(Diet=="Highfibre",Probiotic=="Yes"))  
# 30 obs

str(raw2clean %>%
      filter(Diet=="Normal",Probiotic=="Yes"))  
# 30 obs

str(raw2clean %>%
      filter(Diet=="Highfibre",Probiotic=="No"))  
# 29 obs, not 30, due to removal of erroneous obs

str(raw2clean %>%
      filter(Diet=="Normal",Probiotic=="No"))  
# 30 obs

# Boxplot for each sub-group
boxplot(raw2clean$Shannon.diversity~raw2clean$Probiotic*
          raw2clean$Diet)
# All roughly norm dist. ProbioticYes gives much higher 
# Shannon diversity than ProbioticNo. DietHighfibre 
# gives noticeably higher Shannon diversity than
# DietNormal.
# Possible evidence of synergy between the two.

##### CHOICE OF MODEL #####

# As the data appears to roughly follow a normal
# distribution, the response variable is continuous,
# and there are two explanatory variables, both of
# which are factors, use two-way ANOVA

# Fit model. Include interaction term in first attempt.
# Use lm() function
ANOVA_2_1 <- lm(raw2clean$Shannon.diversity~
                  raw2clean$Probiotic*raw2clean$Diet)
# Also use aov() function to allow post hoc Tukey test
ANOVA_2_2 <- aov(raw2clean$Shannon.diversity~
                   raw2clean$Probiotic*raw2clean$Diet)
# Object naming format:
# Model type ANOVA
# 2 = Dataset 2
# 1 = First attempt at fitting a model

# Look at diagnostic plots:
plot(ANOVA_2_1) # Residuals vs Fitted: variance looks equal
# QQ plot looks good
# No influential outliers according to Cook's distance

# Look at distribution of errors as histogram
hist(ANOVA_2_1$residuals) # Fairly normal

# Check whether interaction term is significant

# Look at summary table from lm object to:
# Check whether interaction terms is significant
# get group means 
# get F statistic:
summary(ANOVA_2_1)

# Interaction term is significant
# p-value is <0.001
# Minimum adequate model has been reached
# Retain main effects terms in model, but don't report 

# "Intercept" must be "DietHighfibre:ProbioticNo". 
# "DietHighfibre:ProbioticNo" Group mean = 2.8069
# "DietHighfibre:ProbioticYes": Group mean = 6.2757 = 2.8069+3.4688
# "DietNormal:ProbioticNo": Group mean = 2.1047 = 2.8069+-0.7022
# "DietNormal:ProbioticYes": Group mean = 4.4391 = 2.8069+-1.1344+3.4688+-0.7022

# Double check means using dplyr:

raw2means <- raw2clean %>%
  filter(Diet=="Highfibre",Probiotic=="No")%>%
  summarize(DietHi_ProNo=mean(Shannon.diversity))
# 2.807 - Correct

raw2means[2] <- raw2clean %>%
  filter(Diet=="Highfibre",Probiotic=="Yes")%>%
  summarize(DietHi_ProYes=mean(Shannon.diversity))
# 6.276 - Correct

raw2means[3] <- raw2clean %>%
  filter(Diet=="Normal",Probiotic=="No")%>%
  summarize(DietNorm_ProNo=mean(Shannon.diversity))
# 2.105 - Correct

raw2means[4] <- raw2clean %>%
  filter(Diet=="Normal",Probiotic=="Yes")%>%
  summarize(DietNorm_ProYes=mean(Shannon.diversity))
# 4.439 - Correct

# Get F value for interaction term
drop1(ANOVA_2_1,test="F") # F=15.287

# Run post-hoc Tukey test to look for significant
# differences in group means
TukeyHSD(ANOVA_2_2)
# Results:
#                               diff       lwr        upr     p adj
# Yes:Highfibre-No:Highfibre  3.4687701  2.931626  4.0059142 0.0000000
# No:Normal-No:Highfibre     -0.7022299 -1.239374 -0.1650858 0.0049332
# Yes:Normal-No:Highfibre     1.6321034  1.094959  2.1692475 0.0000000
# No:Normal-Yes:Highfibre    -4.1710000 -4.703573 -3.6384275 0.0000000
# Yes:Normal-Yes:Highfibre   -1.8366667 -2.369239 -1.3040941 0.0000000
# Yes:Normal-No:Normal        2.3343333  1.801761  2.8669059 0.0000000

# All pairwise comparisons between group means show
# significant differences

# Participants on probiotics and a high fibre diet
# have a Shannon score...
### ... 3.47 units higher than those WITHOUT probiotics
### but who are on a high fibre diet
### ... 4.17 units higher than those WITHOUT probiotics
### and who are NOT on a high fibre diet
### ... 1.84 units higher than those on probiotics
### but who are NOT on a high fibre diet

# In order to have a high Shannon score, take probiotics
# and follow a high fibre diet, but the probiotic has
# a larger effect

##### PLOT GENERATION #####
# Draw an interaction plot to show differences in
# group means

interaction.plot(raw2clean$Diet,raw2clean$Probiotic,
                 raw2clean$Shannon.diversity,
                 type="b",trace.label="Treatment",
                 xlab="Diet",
                 ylab="Group mean for mean Shannon diversity",
                 lty=1,
                 col=c("red","royalblue"),
                 pch=c(18,18),
                 ylim=c(1,7))

# Add group means to plot
text(x=1,y=2.8069+0.2,labels="2.81")
text(x=1,y=6.2757+0.2,labels="6.28")
text(x=2,y=4.439+0.2,labels="4.44")
text(x=2,y=2.1047+0.2,labels="2.10")

# Draw group mean bar plot with asterisks
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)

# Make "Groups" column for "toplot" dataframe
Groups <- c(rep("Hi,Yes",30),rep("Nor,Yes",30),rep("Hi,No",29),rep("Nor,No",30))

# Create "toplot" dataframe:
# Combine Shannon diversity raw results with "Groups"
toplot <- data.frame(raw2clean$Shannon.diversity,Groups)

# Name columns of "toplot" dataframe
colnames(toplot) <- c("Mean_Shannon_diversity","Group")

# Define pariwise comparisons
my_comparisons = list( c("Hi,Yes", "Hi,No"), c("Nor,No", "Hi,No"), 
                       c("Nor,Yes", "Hi,No"), c("Nor,No", "Hi,Yes"),
                       c("Nor,Yes", "Hi,Yes"), c("Nor,Yes","Nor,No"))

# Plot boxplot with pairwise comparisons
ggboxplot(toplot, x = "Group", y ="Mean_Shannon_diversity")+ 
  stat_compare_means(comparisons = my_comparisons, 
                     label="p.signif",label.y = c(10, 11, 12, 13, 14, 15, 16))+
  stat_compare_means(label.y=18,method="anova")
