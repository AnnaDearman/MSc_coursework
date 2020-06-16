# Assessment 2
# Dataset 3: HIV viral load and within-patient 
# population dynamics

library(dplyr)
library(ggplot2)

##### EXPLORATORY DATA ANALYSIS #####

# Read in raw data file
raw3 <- read.table("Experiment3_data_set29.txt")

# Check structure
str(raw3) # As expected, 40 obs of 5 variables,
# Variables are "num" or "Factor", as appropriate.
# Further checking to follow:

# Describe variables in relation to present study,
# and sanity-check raw data:

# There are 40 observations from one HIV+ 
# participant over 40 weeks, to monitor CNS viral
# load. The observations are in chronological
# order.

# RESPONSE VARIABLE, continuous: Viral load.
# Measured in log10(number of viral particles per ml)
# to reflect the abundance of virus in the sample 
# (taken from brain or spinal cord)
raw3$viral_load # numbers, all to 7dp, one negative
length(raw3$viral_load) # 40 obs as expected
range(raw3$viral_load) # Range -0.54-23.13
hist(raw3$viral_load,breaks=20) # normal-ish, with two
# much higher results

# EXPLANATORY VARIABLE 1, factor, 2 levels:

# CD4+ cell counts were taken, from blood, as a 
# snapshot of immune function at the time of 
# sampling. This was recorded as simply "lo" or
# "hi"
raw3$CD4 # 2 levels as expected
length(raw3$CD4) # 40 obs as expected
plot(raw3$CD4) # 20 of each

# Compare to response variable
boxplot(raw3$viral_load~raw3$CD4)
# Norm dist, much longer whiskers for low CD4
# Median viral load is higher for low CD4, as
# one might expect in HIV
# Keep in model, as it seems reasonable that an
# immune cell could play a role in immunoprotection
# of the blood-brain barrier, and HIV attacks
# the CD4+ cells, so it could have some biological 
# significance in this model

# Compare CD4 and covariate Shannon diversity:
boxplot(raw3$Shannon_diversity~raw3$CD4)
# Spread quite different for hi and low
# This might signify that it should be included
# in the model, although means are similar

# Compare response variable and covariate Shannon,
# with CD4 levels different colours
ggplot(raw3,aes(x=Shannon_diversity,y=viral_load,col=CD4))+
  geom_point()
# No visible clustering effect by CD4 level here

# Compare CD4 and covariate distance
boxplot(raw3$Distance~raw3$CD4)
# Bigger genetic distance IQR for high CD4.
# Mean genetic distance is higher for low CD4.
# Retain in model as group medians are different.

# Compare response variable and covariate Distance,
# with CD4 levels different colours
ggplot(raw3,aes(x=Distance,y=viral_load,col=CD4))+
  geom_point()
# No visible clustering effect by CD4 level here

# EXPLANATORY VARIABLE 2, continuous: Distance.

# This is the mean "genetic distance" after pairwise
# comparisons of each virus in the sample to a
# reference sample. aka "evolutionary distance".
raw3$Distance # numbers, all to 6dp
length(raw3$Distance) # 40 obs as expected
range(raw3$Distance) # Range 4.20-13.37
hist(raw3$Distance,breaks=20) # haphazard distribution,
# almost a negatively skewed normal

# Compare to response variable
plot(raw3$viral_load~raw3$Distance)
# Weak positive correlation, keep in model

# Compare to other covariate
plot(raw3$Distance~raw3$Shannon_diversity)
# No correlation

# EXPLANATORY VARIABLE 3, continuous: 

# Shannon diversity to represent diversity of 
# HIV strains detectable by PCR of Env sequence
raw3$Shannon_diversity # numbers, all to 7dp
length(raw3$Shannon_diversity) # 40 obs as expected
range(raw3$Shannon_diversity) # Range 0.16-9.54, no 
# negative values
hist(raw3$Shannon_diversity,breaks=20) # haphazard
# distribution

# Compare to response variable
plot(raw3$viral_load~raw3$Shannon_diversity)
# Possibly a very weak curvilinear relationship
# but with two high values for both variables
# Try transforming one variable to make relationship
# linear
plot(sqrt(raw3$viral_load)~raw3$Shannon_diversity)
plot(raw3$viral_load~sqrt(raw3$Shannon_diversity))
plot(log(raw3$viral_load)~raw3$Shannon_diversity)
plot(raw3$viral_load~log(raw3$Shannon_diversity))
plot(raw3$viral_load^2~raw3$Shannon_diversity)
plot(raw3$viral_load~raw3$Shannon_diversity^2)
# None of these make much difference, try
# polynomial
test1 <- lm(raw3$viral_load~
              poly(raw3$Shannon_diversity,
                   degree=3))

plot(raw3$viral_load~raw3$Shannon_diversity)
lines(smooth.spline(raw3$Shannon_diversity,predict(test1)),col="red")

# This line seems to fit the trend, and thus
# suggests I've found a potentially 
# useful transformation for later

# EXPLANATORY VARIABLE 4, factor, 2 levels:

# Tissue from which sample was taken for viral load,
# Shannon diversity and genetic distance testing.
# Brain or spinal cord.
raw3$Tissue # 2 levels as expected
length(raw3$Tissue) # 40 obs as expected
plot(raw3$Tissue) # 20 of each
boxplot(raw3$viral_load~raw3$Tissue)
# Median viral load is about the same for
# each tissue
# Data looks more spread for spinalcord
# Check spread using SDs in dplyr
raw3 %>%
  filter(Tissue=="spinalcord") %>%
  summarize(SD=sd(viral_load)) # SD 5.78
raw3 %>%
  filter(Tissue=="brain") %>%
  summarize(SD=sd(viral_load)) # SD 5.21
# Not a huge difference in SD

# Omit from model? Check against other continuous
# variables

boxplot(raw3$Shannon_diversity~raw3$Tissue)
# Tissue doesn't have much effect

boxplot(raw3$Distance~raw3$Tissue)
# Tissue doesn't have much effect

ggplot(raw3,aes(x=Shannon_diversity,y=viral_load,col=Tissue))+
  geom_point()
# Tissue doesn't have much effect

ggplot(raw3,aes(x=Distance,y=viral_load,col=Tissue))+
  geom_point()
# Tissue doesn't have much effect

# Omit Tissue from model

# Re-check plots simultaneously to rule out colinear 
# variables
pairs(raw3)
# No particularly colinear variables

# GOING FORWARD...

# Retain all data points, as none are wildly out of
# range so as to suggest a recording error.

# Retain all terms except Tissue

##### CHOICE OF MODEL #####

# Viral load is response variable
# Keep CD4, Shannon and Distance in model
# shannon needs to be transformed
# poly(raw3$Shannon_diversity,degree=3)
# Try interaction term for CD4 and Shannon
# Try interaction term for CD4 and Distance

# Fit model 1

GLM_3_1 <- lm(viral_load~
                CD4+poly(Shannon_diversity,degree=3)+Distance+
                CD4:Shannon_diversity+
                CD4:Distance,
              data=raw3)
# Object naming format:
# Model type GLM
# 3 = Dataset 3
# 1 = First attempt at fitting a mode

# Check diagnostic plots

plot(GLM_3_1) # Not good Fitted vs Redisuals
# QQ not good or terrible
# Cook's distance good, no influential outliers
# which was a concern
hist(GLM_3_1$residuals) # Quite normally dist
drop1(GLM_3_1,test="F")
# It appears that non-significant terms include:
#                                      Df Sum of Sq    RSS    AIC F value  Pr(>F)  
# poly(Shannon_diversity, degree = 3)  2     8.485 605.51 120.69  0.2274 0.79789  
# CD4:Shannon_diversity                1    12.726 609.75 122.97  0.6821 0.41498  
# CD4:Distance                         1    54.553 651.58 125.62  2.9240 0.09695 .

# Remove non-significant terms
# Fit model 2

GLM_3_2 <- lm(viral_load~
                CD4+Shannon_diversity+Distance,
              data=raw3)

plot(GLM_3_2) # Residuals v fitted still not good
# QQ looks better in the middle
# Still no unduly influential outliers
hist(GLM_3_2$residuals) # Looks more skewed than before
drop1(GLM_3_2,test="F")
# It appears that non-significant terms include:
#                    Df Sum of Sq     RSS    AIC F value    Pr(>F)    
# CD4                1     21.87  734.20 122.40  1.1052 0.3001309    
# Shannon_diversity  1     30.79  743.12 122.88  1.5561 0.2202844    

# which would leave us with just

# Distance           1    375.55 1087.88 138.12 18.9796 0.0001052

# Fit model 3

GLM_3_3 <- lm(viral_load~Distance,data=raw3)

plot(GLM_3_3) # Residuals v fitted slightly better
# but not ideal
# QQ plot still not ideal
# Cook's distance plot still good
hist(GLM_3_3$residuals,breaks=10) # not great or terrible
drop1(GLM_3_3,test="F") # significant!

#           Df Sum of Sq     RSS    AIC F value   Pr(>F)    
# Distance  1    386.26 1153.18 136.46  19.139 9.15e-05 ***

##### USE GLM_3_3 <- lm(viral_load~Distance,data=raw3)

anova(GLM_3_3) # Get statistics for reporting
# F(1,38) = 19.139, p<0.001

summary(GLM_3_3) # Get adjusted R-squared
# 0.318
# and intercept: -2.4325

# Make presentable graph

ggplot(data=raw3,aes(x=Distance,y=viral_load))+
  geom_point(pch=17)+
  geom_smooth(method="lm",col="black")+
  labs(y="Viral load (log10(number of viral particles per ml)",
       x="Mean pairwise genetic distance from reference sequence")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(panel.background=element_rect(fill="white"))+
  ylim(c(-1,25))+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))

# Test for model assumptions as fitted v residuals 
# plot was not ideal
install.packages("gvlma")
library(gvlma)
gvmodel <- gvlma(GLM_3_3)
summary(gvmodel)

#                    Value p-value                Decision
# Global Stat        6.8732  0.1427 Assumptions acceptable.
# Skewness           0.1914  0.6617 Assumptions acceptable.
# Kurtosis           1.9317  0.1646 Assumptions acceptable.
# Link Function      2.2348  0.1349 Assumptions acceptable.
# Heteroscedasticity 2.5153  0.1127 Assumptions acceptable.