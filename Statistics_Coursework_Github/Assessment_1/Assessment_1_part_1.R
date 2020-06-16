# Assessment 1: Probiotics and synbiotics
# Investigation 1

##### EXPLORATORY DATA ANALYSIS #####

# Read in raw data file
raw1 <- read.table("Experiment1_data_set 29 .txt")

# Check structure
str(raw1) # As expected, 80 obs of 3 variables, coded
# as Probiotic: Factor w/ 2 levels "No", "Yes",
# Days: int, and "AMP": num. Further checking
# to follow:

# Describe variables in relation to present study,
# and sanity-check raw data:

# RESPONSE VARIABLE, continuous (num):
# Participants' recovery from GI infection was measured 
# by quantifying the antimicrobial pepties (AMP) in five
# intestinal mucus samples and taking a MEAN. Units 
# = mg/ml.
# View all data points and check length:
raw1$AMP # All numbers recorded to 2dp
length(raw1$AMP) # 80 obs, as expected
# Have a quick look at the distribution of AMP means:
range(raw1$AMP) # Range is 55.84-107.79
hist(raw1$AMP,breaks=50) # Positively skewed norm dist

# EXPLANATORY VARIABLE 1, Factor, 2 levels:
# Participants were given a probiotic or a placebo.
# This is indicated in the raw data file as factor
# "Probiotic" which has two levels "No" and "Yes".
# View all data points, levels and check length
raw1$Probiotic # Two levels, as expected
length(raw1$Probiotic) # 80 obs, as expected
# Check that there are 40 "Yes" and 40 "No"
plot(raw1$Probiotic) # Both 40, as expected

# EXPLANATORY VARIABLE 2, continuous (int):
# Participants' illness severity was indicated by the
# pre-treatment duration of diarrhoea and vomiting symptoms,
# measured in number of days. This data is recorded in the 
# raw data file as variable "Days"; data type is "integer" 
# / "int".
# View all data points and check length
raw1$Days # All integers, some 0s
length(raw1$Days) # 80 obs, as expected
# Have a quick look at the distribution of illness duration
range(raw1$Days) # Range is 0-19
hist(raw1$Days, breaks=50) # Fairly haphazard distribution

# No missing values, no "NA"s, for any variable

# Compare response variable to each explanatory variable:

### PROBIOTIC
boxplot(raw1$AMP~raw1$Probiotic)
# AMP was higher in participants who took probiotic
# Distributions are fairly normal

### SEVERITY OF ILLNESS
plot(raw1$AMP~raw1$Days)
# At first glance, it looks like there is a weak positive
# correlation between length of illness in days and AMP
# levels. It looks as though the error variance might increase
# as "Days" increases, but this will be checked later.

##### CHOICE OF MODEL #####

# As the response variable is continuous and broadly normally
# distributed, one explanatory variable is continuous, and the 
# other explanatory variable is a factor, use ANCOVA.

# Fit model. Include interaction term in first attempt.
ANCOVA_1_1 <- lm(raw1$AMP~raw1$Days*raw1$Probiotic)
# Object naming format:
# Model type ANCOVA
# 1 = Dataset 1
# 1 = First attempt at fitting a model

# Look at diagnostic plots:
plot(ANCOVA_1_1) # Residuals vs Fitted: homoscedastic
# QQ plot looks OK
# No influential outliers according to Cook's distance

# Look at distribution of errors as histogram
hist(ANCOVA_1_1$residuals) # Fairly normal

# Check whether interaction term is significant
drop1(ANCOVA_1_1,test="F")
# Interaction term is significant
# F value is 22.036, p-value is 1.162e-05
# Minimum adequate model has been reached
# Retain main effects terms in model, but don't include 
# in reporting statement
# Omit DF from statement as partial F test was used

# Look at summary table to get intercepts and slopes:
summary(ANCOVA_1_1)

# "Intercept" must be "ProbioticNo"
# "ProbioticNo": Intercept is 65.9045
# "ProbioticNo": slope is 0.3337
# "ProbioticYes": Intercept is 70.6042 (65.9045+4.6997)
# "ProbioticYes": slope is 1.4876 (0.3337+1.1539)

# Adjusted R-squared:  0.7018
# Multiple R-squared: 0.7131

##### PLOT GENERATION #####
# Install and load ggplot2
install.packages("ggplot2")
library(ggplot2)

# Plot x axis (Days), y axis (AMP), slopes (ProbioticNo
# and ProbioticYes)

ggplot(raw1,aes(x=Days,y=AMP,color=Probiotic))+
  # Plot raw1 data, use different colours for factor levels
  geom_point()+ # Plot data points
  # Label the axes
  labs(y="Mean AMP (mg/ml)",x="Symptom duration (days)")+
  # Legend settings:
  scale_color_manual(name="",
                     # set legend name
                     labels=c("Placebo","Probiotic"),
                     # set key names
                     values=c(Yes="#9400D3",No="#006400"),
                     # Set key colours
                     guide=guide_legend(reverse=TRUE))+
                     # Place "Probiotic" first in legend
  # Colour settings:
  theme(panel.background=element_rect(fill="white"))+
  # change background colour
  theme(axis.title.x=element_text(size=15))+
  # change x axis label size
  theme(axis.title.y=element_text(size=15))+
  # change y axis label size
  theme(legend.title=element_text(size=15))+
  # change legend title size
  theme(legend.text=element_text(size=12))+
  # change key name sizes
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  # remove gridlines
  geom_smooth(method="lm", # Add regression lines
              alpha=0.3,show.legend=FALSE)+
  # Include visual representation of 95% confidence intervals 
  # for predicted AMP values based on this linear model
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20))+
  scale_y_continuous(breaks=c(0,55,60,65,70,75,80,
                              85,90,95,100,105,110))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))
  # Customise ticks on axes