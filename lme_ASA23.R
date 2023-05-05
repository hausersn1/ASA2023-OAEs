# Linear mixed models

## Example 1, trying to use random rats data
dat <- read.csv("/Volumes/SNH/Projects/ASA/data.csv")
str(dat)

library(ggplot2)
# factors set for categorical variables
dat$OAEtype <- as.factor(dat$OAEtype)
dat$Group <- as.factor(dat$Group)
dat$chin <- as.factor(dat$chin)
dat$Frequency <- as.factor(dat$Frequency)

p <- ggplot(aes(x=OAEtype, y=dBSPL, color=Group), data=dat) + geom_boxplot()
p + xlab('Frequency') + ylab('Amplitude')

library(lme4) # load in lme
library(car)

m <- lmer(dBSPL ~ OAEtype*Group + Frequency + (1|chin), data=dat) 

# after you fit a model
summary(m) # will give you more info about what your model has done
Anova(m, test.statistic='F')

# make a model with only frequency so that it can be regressed out
m1 <- lm(dBSPL ~ Frequency, data=dat)

dat$dBSPL_freq2 <- residuals(m1)

# Plot data after regressing out frequency
p <- ggplot(aes(x=OAEtype, y=dBSPL_freq2, color=Group), data=dat) + geom_boxplot()
p + xlab('OAE type') + ylab('Amplitude (dB SPL)')
