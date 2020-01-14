# This code runs the regression analyses on Cam-CAN motor adaptation and measures of interest
# and plots figures 
# written by Noham Wolpe, nw305@medschl.cam.ac.uk

rm(list=ls(all=TRUE))
library(lavaan)

load('~/Dropbox/MotorLearning/VMLdata.Rda')
mydata <- VMLdata
# scaling all variables
adaptation <- scale(as.numeric(mydata$adaptation))[,1]
age <- scale(as.numeric(mydata$age))[,1]
vstmRMS <- scale(as.numeric(mydata$RMS_mean))[,1]
vstmRMS_int <- scale(as.numeric(mydata$RMS_int))[,1]
education <- scale(as.numeric(mydata$education))[,1]
bias <- scale(as.numeric(mydata$preexposure_epoch))[,1]
handedness <- scale(as.numeric(mydata$handedness))[,1]
gender <- as.factor(mydata$gender)
DLTM <- scale(as.numeric(mydata$DLTM))[,1]
DLTM_int <- scale(as.numeric(mydata$DLTM_int))[,1]

##### run main behavioural regressions #####

dataForReg <- data.frame(education, bias, handedness, gender, adaptation, age, 
                         DLTM, DLTM_int, vstmRMS, vstmRMS_int)

reg_model_DLTM <- 'adaptation ~ bias + education + handedness + gender + age + DLTM + DLTM_int'
fit_reg_model_DLTM <- sem(reg_model_DLTM, data = dataForReg, estimator = "ml",  missing = 'fiml')
summary(fit_reg_model_DLTM)
inspect(fit_reg_model_DLTM, 'r2')

reg_model_STM <- 'adaptation ~ bias + education + handedness + gender + age + vstmRMS + vstmRMS_int'
fit_reg_model_STM <- sem(reg_model_STM, data = dataForReg, estimator = "ml",  missing = 'fiml')
summary(fit_reg_model_STM)
inspect(fit_reg_model_STM, 'r2')


# run supplementary regression analyses
FOC <- scale(as.numeric(mydata$FOC))[,1]
FOC_int <- scale(mydata$FOC_int)[,1]

emoMem <- scale(as.numeric(mydata$expDetSum))[,1]
emoMem_int <- scale(as.numeric(mydata$emoMem_int))[,1] 

cattell <- scale(as.numeric(mydata$cattellTotal))[,1]
cattell_int <- scale(as.numeric(mydata$cattell_int))[,1]

dataForSuppReg <- data.frame(education, bias, handedness, gender, adaptation, age,
                             FOC, FOC_int, emoMem, emoMem_int, cattell, cattell_int)

reg_model_emo <- 'adaptation ~ bias + education + handedness + gender + age + emoMem + emoMem_int'
fit_reg_model_emo <- sem(reg_model_emo, data = dataForSuppReg, estimator = "ml",  missing = 'fiml')
summary(fit_reg_model_emo)

reg_model_cattell <- 'adaptation ~ bias + education + handedness + gender + age + cattell + cattell_int'
fit_reg_model_cattell <- sem(reg_model_cattell, data = dataForSuppReg, estimator = "ml",  missing = 'fiml')
summary(fit_reg_model_cattell)

reg_model_FOC <- 'adaptation ~ bias + education + handedness + gender + age + FOC + FOC_int'
fit_reg_model_FOC <- sem(reg_model_FOC, data = dataForSuppReg, estimator = "ml",  missing = 'fiml')
summary(fit_reg_model_FOC)


##### Plot Figures ####
library(ggplot2)

# Fig 1

ageFactorStruc = vector(mode="double", length=length(mydata$age))
for (i in 1:length(mydata$age)){
  if (mydata$age[i] < 46)
    ageFactorStruc[i] = "Young"
  else if (mydata$age[i] >= 66)
    ageFactorStruc[i] = "Old"
  else ageFactorStruc[i] = "Middle"
}
ageFactorStruc = factor(ageFactorStruc, levels = c("Young", "Middle", "Old"))
length(which(ageFactorStruc=="Young"))

trajErrorCycle = read.csv('~/Dropbox/MotorLearning/TrajErrorCycle.csv', header = FALSE)

youngMean = apply(trajErrorCycle[ageFactorStruc=='Young',], 2, mean)
youngSE = apply(trajErrorCycle[ageFactorStruc=='Young',], 2, sd)/sqrt(length(youngMean))
middleMean = apply(trajErrorCycle[ageFactorStruc=='Middle',], 2, mean)
middleSE = apply(trajErrorCycle[ageFactorStruc=='Middle',], 2, sd)/sqrt(length(middleMean))
oldMean = apply(trajErrorCycle[ageFactorStruc=='Old',], 2, mean)
oldSE = apply(trajErrorCycle[ageFactorStruc=='Old',], 2, sd)/sqrt(length(oldMean))
cycle = as.numeric(rep(c(1:48), 3))

allMean = matrix(data = c(cycle, 
                          as.numeric(youngMean), as.numeric(middleMean), as.numeric(oldMean), 
                          as.numeric(youngSE), as.numeric(middleSE), as.numeric(oldSE), 
                          rep('Young', 48), rep('Middle', 48), rep('Old', 48)), ncol = 4)

colnames(allMean) <- c('Cycle', 'TrajErrorMean', 'TrajErrorSE', 'Age')

df.fig1 <- data.frame(allMean)
df.fig1$Age = factor(df.fig$Age, levels = c("Young", "Middle", "Old"))
fig1 <- ggplot(data = df.fig1, aes(x = as.numeric(paste(Cycle)), y = as.numeric(paste(TrajErrorMean)), color = Age)) + 
  geom_ribbon(data = df.fig1, aes(ymin = as.numeric(paste(TrajErrorMean)) - as.numeric(paste(TrajErrorSE)), 
                                 ymax = as.numeric(paste(TrajErrorMean)) + as.numeric(paste(TrajErrorSE)), fill = Age), alpha=0.5, linetype = 0.1)
fig1 <- fig1 + geom_line() + geom_point()
fig1 <- fig1 + xlab('Cycle') + ylab('Mean Trajectory Error (degrees)')
print(fig1)

ggsave('~/Dropbox/MotorLearning/Figures/Fig1.pdf', plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300)

# now fig 2

library("ggstatsplot")
fig2 <- ggstatsplot::ggscatterstats(
  data = mydata, 
  x = age, 
  y = adaptation,
  title = "Final adaptation by age",
  messages = FALSE
)

ggsave('~/Dropbox/MotorLearning/Figures/Fig2.pdf', plot = fig2, device = "pdf", path = NULL,
       scale = 1, width = 7.44, height = 7.5, units = c("in"),
       dpi = 300)

# fig 3
setwd('~/Dropbox/MotorLearning/Figures/')
load('GMtoPlot.Rda')
names(GMtoPlot) <- c("adaptation", "age", "GM")

ageFactorStruc_imaging = vector(mode="double", length=length(GMtoPlot$age))
for (i in 1:length(GMtoPlot$age)){
  if (GMtoPlot$age[i] <= 46)
    ageFactorStruc_imaging[i] = "Young"
  else if (GMtoPlot$age[i] >= 66)
    ageFactorStruc_imaging[i] = "Old"
  else ageFactorStruc_imaging[i] = "Middle"
}
ageFactorStruc_imaging = factor(ageFactorStruc_imaging, levels = c("Young", "Middle", "Old"))
length(which(ageFactorStruc_imaging=="Old"))
GMtoPlot$ageFactor <- ageFactorStruc_imaging

fig3c <- ggplot(GMtoPlot, aes(x = GM, y = adaptation, colour = ageFactor))
fig3c <- fig3c + geom_point() + stat_smooth(method = "lm", se = TRUE)
fig3c <- fig3c + labs(x = "Mean Grey Matter Volume", y = "Final adaptation (degrees)", color = "Age Group\n")
fig3c <- fig3c + theme(plot.title = element_text(hjust = 0.5))
print(fig3c)

ggsave('~/Dropbox/MotorLearning/Figures/Fig3C.pdf', plot = fig3c, device = "pdf", path = NULL,
       scale = 1, width = 7.44, height = 7, units = c("in"),
       dpi = 300)

# and fig 4
ageRaw <- mydata$age
DLTMraw <- mydata$DLTM
vstmRaw <- mydata$RMS_mean
adaptationRaw <- mydata$adaptation
dataFig4 = data.frame(adaptationRaw, ageRaw, DLTMraw, vstmRaw, ageFactorStruc)

setwd('~/Dropbox/MotorLearning/Figures/')
tiff("fig4a.tiff", units="in", width=5, height=5, res=150)
fig4a <- ggplot(dataFig4, aes(x = DLTMraw, y = adaptationRaw, colour = ageFactorStruc))
fig4a <- fig4a + geom_point() + stat_smooth(method = "lm", se = TRUE)
fig4a <- fig4a + labs(x = "Declarative long-term memory", y = "Final adaptation (degrees)", color = "Age Group\n")
fig4a <- fig4a + theme(plot.title = element_text(hjust = 0.5))
fig4a
dev.off()

# ggsave('~/Dropbox/MotorLearning/Figures/Fig4a', plot = fig4a, device = "pdf", path = NULL,
#        scale = 1, width = 7, height = 7.2, units = c("in"),
#        dpi = 300)

# pdf("fig4b.pdf", units="in", width=5, height=5, res=150)
fig4b <- ggplot(dataFig4, aes(x = vstmRaw, y = adaptationRaw, colour = ageFactorStruc))
fig4b <- fig4b + geom_point() + stat_smooth(method = "lm", se = TRUE)
fig4b <- fig4b + labs(x = "Visual short-term memory error", y = "Final adaptation (degrees)", color = "Age\n")
fig4b <- fig4b + theme(plot.title = element_text(hjust = 0.5))
fig4b
dev.off()

ggsave('~/Dropbox/MotorLearning/Figures/Fig4b.pdf', plot = fig4b, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 7.2, units = c("in"),
       dpi = 300)
