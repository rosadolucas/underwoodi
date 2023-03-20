###  Script by Lucas Rosado for master's analysis
###  March 2021 to March 2023
###  Mendonça et al. In Prep. Effects of the urban environment and climate change on a neotropical parthenogenetic lizard.

# 0. Main Directory and Packages ------------------------------------------
{ 
  if(!require(ggpubr)) install.packages("ggpubr", dependencies = T)
  if(!require(tidyverse)) install.packages('tidyverse', dependencies = T)
  if(!require(DHARMa)) install.packages('DHARMa', dependencies = T)
  if(!require(cowplot)) install.packages('cowplot', dependencies = T)
  if(!require(remotes)) install.packages('remotes', dependencies = T)
  if(!require(gratia)) remotes::install_github('gavinsimpson/gratia')
  if(!require(lme4)) install.packages('lme4', dependencies = T)
  if(!require(vroom)) install.packages('vroom', dependencies = T)
  if(!require(rlang)) install.packages('rlang', dependencies = T)
  if(!require(forcats)) install.packages('forcats', dependencies = T)
  if(!require(regclass)) install.packages('regclass', dependencies = T)
  if(!require(RRPP)) install.packages('RRPP', dependencies = T)
  if(!require(ggplot2)) install.packages('ggplot2', dependencies = T)
  if(!require(ggbiplot)) install.packages('ggbiplot', dependencies = T)
  if(!require(heatmaply)) install.packages('heatmaply', dependencies = T)
  if(!require(Rmisc)) install.packages('Rmisc', dependencies = T)
  if(!require(car)) install.packages('car', dependencies = T)
  if(!require(factoextra)) install.packages('factorextra', dependencies = T)
  if(!require(FactoMineR)) install.packages('FactorMineR', dependencies = T)
  if(!require(vegan)) install.packages('vegan', dependencies = T)
  if(!require(sjPlot)) install.packages('sjPlot', dependencies = T)
  if(!require(lmtest)) install.packages('lmtest', dependencies = T)
  if(!require(mgcv)) install.packages('mgcv', dependencies = T)
  if(!require(itsadug)) install.packages('itsadug', dependencies = T)
  if(!require(wiqid)) install.packages('wiqid', dependencies = T)
  if(!require(GGally)) install.packages('GGally', dependencies = T)
  if(!require(glmmTMB)) install.packages('glmmTMB', dependencies = T)
} # Packages

mainDIR <- "D:/Lucas/Underwoodi"

setwd(mainDIR)

rm(list=ls()) # clean environment

# 1. Morphology >>>>> Neonative x Native ---------------------------------------------
Morph <- read.csv(paste0(mainDIR, "/Underwoodi_morpho.csv"), header = T, sep = ",")
head(Morph)

Morph <- Morph %>% 
  select(-(BTAL), -(OBS), -(Maturity)) %>% 
  mutate_if(is.character, as.factor)

summary(Morph)
head(Morph)

### Removing NA's ---------------------------------------------------------------#
Morph <- Morph[-c(18, 21),] # LRM40 and 43
summary(Morph)

### Completing NA's with lm() ---------------------------------------------------#
{
  HH.M <- lm(HH ~ SVL + HL + HW, data = Morph)
  SG.M <- lm(SG ~ SVL + INF, data  = Morph)
  PG.M <- lm(PG ~ SVL + TRL, data = Morph)
  FTL_L.M <- lm(FTL_L ~ SVL + FL_L + TL_L, data = Morph)
  FTL_R.M <- lm(FTL_R ~ SVL + FL_R + TL_R, data = Morph)
  HUL_L.M <- lm(HUL_L ~ SVL, data = Morph)
  FAL_L.M <- lm(FAL_L ~ SVL, data = Morph)
  FAL_R.M <- lm(FAL_L ~ SVL, data = Morph)
  TAL.M <- lm(TAL ~ SVL, data = Morph)
  LFT_L.M <- lm(LFT_L ~ SVL + FL_L + TL_L + FTL_L, data = Morph)
  LFT_R.M <- lm(LFT_R ~ SVL + FL_R + TL_R + FTL_R, data = Morph)
  LFF_L.M <- lm(LFF_L ~ SVL + HUL_L + FAL_L, data = Morph)
  LFF_R.M <- lm(LFF_R ~ SVL + HUL_R + FAL_R, data = Morph)
} # lm()

{
  KnownFPW03738 <- Morph %>% 
    filter(Specimen == 'FPWerneck03738') %>% 
    select(-(Specimen:distribution), -(HH)) 
  
  KnownINPA_H040754_SG <- Morph %>% 
    filter(Specimen == 'INPA-H040754') %>% 
    select(-(Specimen:distribution), -(SG))
  
  KnownINPA_H040754_PG <- Morph %>% 
    filter(Specimen == 'INPA-H040754') %>% 
    select(-(Specimen:distribution), -(PG))
  
  KnownLRM43_FTL <- Morph %>% 
    filter(Specimen == 'LRM43') %>% 
    select(-(Specimen:distribution), -(FTL_L))
  
  KnownLRM43_FAL <- Morph %>% 
    filter(Specimen == 'LRM43') %>% 
    select(-(Specimen:distribution), -(FAL_L))
  
  KnownLRM43_FAL_R <- Morph %>% 
    filter(Specimen == 'LRM43') %>% 
    select(-(Specimen:distribution), -(FAL_R))
  
  KnownFPW03730_FTL <- Morph %>% 
    filter(Specimen == 'FPWerneck03730') %>% 
    select(-(Specimen:distribution), -(FTL_L))
  
  KnownFPW03730_FAL <- Morph %>% 
    filter(Specimen == 'FPWerneck03730') %>% 
    select(-(Specimen:distribution), -(FAL_L))
  
  KnownINPA_H031307 <- Morph %>% 
    filter(Specimen == 'INPA-H031307') %>% 
    select(-(Specimen:distribution), -(FTL_R))
  
  KnownLRM43_FTL_R <- Morph %>% 
    filter(Specimen == 'LRM43') %>% 
    select(-(Specimen:distribution), -(FTL_R))
  
  KnownLRM48_HUL <- Morph %>% 
    filter(Specimen == 'LRM48') %>% 
    select(-(Specimen:distribution), -(HUL_L))
  
  KnownLRM48_FAL <- Morph %>% 
    filter(Specimen == 'LRM48') %>% 
    select(-(Specimen:distribution), -(FAL_L))
  
  KnownFPWerneck03728 <- Morph %>% 
    filter(Specimen == 'FPWerneck03728') %>% 
    select(-(Specimen:distribution), -(TAL))
  
  KnownFPWerneck03728 <- Morph %>% 
    filter(Specimen == 'FPWerneck03728') %>% 
    select(-(Specimen:distribution), -(TAL))
  
} 

Morph_full <- Morph %>% 
  mutate(HH = ifelse(Specimen == 'FPWerneck03738', predict(HH.M, KnownFPW03738), HH),
         SG = ifelse(Specimen == 'INPA-H040754', predict(SG.M, KnownINPA_H040754_SG), SG),
         PG = ifelse(Specimen == 'INPA-H040754', predict(PG.M, KnownINPA_H040754_PG), PG),
         FTL_L = ifelse(Specimen == 'LRM43', predict(FTL_L.M, KnownLRM43_FTL), FTL_L),
         FTL_R = ifelse(Specimen == 'LRM43', predict(FTL_R.M, KnownLRM43_FTL_R), FTL_R),
         FAL_L = ifelse(Specimen == 'LRM43', predict(FAL_L.M, KnownLRM43_FAL), FAL_L),
         FAL_R = ifelse(Specimen == 'LRM43', predict(FAL_R.M, KnownLRM43_FAL_R), FAL_R),
         FTL_R = ifelse(Specimen == 'INPA-H031307', predict(FTL_R.M, KnownINPA_H031307), FTL_R),
         FTL_L = ifelse(Specimen == 'FPWerneck03730', predict(FTL_L.M, KnownFPW03730_FTL), FTL_L),
         FAL_L = ifelse(Specimen == 'FPWerneck03730', predict(FAL_L.M, KnownFPW03730_FAL), FAL_L),
         HUL_L = ifelse(Specimen == 'LRM48', predict(HUL_L.M, KnownLRM48_HUL), HUL_L),
         FAL_L = ifelse(Specimen == 'LRM48', predict(FAL_L.M, KnownLRM48_FAL), FAL_L), 
         FAL_L = ifelse(Specimen == 'FPWerneck03722', predict(FAL_L.M, KnownFPW03722_FAL), FAL_L),
         TAL = ifelse(Specimen == 'FPWerneck03728', predict(TAL.M, KnownFPWerneck03728), TAL))

{KnownINPAH031306_LFT_L <- Morph_full %>% 
    filter(Specimen == 'INPA-H031306') %>% 
    select(-(Specimen:distribution), -(LFT_L))
  
  KnownINPAH031306_LFT_R <- Morph_full %>% 
    filter(Specimen == 'INPA-H031306') %>% 
    select(-(Specimen:distribution), -(LFT_R))
  
  KnownINPAH031308_LFT_L <- Morph_full %>% 
    filter(Specimen == 'INPA-H031308') %>% 
    select(-(Specimen:distribution), -(LFT_L))
  
  KnownINPAH031308_LFT_R <- Morph_full %>% 
    filter(Specimen == 'INPA-H031308') %>% 
    select(-(Specimen:distribution), -(LFT_R))
  
  KnownINPAH031310_LFT_L <- Morph_full %>% 
    filter(Specimen == 'INPA-H031310') %>% 
    select(-(Specimen:distribution), -(LFT_L))
  
  KnownINPAH031310_LFT_R <- Morph_full %>% 
    filter(Specimen == 'INPA-H031310') %>% 
    select(-(Specimen:distribution), -(LFT_R))
  
  KnownLRM48_LFF_L <- Morph_full %>% 
    filter(Specimen == 'LRM48') %>% 
    select(-(Specimen:distribution), -(LFF_L))
  
  KnownFPWerneck03730_LFT_L <- Morph_full %>% 
    filter(Specimen == 'FPWerneck03730') %>% 
    select(-(Specimen:distribution), -(LFT_L))
  
  KnownFPWerneck03730_LFF_L <- Morph_full %>% 
    filter(Specimen == 'FPWerneck03730') %>% 
    select(-(Specimen:distribution), -(LFF_L))
  }

Morph_full <- Morph_full %>% 
  mutate(LFT_L = ifelse(Specimen == 'INPA-H031306', predict(LFT_L.M, KnownINPAH031306_LFT_L), LFT_L),
         LFT_R = ifelse(Specimen == 'INPA-H031306', predict(LFT_R.M, KnownINPAH031306_LFT_R), LFT_R),
         LFT_L = ifelse(Specimen == 'INPA-H031308', predict(LFT_L.M, KnownINPAH031308_LFT_L), LFT_L),
         LFT_R = ifelse(Specimen == 'INPA-H031308', predict(LFT_R.M, KnownINPAH031308_LFT_R), LFT_R),
         LFT_L = ifelse(Specimen == 'INPA-H031310', predict(LFT_L.M, KnownINPAH031310_LFT_L), LFT_L),
         LFT_R = ifelse(Specimen == 'INPA-H031310', predict(LFT_R.M, KnownINPAH031310_LFT_R), LFT_R),
         LFF_L = ifelse(Specimen == 'LRM48', predict(LFF_L.M, KnownLRM48_LFF_L), LFF_L),
         LFT_L = ifelse(Specimen == 'FPWerneck03730', predict(LFT_L.M, KnownFPWerneck03730_LFT_L), LFT_L),
         LFF_L = ifelse(Specimen == 'FPWerneck03730', predict(LFF_L.M, KnownFPWerneck03730_LFF_L), LFF_L))

summary(Morph_full)
Morph_2 <- Morph_full[,-c(6)]

### BM ------------------------------------#
Morph_BM <- Morph_full %>% 
  drop_na()
head(Morph_BM)      

summary(lm(BM ~ SVL, data = Morph_BM)) # *

rbm1 <- lm(Morph_BM$BM ~ Morph_BM$SVL)
Morph_BM$rbm1 <- residuals(rbm1)

# Normality
qqPlot(rbm1)
# Homocedacidade
shapiro.test(residuals(rbm1)) 
leveneTest(BM ~ distribution, data = Morph_BM) 

t.test(rbm1 ~ distribution, data = Morph_BM, var.equal = T)

ggplot(data = Morph_BM, aes(x = distribution, y = rbm1, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("Body mass (g) - ", italic("G. underwoodi")))) +
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.75) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.7) +
  scale_color_manual(values = c("black", "black")) +
  theme_light() +
  theme(legend.position = "none")

summarySE(Morph_BM, groupvars = 'distribution', measurevar = 'BM')

### Continuos
Morph_continuos <- Morph_2[,-(24:29)] 

### Relation SVL ~ traits ------------------------------------------------#
list <- colnames(Morph_continuos); list

for(i in 6:23){
  print(ggplot(data = Morph_continuos, aes(x = SVL, y = Morph_3[, i]))+
          geom_smooth(method = "lm", formula = y ~ x, se = T, 
                      fill = "red", alpha = 0.2, col = "red")+
          theme_classic()+
          geom_jitter(size = 2, alpha = .5)+
          labs(title = paste0(list[i])))
} 

summary(lm(Morph_continuos$TAL ~ Morph_continuos$SVL))

### t-test ---------------------------------------------------------------#
# Normality
residual <- lm(SVL ~ distribution, data = Morph_continuos)
qqPlot(residual)

# Homocedacidade
shapiro.test(residuals(residual))
leveneTest(SVL ~ distribution, data = Morph_continuos) 

hist(Morph_3$SVL)

# T-test
t.test(SVL ~ distribution, data = Morph_continuos, var.equal = T) 

ggplot(data = Morph_continuos, aes(x = distribution, y = SVL, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("SVL (mm) - ", italic("G. underwoodi")))) +
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.75) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.7) +
  scale_color_manual(values = c("black", "black")) +
  theme_light() +
  theme(legend.position = "none")


### PCA -------------------------------------------------------------------#
# Correlation 
Morph_3 <- Morph_continuos
Morph_4 <- drop_na(Morph_3)

summary(Morph_4)

Morph_5 <- Morph_4[,-4]
head(Morph_5)

Morph_residuals <- matrix(0, nrow=nrow(Morph_5), ncol=ncol(Morph_5)) 

for (y in 4:ncol(Morph_5)) {
  rl.f <- lm(Morph_5[,y] ~ Morph_5[,3])
  Morph_residuals[,y] <- as.array(rl.f$residuals)
}

head(Morph_residuals)

Morph_residuals_2 <- cbind(Morph_5[,1:3], Morph_residuals[,4:20]) 
head(Morph_residuals_2)

colnames(Morph_residuals_2) <- colnames(Morph_5[1:20]) 
head(Morph_residuals_2)

Morph_residuals_2[,c(21:24)] <- Morph_4[,c(26:29)] 
head(Morph_residuals_2)

Morph_residuals_2$TAL <- Morph_4$TAL 

Y_Residual <- Morph_residuals_2  %>% 
  select(-(Specimen:SVL)) %>%
  scale() %>% 
  data.matrix()

head(Y_Residual)

Morph_6 <- drop_na(Morph_3) 
summary(Morph_6)
summarySE(Morph_6, groupvars = 'distribution', measurevar = 'TAL')

# PCA
pca.res <- prcomp(Y_Residual, center = T, scale = T) 

ggbiplot(pca.res, obs.scale = 1, var.scale = 1,
         groups = Morph_5$distribution, ellipse = T,
         circle = F)+
  scale_color_discrete(name = '')+
  geom_hline(yintercept = 0, lty = 3, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 3, color = "grey", alpha = 0.9) +
  theme_classic()

var_env <- get_pca_var(pca.res)
var_env$coord 
var_env$contrib 
var_env$cor 

### NP-MANOVA (RRPP - Following Telemeco  & Gangloff 2020 ) -----------------#
# Using response matrix (Y_Residual)

# 1. rrpp.frame
MorphResidual_Standardized <- rrpp.data.frame(Id = Morph_residuals_2$Specimen,
                                              distribution = Morph_residuals_2$distribution,
                                              Y = Y_Residual)

head(MorphResidual_Standardized)

# 2. Making models and testing hypothesis
M1.1 <- lm.rrpp(Y ~ distribution, SS.type = 'I', data = MorphResidual_Standardized, print.progress = F, iter = 99999)

anova(M1.1, print.progress = T)
summary(M1.1, formula = T)

coef(M1.1, test = T) # examine model coefficients and compare their effect sizes

# Creating a dataframe with all combinations to collect the least-square means from model, for predict values
M1.Predictions <- expand.grid(distribution = levels(factor(MorphResidual_Standardized$distribution)))
rownames(M1.Predictions) <- M1.Predictions$distribution

M1.Predictions <- predict(M1.1, M1.Predictions, confidence = 0.95)

summary(M1.Predictions) # Least-square means for each response variable

# Acessing PCA 
M1.Predictions[['pca']]

# Least-squares mean 
PredMean <- data.frame(M1.Predictions$mean) %>%
  mutate(Treatment = row.names(M1.Predictions$mean),
         statistic = "LSmean")
PredLCL <- data.frame(M1.Predictions$lcl) %>%
  mutate(Treatment = row.names(M1.Predictions$mean),
         statistic = "LCL")
PredUCL <- data.frame(M1.Predictions$ucl) %>%
  mutate(Treatment = row.names(M1.Predictions$mean),
         statistic = "UCL")

PredVal <- rbind(PredMean, PredLCL, PredUCL) %>%
  gather(key = "Variable", value = "Value", -(Treatment:statistic)) %>%
  spread(key = statistic, value = Value)

ggplot(data = PredVal, aes(x = Variable, y = LSmean, fill = Treatment,
                           shape = Treatment, col = Treatment)) +
  rory_theme +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge(width = 0.5)) +
  geom_point(size = 5, col = "black", position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  scale_fill_manual(values = c("#A86D00", "#509AF9")) +
  scale_shape_manual(values = c(21,22,24)) +
  coord_flip() +
  theme(plot.title=element_text(hjust = 0.5, size = 20)) +
  labs(y = "Least-squares mean",
       x = "Morphological measurement",
       title = "Predicted values and 95% confidence intervals for
each morphological measurement")

# 2. Asymmetry >>>>> Neonative x Native ------------------------------------------------------------
asymmetry <- read.csv(paste0(mainDIR, "/Asymmetry.csv"), header = T, sep = ",")

head(asymmetry)

asymmetry <- asymmetry %>% 
  mutate_if(is.character, as.factor) %>% 
  filter(Specimen != 'INPA-H031306') %>% 
  filter(Specimen != 'INPA-H031307') %>%
  filter(Specimen != 'INPA-H031308') %>%
  filter(Specimen != 'INPA-H031310') %>%
  filter(Specimen != 'FPWerneck03730') %>%
  filter(Specimen != 'LRM43') %>%
  filter(Specimen != 'LRM48')

# Making a asymmetry index (AI) table 
Left <- filter(asymmetry, Side == 'Left') 
Right <- filter(asymmetry, Side == 'Right')

{AI <- matrix(0, nrow = nrow(Left), ncol = ncol(Left))
  AI[,5] <- (Right$FL - Left$FL)
  AI[,6] <- (Right$TL - Left$TL)
  AI[,7] <- (Right$FTL - Left$FTL)
  AI[,8] <- (Right$HUL - Left$HUL)
  AI[,9] <- (Right$FAL - Left$FAL)
  AI[,10] <- (Right$LFT - Left$LFT)
  AI[,11] <- (Right$LFF - Left$LFF)
}
AI <- as.data.frame(AI)

AI[,1:4] <- Left[,1:4]
colnames(AI) <- colnames(Left)

AI <- AI %>% 
  mutate_if(is.character, as.factor)

AI <- AI[,-4]

summary(AI)
head(AI)

# testing size dependence on body size of each variable
anova(lm(LFT ~ SVL, data = AI)) # all values has p > 0.05  

# testing size dependence on traid size of each variable
{RL2 <- matrix(0, nrow = nrow(Left), ncol = ncol(Left))
  RL2[,5] <- (Right$FL + Left$FL)/2
  RL2[,6] <- (Right$TL + Left$TL)/2
  RL2[,7] <- (Right$FTL + Left$FTL)/2
  RL2[,8] <- (Right$HUL + Left$HUL)/2
  RL2[,9] <- (Right$FAL + Left$FAL)/2
  RL2[,10] <- (Right$LFT + Left$LFT)/2
  RL2[,11] <- (Right$LFF + Left$LFF)/2
}
RL2 <- as.data.frame(RL2)

RL2[,1:4] <- Left[,1:4]
colnames(RL2) <- colnames(Left)

RL2 <- RL2[,-4]
summary(RL2)
head(RL2)

write.csv(AI, paste0(mainDIR, '/RL2.csv'))

anova(lm(AI$LFF ~ RL2$LFF)) # all values has p > 0.05

# Assessing FA and DA 
asymmetry_RR <- filter(asymmetry, distribution == 'Native')
asymmetry_AM <- filter(asymmetry, distribution == 'Neonative')

summary(asymmetry_RR)

lmm1 <- nlme::lme(LFT ~ Side + Specimen:Side,
                  random = ~ 1|Specimen, method = "REML",
                  data = asymmetry_RR)

anova(lmm1)
summary(lmm1)

qqnorm(lmm1$residuals)
qqline(lmm1$residuals)
hist(lmm1$residuals)


lmm2 <- glmmTMB(LFF ~ Side + Side:Specimen  
                + (1|Specimen),
                family = poisson,
                data = asymmetry_AM)

summary(lmm2)
car::Anova(lmm2, type = 'III')

# Ploting DA or FA into the pop
if(!require(visreg)) install.packages('visreg', dependencies = T)
visreg(lmm1, 'Side', by = 'Specimen', overlay = F ) # to see how the model perform

ggplot(data = AI, aes(x = LFT, group = distribution, fill = distribution))+
  geom_density(adjust = 1.5, alpha = 0.4)+
  scale_fill_manual(values = c("#A86D00", "#509AF9"))+
  geom_vline(xintercept = 0.0, size = 1.5, color = 'black')+
  ggtitle('DA LFT')+
  labs(x = "AI LFT")+
  theme_light()

ggplot(data = asymmetry_AM, aes(x = Side, y = LFT, color = Side)) + 
  labs(x = "distribution", 
       y = expression(paste("FL (mm) - ", italic("G. underwoodi")))) +
  ggtitle('FL AI')+
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.75) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.7) +
  stat_boxplot(geom = 'errorbar')+
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  theme_light() +
  theme(legend.position = "none")

ggplot(data = asymmetry, aes(x = distribution, y = LFT, color = Side)) + 
  ggtitle('LFF DA')+
  geom_violin(aes(fill = Side), alpha = 0.75, color = '#767A79', outlier.shape = NA)+
  scale_fill_manual(values = c('#C183FF', '#5EEB8B'))+
  scale_color_manual(values = c('#C183FF', '#5EEB8B'))+
  facet_wrap(~distribution, scale = 'free')+
  geom_point(position = position_jitterdodge(), shape = 16, cex = 3, alpha = 0.9)+
  stat_boxplot(geom = 'errorbar')+
  theme_light()

# degree of FA for each trait
Left_AM <- filter(Left, distribution == 'Neonative')
Left_RR <- filter(Left, distribution == 'Native')
Right_AM <- filter(Right, distribution == 'Neonative')
Right_RR <- filter(Right, distribution == 'Native')

Left_AM_mean <- aggregate(.~Specimen, FUN = mean, data = Left_AM[, -2])
Right_AM_mean <- aggregate(.~Specimen, FUN = mean, data = Right_AM[, -2])
Left_RR_mean <- aggregate(.~Specimen, FUN = mean, data = Left_RR[, -2])
Right_RR_mean <- aggregate(.~Specimen, FUN = mean, data = Right_RR[, -2])


RR_mean <- rbind(Left_RR_mean, Right_RR_mean)
RR_mean['Side'][RR_mean['Side'] == '1'] <- 'Left' 
RR_mean['Side'][RR_mean['Side'] == '2'] <- 'Right' 
head(RR_mean)

AM_mean <- rbind(Left_AM_mean, Right_AM_mean)
AM_mean['Side'][AM_mean['Side'] == '1'] <- 'Left' 
AM_mean['Side'][AM_mean['Side'] == '2'] <- 'Right' 
head(AM_mean)

{AI_AM <- matrix(0, nrow = nrow(Left_AM_mean), ncol = ncol(Left_AM_mean))
  AI_AM[,4] <- abs(log(Right_AM_mean$FL) -  log(Left_AM_mean$FL))
  AI_AM[,5] <- abs(log(Right_AM_mean$TL) -  log(Left_AM_mean$TL))
  AI_AM[,6] <- abs(log(Right_AM_mean$FTL) - log(Left_AM_mean$FTL))
  AI_AM[,7] <- abs(log(Right_AM_mean$HUL) - log(Left_AM_mean$HUL))
  AI_AM[,8] <- abs(log(Right_AM_mean$FAL) - log(Left_AM_mean$FAL))
  AI_AM[,9] <- abs(log(Right_AM_mean$LFT) - log(Left_AM_mean$LFT))
  AI_AM[,10] <- abs(log(Right_AM_mean$LFF) - log(Left_AM_mean$LFF))
}

{AI_RR <- matrix(0, nrow = nrow(Left_RR_mean), ncol = ncol(Left_RR_mean))
  AI_RR[,4] <- abs(log(Right_RR_mean$FL) -  log(Left_RR_mean$FL))
  AI_RR[,5] <- abs(log(Right_RR_mean$TL) -  log(Left_RR_mean$TL))
  AI_RR[,6] <- abs(log(Right_RR_mean$FTL) -  log(Left_RR_mean$FTL))
  AI_RR[,7] <- abs(log(Right_RR_mean$HUL) -  log(Left_RR_mean$HUL))
  AI_RR[,8] <- abs(log(Right_RR_mean$FAL) -  log(Left_RR_mean$FAL))
  AI_RR[,9] <- abs(log(Right_RR_mean$LFT) -  log(Left_RR_mean$LFT))
  AI_RR[,10] <- abs(log(Right_RR_mean$LFF) - log(Left_RR_mean$LFF))
}

AI_AM <- AI_AM %>% 
  as.data.frame()

AI_AM[,1:3] <- Left_AM_mean[,1:3]

colnames(AI_AM) <- colnames(Left_AM_mean)

AI_AM$distribution<- sample('AM', nrow(AI_AM), replace = T)

AI_RR <- AI_RR %>% 
  as.data.frame() 

AI_RR[,1:3] <- Left_RR_mean[,1:3]

colnames(AI_RR) <- colnames(Left_RR_mean)

AI_RR$distribution<- sample('RR', nrow(AI_RR), replace = T)

AI_abs <- rbind.data.frame(AI_AM, AI_RR) 

AI_abs <- AI_abs %>% 
  as.data.frame()

AI_abs <- AI_abs %>% 
  mutate_if(is.character, as.factor) 

AI_abs <- AI_abs[,-3]

head(AI_abs)
summary(AI_abs)

write.csv(AI_abs, paste0(mainDIR, '/AI_abs.csv'))

# AI average of two measures
{AI_AM <- matrix(0, nrow = nrow(Left_AM_mean), ncol = ncol(Left_AM_mean))
  AI_AM[,4] <- (Right_AM_mean$FL) -  (Left_AM_mean$FL)
  AI_AM[,5] <- (Right_AM_mean$TL) -  (Left_AM_mean$TL)
  AI_AM[,6] <- (Right_AM_mean$FTL) - (Left_AM_mean$FTL)
  AI_AM[,7] <- (Right_AM_mean$HUL) - (Left_AM_mean$HUL)
  AI_AM[,8] <- (Right_AM_mean$FAL) - (Left_AM_mean$FAL)
  AI_AM[,9] <- (Right_AM_mean$LFT) - (Left_AM_mean$LFT)
  AI_AM[,10] <- (Right_AM_mean$LFF) - (Left_AM_mean$LFF)
}

{AI_RR <- matrix(0, nrow = nrow(Left_RR_mean), ncol = ncol(Left_RR_mean))
  AI_RR[,4] <- (Right_RR_mean$FL) -  (Left_RR_mean$FL)
  AI_RR[,5] <- (Right_RR_mean$TL) -  (Left_RR_mean$TL)
  AI_RR[,6] <- (Right_RR_mean$FTL) -  (Left_RR_mean$FTL)
  AI_RR[,7] <- (Right_RR_mean$HUL) -  (Left_RR_mean$HUL)
  AI_RR[,8] <- (Right_RR_mean$FAL) -  (Left_RR_mean$FAL)
  AI_RR[,9] <- (Right_RR_mean$LFT) -  (Left_RR_mean$LFT)
  AI_RR[,10] <- (Right_RR_mean$LFF) - (Left_RR_mean$LFF)
}
AI_AM <- AI_AM %>% 
  as.data.frame()
AI_AM[,1:3] <- Left_AM_mean[,1:3]
colnames(AI_AM) <- colnames(Left_AM_mean)
AI_AM$distribution<- sample('AM', nrow(AI_AM), replace = T)

AI_RR <- AI_RR %>% 
  as.data.frame() 

AI_RR[,1:3] <- Left_RR_mean[,1:3]

colnames(AI_RR) <- colnames(Left_RR_mean)

AI_RR$distribution<- sample('RR', nrow(AI_RR), replace = T)

AI_2 <- rbind.data.frame(AI_AM, AI_RR) 

AI_2 <- AI_2 %>% 
  as.data.frame()

AI_2 <- AI_2 %>% 
  mutate_if(is.character, as.factor) 

AI_2 <- AI_2[,-3]

head(AI_2)

write.csv(AI_2, paste0(mainDIR, '/AI.csv'))

summarySE(AI_2, measurevar = 'LFF', groupvars = 'distribution')

# testing between population
m1 <- lm(LFT ~ distribution, data = AI_abs)

anova(m1)
summary(m1)

plot(FTL ~ distribution, data = AI_abs)

# Ploting assimetry
ggplot(data = AI_abs, aes(x = distribution, y = LFT, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("LogFA LFF - ", italic("G. underwoodi")))) +
  ggtitle('LFF AI')+
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.75) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.7) +
  stat_boxplot(geom = 'errorbar')+
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  theme_light() +
  theme(legend.position = "none")

# 3. Thermal Preferences >>>> Neonative x Native -------------------------------------

ThermalPref <- read.csv(paste0(mainDIR, "/Underwoodi_thermal.csv"), header = T, sep = ",")

ThermalPref <- ThermalPref %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-(Sp:City), -(Lat:Lon))

ThermalPref <- drop_na(ThermalPref)

head(ThermalPref)
summary(ThermalPref)

g1 <- ggboxplot(ThermalPref, x = "treatment", y = "speed_c_mean", color = "distribution",
                palette = c("#A86D00", "#509AF9"))
g2 <- ggboxplot(ThermalPref, x = "treatment", y = "speed_e_mean", color = "distribution",
                palette = c("#A86D00", "#509AF9"))
g3 <- ggboxplot(ThermalPref, x = "treatment", y = "speed_h_mean", color = "distribution",
                palette = c("#A86D00", "#509AF9"))
g4 <- ggboxplot(Speed, x = "treatment", y = "Speed_mean", color = "distribution",
                palette = c("#A86D00", "#509AF9"))

cowplot::plot_grid(g1, g2, g3, g4,
                   nrow = 2,
                   ncol = 2,
                   labels = "AUTO",
                   label_size = 14,
                   scale = 0.95,
                   align = "hv")

ggplot(data = ThermalPref, aes(x = distribution, y = speed, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("Speed mean (cm/s) -", italic("G. underwoodi")))) +
  ggtitle('Speed temperature')+
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.75) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.7) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  theme_light() +
  theme(legend.position = "none")


## Tpref ---------------------------------------------------------#
# t-test
residual <- lm(tpref ~ distribution, data = ThermalPref)
qqPlot(residual)

shapiro.test(residuals(residual))
leveneTest(tpref ~ distribution, data = ThermalPref)

t.test(tpref ~ distribution, data = ThermalPref, var.equal = T) 

ggplot(data = ThermalPref, aes(x = distribution, y = tpref, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("Temperature Preferential (°C) - ", italic("G. underwoodi")))) +
  ggtitle('Temperature preferential')+
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.6) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.8) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  theme_light() +
  theme(legend.position = "none")

summarySE(ThermalPref, measurevar = 'tpref', groupvars = 'distribution')

## CTmax ---------------------------------------------------------#

head(ThermalPref)

# Nomality
residual <- lm(ctmax ~ distribution, data = ThermalPref)
qqPlot(residual)

# Homocedacidade
shapiro.test(residuals(residual))
leveneTest(ctmax ~ distribution, data = ThermalPref)

Max1.1 <- aov(ctmax ~ distribution, data = ThermalPref)

anova(Max1.1)
summary(Max1.1)

plot_grid(plot_model(Max1.1, type = 'diag')) # normalidade e homogeneidade da variancia

t.test(ctmax ~ distribution, data = ThermalPref, var.equal = T) 


ggplot(data = ThermalPref_2, aes(x = distribution, y = ctmax, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("Critical Temperature Maximum (°C) - ", italic("G. underwoodi")))) +
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.6) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.8) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  theme_light() +
  theme(legend.position = "none")

## CTmin ---------------------------------------------------------#

head(ThermalPref)

# Nomality
residual <- lm(ctmin ~ distribution, data = ThermalPref)
qqPlot(residual)

# Homocedacidade
shapiro.test(residuals(residual))
leveneTest(ctmin ~ distribution, data = ThermalPref)

Min1.1 <- aov(ctmin ~ distribution, data = ThermalPref)

anova(Min1.1)
summary(Min1.1)

plot_grid(plot_model(Min1.1, type = 'diag')) # normalidade e homogeneidade da variancia

t.test(ctmin ~ distribution, data = ThermalPref, var.equal = T)

ggplot(data = ThermalPref_2, aes(x = distribution, y = ctmin, color = distribution)) + 
  labs(x = "distribution", 
       y = expression(paste("Critical Temperature Minimum (°C) - ", italic("G. underwoodi")))) +
  geom_boxplot(fill = c("#A86D00", "#509AF9"), color = "black", 
               outlier.shape = NA, alpha = 0.6) +
  geom_jitter(shape = 16, position = position_jitter(0.1), 
              cex = 5, alpha = 0.8) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  theme_light() +
  theme(legend.position = "none")


# 4. TPC ------------------------------------------------------------------
setwd(mainDIR)
perf_all <- read.csv("TPC.csv", header = T, sep = ",", fileEncoding="latin1")

perf_all <- perf_all %>% 
  mutate_if(is.character, as.factor)

head(perf_all)
summary(perf_all)

### Plot per indv ---------------------------------------------------------------#
ggplot(perf_all, aes(x = temp, y = speed_mean,
                       group = id, colour = distribution)) +
  geom_line() +
  theme_classic()+
  facet_wrap(~ treatment, ncol = 2) + # "curva" por indivíduo 
  theme(legend.title = element_blank()) +
  xlab('Body temperature - ºC') +
  ylab('Speed mean - cm/s')

ggplot(perf_all, aes(x = temp, y = speed_mean,
                       group = id, colour = distribution)) +
  geom_line() +
  theme_classic()+
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  facet_wrap(~ treatment, ncol = 2) + # "curva" por indivíduo 
  theme(legend.title = element_blank()) +
  xlab('Body temperature - ºC') +
  ylab('Speed mean - cm/s')

### GAMM ---------------------------------------------------------------#
### Model construction
if(!require(MuMIn)) install.packages('MuMIn', dependencies = T)

M.1 <- uGamm(speed_mean ~ s(temp, by = distribution)
             + distribution
             + treatment
             + SVL
             + BM
             + TAL
             + HH
             + TRL
             + INF
             + FL_L
             + FL_R
             + TL_L
             + TL_R
             + FTL_L
             + FTL_R
             + FAL_L
             + FAL_R
             + SFT_L
             + SFF_R,
             family = gaussian,
             random = list(id = ~1),
             data = perf_all)

dd2 <- MuMIn::dredge(global.model = M.1, evaluate = T, rank = 'AICc')

#write.csv(dd2, paste0(mainDIR, '/AUCc_TPC_model.csv'))

dd3 <- get.models(dd2, subset = delta < 2) # getting models with delta < 2

get.models(dd2, subset = 1)[[1]] # getting the top model

dd4 <- model.avg(dd3, fit = T) 

summary(dd3)
head(dd2)

best_model <- gamm(speed_mean ~ s(temp, by = distribution, bs = 'cs')
                   + distribution
                   + FAL_L
                   + FTL_L
                   + SVL,
                   random = list(id = ~ 1),
                   data = perf_all)

anova(best_model$gam)
summary(best_model$gam)

### Ploting 
ggTPC <- ggplot(perf_all, aes(x = temp, y = speed_mean, colour = distribution)) +
  geom_point(size = 2, alpha = 0.35) +
  scale_color_manual(values = c("#509AF9", "#A86D00")) +
  scale_fill_manual(values = c("#509AF9", "#A86D00")) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"),  
              aes(fill = distribution), alpha = 0.2, size = 1) +
  scale_y_continuous(limits = c(0,9))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size = 20)) +
  labs(y = "Sprint speed (cm/s)",
       x = "Body temperature (ºC)",
       title = expression(paste("Thermal Performance Curve of ", italic("G. underwoodi"))))+
  annotate('text', label = 'Ctmin', x = 13, y = 0.2, vjust = 0)+
  annotate('text', label = 'Ctmax', x = 47.8, y = 0.2, vjust = 0)+
  scale_x_continuous(breaks = seq(10, 50, 10), limits = c(10,50))


ggplot(perf_all, aes(x = temp, y = speed_mean, colour = distribution)) +
  geom_point(size = 2, alpha = 0.3) +
  scale_color_manual(values = c("#A86D00", "#509AF9")) +
  scale_fill_manual(values = c("#A86D00", "#509AF9")) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"),  
              aes(fill = distribution), alpha = 0.2, size = 1) +
  scale_y_continuous(limits = c(0,9))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size = 20)) +
  facet_wrap(~distribution, ncol = 3) + 
  labs(y = "Sprint speed (cm/s)",
       x = "Body temperature (ºC)",
       title = expression(paste("Thermal Performance Curve of ", italic("G. underwoodi"))))

perf_all_AM <- filter(perf_all, distribution!= 'Neonative')
perf_all_RR <- filter(perf_all, distribution!= 'Native')

ggnative <- ggplot(perf_all_RR, aes(x = temp, y = speed_mean, colour = "#509AF9")) +
  geom_point(size = 2, alpha = 0.35) +
  scale_color_manual(values = c("#509AF9")) +
  scale_fill_manual(values = c("#509AF9")) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"),  
              aes(fill = '#509AF9'), alpha = 0.2, size = 1) +
  scale_y_continuous(limits = c(0,9))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size = 20),legend.position = 'none') +
  labs(y = "Sprint speed (cm/s)",
       x = "Body temperature (ºC)",
       title = expression(paste("Thermal Performance Curve of ", italic("G. underwoodi"), '- RR')))+
  geom_vline(xintercept = c(21.95, 40.31), size = 1, color = 'black', linetype = 'dashed', alpha = .6)+ # B80
  geom_vline(xintercept = c(29.43, 35.08), size = 1.6, color = 'red', linetype = 'dotted', alpha = .8)+ # VT
  geom_vline(xintercept = 32.29, size = 1, color = 'purple', alpha = .5)+
  annotate('text', label = 'Ctmin', x = 15.72, y = 0.5, vjust = 0)+
  annotate('text', label = 'Ctmax', x = 45.59, y = 0.5, vjust = 0)+
  scale_x_continuous(breaks = seq(15, 50, 5))

ggneonative <- ggplot(perf_all_AM, aes(x = temp, y = speed_mean, colour = "#A86D00")) +
  geom_point(size = 2, alpha = 0.35) +
  scale_color_manual(values = c("#A86D00")) +
  scale_fill_manual(values = c("#A86D00")) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"),  
              aes(fill = "#A86D00"), alpha = 0.2, size = 1) +
  scale_y_continuous(limits = c(0,9))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size = 20), legend.position = 'none') +
  labs(y = "Sprint speed (cm/s)",
       x = "Body temperature (ºC)",
       title = expression(paste("Thermal Performance Curve of ", italic("G. underwoodi"), '- AM')))+
  geom_vline(xintercept = c(27.29, 40.58), size = 1, color = 'black', linetype = 'dashed', alpha = .6)+
  geom_vline(xintercept = c(30.171, 36.129), size = 1.6, color = 'red', linetype = 'dotted', alpha = .8)+
  geom_vline(xintercept = 38.21, size = 1, color = 'purple', alpha = .5)+
  annotate('text', label = 'Ctmin', x = 15.2, y = 0.2, vjust = 0)+
  annotate('text', label = 'Ctmax', x = 46.5, y = 0.2, vjust = 0)+
  scale_x_continuous(breaks = seq(15, 50, 10), limits = c(12, 48))

cowplot::plot_grid(ggTPC, ggnative, ggneonative,
                   nrow = 1,
                   ncol = 3,
                   labels = "AUTO",
                   label_size = 14,
                   scale = 0.95,
                   align = "h")

# Splitting data into the AM and RR to calculate the B80% --------------------------------------------# 
perf_all_3 <- drop_na(perf_all)

distribution_data <- perf_all_3[perf_all_3$distribution == unique(perf_all_3$distribution)[2],]

temp_vals <- unique(distribution_data$temp)
SVL_vals <- unique(distribution_data$SVL)
FAL_L_vals <- unique(distribution_data$FAL_L)
FTL_L_vals <- unique(distribution_data$FTL_L)

# Create a new data frame with the predicted values and the corresponding predictor variables
new_data_RR <- expand.grid(temp = seq(min(temp_vals), max(temp_vals), length.out = 100), 
                           distribution = unique(distribution_data$distribution),
                           SVL = seq(min(SVL_vals), max(SVL_vals), length.out = 100),
                           FAL_L = seq(min(FAL_L_vals), max(FAL_L_vals), length.out = 100),
                           FTL_L = seq(min(FTL_L_vals), max(FTL_L_vals), length.out = 100))


new_data_AM <- expand.grid(temp = seq(min(temp_vals), max(temp_vals), length.out = 100), 
                           distribution = unique(distribution_data$distribution),
                           SVL = seq(min(SVL_vals), max(SVL_vals), length.out = 100),
                           FAL_L = seq(min(FAL_L_vals), max(FAL_L_vals), length.out = 100),
                           FTL_L = seq(min(FTL_L_vals), max(FTL_L_vals), length.out = 100))

# Calculate the 80th percentile of the predicted values
new_data_RR$pred <- predict(best_model, newdata = new_data_RR)
new_data_AM$pred <- predict(best_model, newdata = new_data_AM)
summary(new_data)

# Filter the data frame to include only the values that correspond to the 80th percentile
cutoff_value_RR <- quantile(new_data_RR$pred, 0.95)
cutoff_value_AM <- quantile(new_data_AM$pred, 0.95)

# Calculate the minimum and maximum values of "temp"
selected_data_RR <- new_data_RR %>%
  filter(pred >= cutoff_value_RR)

selected_data_AM <- new_data_AM %>%
  filter(pred >= cutoff_value_AM)

min_temp <- min(selected_data$temp) %>% print
max_temp <- max(selected_data$temp) %>% print

summarySE(selected_data_AM, groupvars = 'distribution', measurevar = 'temp')

# b80 neonative 27.29-40.58
# b80 native 21.95-40.31

# distribution     N     temp       sd          se          ci
#       Native 5e+06 32.29798 2.807465  0.001255537 0.002460807
#    Neonative 5e+06 37.61784 0.8144675 0.0003642409 0.0007138993

# VTmin and VTmax --------------------------------------------# 
# Tpref gradient
underwoodiGradient <- read.csv(paste0(mainDIR, "/Tpref_SDM.csv"), header = T, sep = ",", fileEncoding="latin1")
head(underwoodiGradient)

# remove outliers 
outliers_gr <- boxplot(underwoodiGradient$temp)$out
underwoodiGradient_no <- underwoodiGradient[-which(underwoodiGradient$temp %in% outliers_gr),]
boxplot(underwoodiGradient_no$temp)$out

# calcule VTmax (95% percentil) and VTmin (5% percentil)
vtmin_AM <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'AM'], 0.05) %>% print
vtmax_AM <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'AM'], 0.95) %>% print
vtmin_RR <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'RR'], 0.05) %>% print
vtmax_RR <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'RR'], 0.95) %>% print

