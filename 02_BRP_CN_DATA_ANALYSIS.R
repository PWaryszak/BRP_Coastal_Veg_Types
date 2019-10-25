#Plot NEW BRP DATA (See 01_BRP_CN_DATA.R)========
library("tidyverse")#install.packages() first if not in your local library
library(lme4)
library(lmerTest)
library(sjPlot)
library(nlme)
library(sjmisc)

#LOAD and DATA:
setwd("~/00DeakinUni/R/BCL_R/BCL/BRP/BRP_Coastal_Veg_Types")#Set WD = working directory
NewDATA <- read.csv("CN_BRP_CoastalVeg_NewDATA.csv")#Sampel CV-163 was deemd "NA as C.Percent == 0
NewDATA <- NewDATA [ !is.na(NewDATA$site),] #remove NA-s
NewDATA <- NewDATA [ !is.na(NewDATA$C.percent),] #remove NA-s
NewDATA$C.percent <- ifelse(NewDATA$C.percent == 0, 0.001, NewDATA$C.percent)#convert 0 into 0.001 to run log-models
NewDATA$SliceLength.cm <- (NewDATA$DepthTo_cm - NewDATA$DepthFrom_cm) #round % to full numbers to run Poisson
NewDATA$SampleVolume.cm3 <- (pi*(NewDATA$PipeDiameter_cm/2)^2)*NewDATA$SliceLength.cm  #slice volume
NewDATA$CarbonDensity.gcm3 <- NewDATA$dry_bulk_density.gcm3 * NewDATA$C.percent/100
NewDATA$CarbonStock.Mgha <- (((NewDATA$CarbonDensity.gcm3  / 1000000 ) *100000000) * NewDATA$SliceLength.cm )
NewDATA$Core_in.mm <- (NewDATA$PipeLength_cm *10) - NewDATA$Cmptn_in_mm 
NewDATA$Pipe_in.mm <- (NewDATA$PipeLength_cm *10) - NewDATA$Cmptn_out_mm 
NewDATA$Compaction_Correction_Value<- NewDATA$Core_in.mm/NewDATA$Pipe_in.mm
NewDATA$CarbonStock.Mgha_CORRECTED <- NewDATA$CarbonStock.Mgha * NewDATA$Compaction_Correction_Value
NewDATA$DepthCorrected = NewDATA$DepthTo_cm/ NewDATA$Compaction_Correction_Value
NewDATA$DepthChange = NewDATA$DepthCorrected - NewDATA$DepthTo_cm
NewDATA$DepthChangePerc = 100-(NewDATA$DepthChange/NewDATA$DepthTo_cm *100) #Percent increase in depth
NewDATA$CorrectedCarbonStock.Mgha = NewDATA$CarbonStock.Mgha * NewDATA$DepthChangePerc / 100
  
range(NewDATA$Compaction_Correction_Value)# Check if all values are  below 1 = 0.5342903 0.9970760
range(NewDATA$CarbonStock.Mgha )# 0.004913885 68.634732232
range( NewDATA$CarbonStock.Mgha_CORRECTED )# 0.002771935 62.071168539
range(NewDATA$CorrectedCarbonStock.Mgha)


#CORRECTED C-stock to account for compaction:#################

#Analyze  Coastal Veg DATA from Avalon ======
#Gamma Model:
CN_gamma <- glmer(CarbonStock.Mgha_CORRECTED ~ habitat  +DepthRange.cm +
                    (1|core) + (1|site) ,
                  family = Gamma(link = "inverse"), data=NewDATA)
summary(CN_gamma)
plot(resid(CN_gamma))

#Linear Model:
CN_lm <- lm(CarbonStock.Mgha_CORRECTED ~ habitat  +DepthRange.cm , data = NewDATA)
hist(NewDATA$CarbonStock.Mgha_CORRECTED)
plot(resid(CN_lm))

#Linear Model on Site effect:
CN_lm_site <- lm(CarbonStock.Mgha_CORRECTED ~ site , data = NewDATA)
summary(CN_lm_site)
plot(resid(CN_lm_site))

#Linear Model on SiteNumber effect:
CN_lm_siteNum <- lm(CarbonStock.Mgha_CORRECTED ~ as.factor(SiteNumber) , data = NewDATA)
summary(CN_lm_siteNum)
plot(resid(CN_lm_siteNum))

#Log-Linear Model:
CN_log_lm <- lm(log(CarbonStock.Mgha_CORRECTED) ~ habitat + DepthRange.cm , data = NewDATA)
summary(CN_log_lm)
plot(resid(CN_log_lm))

#Gaussian glmm distribution:=
CN_log_lmer1 <- lmer(log(CarbonStock.Mgha_CORRECTED) ~ habitat+ DepthRange.cm   +
                   (1|site)  , data=NewDATA)

summary(CN_log_lmer1)
plot(resid(CN_log_lmer1))

#TABLE:
tab_model(CN_log_lmer1)
tab_model(CN_lm)


#Gaussian glmm distribution no depths:=
CN_log_lmer2 <- lmer(log(CarbonStock.Mgha_CORRECTED) ~ habitat  +
                   (1|site) +(1|DepthRange.cm ),
                 data=NewDATA)

#Compare all models:
AIC(CN_gamma,CN_lm,CN_log_lm, CN_log_lmer1,CN_log_lmer2 )






#MERGE Coastal Veg DATA from Avalon ======
#Merge top and bottom depthRange.cm into two-level factors.
#Rationale:
#--->Top soil contains roots.
#--->Bottom soil shows how well that carbon is beeing stored.
levels(NewDATA$DepthRange)#"00to05" "05to10" "10to20" "20to30"

NewDATA2 <- NewDATA %>%
  mutate(TwoDepths = ifelse(DepthRange.cm == "00to05" | DepthRange.cm== "05to10", "00to10", "10to30"))
#Gamma Model:
CN_gamma2 <- glmer(CarbonStock.Mgha ~ habitat  + TwoDepths +
                     (1|core) + (1|SiteNumber) ,
                   family = Gamma(link = "inverse"), data=NewDATA2)
summary(CN_gamma2)
plot(resid(CN_gamma2))

#Linear Model:
CN_lm2 <- lm(CarbonStock.Mgha ~ habitat  + TwoDepths , data = NewDATA2)
summary(CN_lm2)
plot(resid(CN_lm2))

#Log-Linear Model:
CN_log_lm2 <- lm(log(CarbonStock.Mgha) ~ habitat +  TwoDepths , data = NewDATA2)
summary(CN_log_lm2)
plot(resid(CN_log_lm2))

#Gaussian glmm:
CN_lmer_log2 <- lmer(log(CarbonStock.Mgha) ~ habitat  +  TwoDepths +
                       (1|SiteNumber)++(1|core) , data=NewDATA2)

summary(CN_lmer_log2)
plot(resid(CN_lmer_log2))


#Gaussian glmm non-transfomred:
CN_lmer2 <- lmer(CarbonStock.Mgha ~ habitat  +  TwoDepths +
                   (1|SiteNumber) +(1|core) , data=NewDATA2)

summary(CN_lmer2)
plot(resid(CN_lmer2))

#Compare all models:
AIC(CN_gamma2,CN_lm2, CN_log_lm2, CN_lmer_log2,CN_lmer2 )


#Draw a Stats Table:
tab_model(CN_lmer_log2)
tab_model(CN_gamma)
tab_model(CN_lm2)


#Plot by by habitat
aa <- ggplot(NewDATA, aes(x = habitat, y = CarbonStock.Mgha)) +
  geom_boxplot() +
  facet_grid(.~ DepthRange.cm)+ geom_jitter( alpha = 0.4)+
  ylab("Organic Carbon Stock (Mg/ha)") + xlab("") +
  theme_bw() +
  coord_flip()+
  ggtitle("Coastal Vegetation (BRP Avalon)")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        legend.text = element_text(size = 16),
        strip.text=element_text(size=16),
        plot.title = element_text(size = 20, face = "bold", vjust = 0.5),
        strip.background =  element_rect(fill = "white"))

NewDATA3 <- NewDATA %>%
  group_by(habitat,core) %>%
  summarise(TotalCarbonStock = sum(CarbonStock.Mgha, na.rm = T))
  
a <- ggplot(NewDATA3, aes(x = "Within Each Core", y = TotalCarbonStock)) +
  geom_boxplot() +
  facet_grid(.~ habitat)+ geom_jitter( alpha = 0.4)+
  ylab("Total Carbon Stock (Mg/ha)") + xlab("") +
  theme_bw() +
  ggtitle("Mean total C-stock per core \nCoastal Vegetation (BRP Avalon)")+
  theme(axis.text.x = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 16, angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.position = "none",
        legend.text = element_text(size = 16),
        strip.text=element_text(size=16),
        plot.title = element_text(size = 18, face = "bold", vjust = 0.5),
        strip.background =  element_rect(fill = "white"))
a
library(grid)
library(gridExtra)
grid.arrange(aa, a)

