##################################################################################
##################################################################################
##################################################################################
# ------------------------------------------------------------------------------ #
# ////// Guinea grass Allelopathy-Competition Experiment ANALYSES (G039)# ////// #
# ------------------------------------------------------------------------------ #
##################################################################################
##################################################################################
##################################################################################

# Colin Richard Morrison 
# PhD Candidate
# The University of Texas at Austin 
# Department of Integrative Biology 
# Graduate Program In Ecology, Evolution and Behavior
# crmorrison@utexas.edu

getwd()
setwd("~/Desktop/guinea.grass.chemistry/allelopathy/analyses")

library(ggplot2)
library(plotrix)
library(multcomp)
library(lme4)
library(plyr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(survival)
library(nlme) ### Mixed effects modelpackage
library(MASS) #for glmmPQL
#library(psych)
#library(MuMIn)  ### This has the r.squaredGLMM command from
#library(DescTools)


### make summary data table with mean, sd, and se for the plot
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#####################################################################################
#####################################################################################
#### --- 2-Hydroxyphenylacetic Acid Concentrations from Rancho La Paloma, TX --- ####
#####################################################################################
#####################################################################################
### GOAL:
# 1.) Analyze differences in 2-HPAA concentration between Guinea grass tissues and soil

dats<-read.csv("LP.2HPAA.data.csv")
lp.chem<-dats[,c(3,6:7)]


### ------------------------------------------------------------ ###
### /// Analyze 2-HPAA concentrations per type data with GLM /// ###
### ------------------------------------------------------------ ###
### test for normality 
hist(dats$pg.mg)
qqnorm(dats$pg.mg) # huge right skew
### do a log10 transformation 
dats$log.pg.mg<-log(dats$pg.mg + 1)
hist(dats$log.pg.mg) # distribution more normal now
qqnorm(dats$log.pg.mg)  
shapiro.test(dats$log.pg.mg) # p-value = 0.1677
boxplot(log.pg.mg~type, data = dats) # variance more homogeneous now

### --- analyze with a GLM
dats$type=as.factor(dats$type)

m=glm(log.pg.mg ~ type,
      data = dats)
anova(m,test="F")
#      Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                    46     29.173                     
# type  3   13.585        43     15.588 12.492 5.232e-06 ***

summary(m)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.3091     0.1670   7.839 8.01e-10 ***
# typelitter    0.3996     0.2410   1.658  0.10462    
# typeroot      0.7210     0.2362   3.053  0.00388 ** 
# typesoil     -0.8108     0.2611  -3.105  0.00336 ** 

### post-hoc test
summary(glht(m, linfct=mcp(type="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|)    
#litter - leaf == 0   0.3996     0.2410   1.658   0.3456    
# root - leaf == 0     0.7210     0.2362   3.053   0.0119 *  
# soil - leaf == 0    -0.8108     0.2611  -3.105   0.0105 *  
# root - litter == 0   0.3214     0.2410   1.333   0.5409    
# soil - litter == 0  -1.2104     0.2655  -4.559   <0.001 ***
# soil - root == 0    -1.5317     0.2611  -5.867   <0.001 *** 


### --------------------------------------------------------- ###
### /// Calculate Gini Coefficient for 2HPAA across space /// ###
### --------------------------------------------------------- ###
# ranges from 0 to 1, 
# with 1 being the highest possible inequality,high concentrations all at one site
# 0 being lowest, concentrations homogenous across space
?Gini

dats<-read.csv('LP_2HPAA_correlates.csv')
#log.pg.mg<-log(dats$pg.mg)
log.pg.mg<-log(dats$pg.mg+1) # add a constant to deal with negative values
dats=cbind.data.frame(dats,log.pg.mg)
dats

# subset dataframe for each substrate
names(dats)
leaf<-dats[which(dats$type=='leaf'),names(dats) %in%
             c("site","longitude","latitude","type","pg.mg","log.pg.mg")]
litter<-dats[which(dats$type=='litter'),names(dats) %in%
               c("site","longitude","latitude","type","pg.mg","log.pg.mg")]
root<-dats[which(dats$type=='root'),names(dats) %in%
             c("site","longitude","latitude","type","pg.mg","log.pg.mg")]
soil<-dats[which(dats$type=='soil'),names(dats) %in%
             c("site","longitude","latitude","type","pg.mg","log.pg.mg")]

# 2HPAA ranges (pg/mg)
range(soil$pg.mg) # 0.3235263 1.9385575
range(leaf$pg.mg) # 0.7267854 12.8231068
range(litter$pg.mg) # 1.749367 31.009361
range(root$pg.mg) # 2.830972 61.956404

### soil Gini
hist(soil$pg.mg)
Gini(soil$pg.mg,n=nrow(soil),R = 1000, type = "bca",conf.level=0.95,unbiased=TRUE) 
#      gini    lwr.ci    upr.ci 
# 0.2852814 0.2233629 0.3462005

### roots Gini
hist(root$pg.mg)
Gini(root$pg.mg,n=nrow(root),R = 1000, type = "bca",conf.level=0.95,unbiased=TRUE) 
# gini    lwr.ci    upr.ci 
# 0.5216274 0.4494924 0.5725235

### leaves Gini
hist(leaf$pg.mg)
Gini(leaf$pg.mg,n=nrow(leaf),R = 1000, type = "bca",conf.level=0.95,unbiased=TRUE) 
#      gini    lwr.ci    upr.ci 
# 0.4181239 0.3787839 0.4588521 

### litter Gini
hist(litter$pg.mg)
Gini(litter$pg.mg,n=nrow(litter),R = 1000, type = "bca",conf.level=0.95,unbiased=TRUE)
#      gini    lwr.ci    upr.ci 
# 0.4664463 0.4029593 0.5206388



####################################################################
####################################################################
### ---  Seedling Survival Across Chemical and Light Levels  --- ###
####################################################################
####################################################################
surv<-read.csv("survival2.csv")
head(surv)

### ------------------------------ ###
### /// Final Percent Survival /// ###
### ------------------------------ ###
### GOAL:
# 1.) Analyze effects of light, chemical and interaction on percent survival with negative binomial GLM.
#   - Analyze pairwise survival across chemical treatments with post-hoc tests for each species.
final = filter(surv,week==5)
nrow(final) # [1] 2352
head(final)

### summary recruitment statistics
stats <- data_summary(final, varname="result", 
                      groupnames=c("chemical","light"))
stats

##########################
### /// M. maximus /// ###
##########################
surv_MM = final[(final$species =="M.maximus"),] #retaining only M.maximus
dotchart(surv_MM$result, xlab = "survival",
         ylab = "Order of the data")
### check distribution
Lm=lm(result~chemical*light,data=surv_MM)
plot(hist(resid(Lm))) #Plotting a histogram of the model residuals shows normally distributed data.

### --- GLM fit with logit link style negative binomial distributions
surv_MM$chemical<-as.factor(surv_MM$chemical)
surv_MM$light<-as.factor(surv_MM$light)

### check homogeneity of variance
Lm1=glm.nb(result~chemical*light, link=log,data=surv_MM)
summary(Lm1) 
anova(Lm1,test="F") # the interaction is significant so keep both terms. Report p and F
#               Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                             587     413.89                      
# chemical        2  27.4068       585     386.49 13.7034 1.119e-06 ***
# light           1  22.3480       584     364.14 22.3480 2.274e-06 ***
# chemical:light  2   7.0625       582     357.08  3.5312   0.02927 * 

### --- post-hoc tests for pairwise comparisons
### sun plants
surv_MM_sun = surv_MM[(surv_MM$light =="open"),] 
surv_MM_sun$chemical<-as.factor(surv_MM_sun$chemical)
Lm1sun=glm(result~chemical,data=surv_MM_sun)
summary(Lm1sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      0.25510    0.06659   3.831 0.000363 ***
# standard - crude == 0 -0.09184    0.06659  -1.379 0.351895    
# standard - none == 0  -0.34694    0.06659  -5.210  < 1e-04 ***
### shade plants
surv_MM_shade = surv_MM[(surv_MM$light =="shade"),] 
Lm1shade=glm(result~chemical,data=surv_MM_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      0.40816    0.06372   6.405   <1e-04 ***
# standard - crude == 0  0.11224    0.06372   1.761    0.183    
# standard - none == 0  -0.29592    0.06372  -4.644   <1e-04 *** 

### no chemical plants
surv_MM_no = surv_MM[(surv_MM$chemical =="none"),]
none=glm(result~light,data=surv_MM_no)
summary(none)
#            Estimate Std. Error t value Pr(>|t|)    
# lightshade  -0.23469    0.06367  -3.686 0.000296 ***

### 2HPAA standard plants
surv_MM_stand = surv_MM[(surv_MM$chemical =="standard"),]
stand=glm(result~light,data=surv_MM_stand)
summary(stand)
#            Estimate Std. Error t value Pr(>|t|)    
# lightshade  -0.18367    0.06835  -2.687  0.00783 ** 

### whole plant extract plants
surv_MM_crude = surv_MM[(surv_MM$chemical =="crude"),]
crude=glm(result~light,data=surv_MM_crude)
summary(crude)
#            Estimate Std. Error t value Pr(>|t|)    
# lightshade  -0.38776    0.06338  -6.118 5.13e-09 ***

###############################
### /// C. fasciculata /// ###
##############################
surv_CF = final[(final$species =="C.fasciculata"),] 
dotchart(surv_CF$result, xlab = "survival",
         ylab = "Order of the data")
### check distribution
Lm=lm(result~chemical*light,data=surv_CF)
plot(hist(resid(Lm))) #Plotting a histogram of the model residuals shows normally distributed data.

### --- GLM fit with logit link style negative binomial distributions
surv_CF$chemical<-as.factor(surv_CF$chemical)
surv_CF$light<-as.factor(surv_CF$light)

Lm1=glm.nb(result~chemical*light,link='log',data=surv_CF)
summary(Lm1) 
anova(Lm1,test="F") # the interaction is significant so keep both terms. Report p and F
#               Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                             587     430.87                      
# chemical        2   23.137       585     407.73 11.5684 9.460e-06 ***
# light           1   33.618       584     374.12 33.6181 6.707e-09 ***
# chemical:light  2    6.393       582     367.72  3.1965   0.04091 * 

### --- post-hoc tests for pairwise comparisons
### sun plants
surv_CF_sun = surv_CF[(surv_CF$light =="open"),] 
surv_CF_sun$chemical<-as.factor(surv_CF_sun$chemical)
Lm1sun=glm(result~chemical,data=surv_CF_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      0.12245    0.06576   1.862     0.15    
# standard - crude == 0  0.45918    0.06576   6.983   <1e-04 ***
# standard - none == 0   0.33673    0.06576   5.121   <1e-04 ***
### shade plants
surv_CF_shade = surv_CF[(surv_CF$light =="shade"),] 
Lm1shade=glm(result~chemical,data=surv_CF_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      0.21429    0.06099   3.513  0.00133 **
# standard - crude == 0  0.14286    0.06099   2.342  0.05023 . 
# standard - none == 0  -0.07143    0.06099  -1.171  0.47046   

### no chemical plants
surv_CF_no = surv_CF[(surv_CF$chemical =="none"),]
none=glm(result~light,data=surv_CF_no)
summary(none)
#             Estimate Std. Error t value Pr(>|t|)    
# lightshade  -0.13265    0.07006  -1.893   0.0548 *

### 2HPAA standard plants
surv_CF_stand = surv_CF[(surv_CF$chemical =="standard"),]
stand=glm(result~light,data=surv_CF_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -0.54082    0.06003  -9.009   <2e-16 ***

### whole plant extract plants
surv_CF_crude = surv_CF[(surv_CF$chemical =="crude"),]
crude=glm(result~light,data=surv_CF_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -0.22449    0.05961  -3.766  0.00022 ***


###############################
### /// S. scheelei /// ###
##############################
surv_SS = final[(final$species =="S.scheelei"),] 
dotchart(surv_SS$result, xlab = "survival",
         ylab = "Order of the data")
### check distribution
Lm=lm(result~chemical*light,data=surv_SS)
plot(hist(resid(Lm))) 

### --- GLM fit with logit link style negative binomial distributions
surv_SS$chemical<-as.factor(surv_SS$chemical)
surv_SS$light<-as.factor(surv_SS$light)
Lm1=glm.nb(result~chemical*light,link=log,data=surv_SS)
summary(Lm1) 
anova(Lm1,test="F") # the interaction is non-significant. 
#               Df Deviance Resid. Df Resid. Dev      F  Pr(>F)  
# NULL                             587     399.72                 
# chemical        2   4.9318       585     394.78 2.4659 0.08493 .
# light           1   0.3268       584     394.46 0.3268 0.56752  
# chemical:light  2   0.6719       582     393.78 0.3359 0.71467 

### --- post-hoc tests for pairwise comparisons
### sun plants
surv_SS_sun = surv_SS[(surv_SS$light =="open"),] 
surv_SS_sun$chemical<-as.factor(surv_SS_sun$chemical)
Lm1sun=glm(result~chemical,data=surv_SS_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      0.05102    0.07123   0.716    0.754
# standard - crude == 0  0.12245    0.07123   1.719    0.198
# standard - none == 0   0.07143    0.07123   1.003    0.575
### shade plants
surv_SS_shade = surv_SS[(surv_SS$light =="shade"),] 
Lm1shade=glm(result~chemical,data=surv_SS_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      0.15306    0.07079   2.162   0.0778 .
# standard - crude == 0  0.19388    0.07079   2.739   0.0170 *
# standard - none == 0   0.04082    0.07079   0.577   0.8326  

### no chemical plants
surv_SS_no = surv_SS[(surv_SS$chemical =="none"),]
none=glm(result~light,data=surv_SS_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade   0.01020    0.07161   0.142    0.887   

### 2HPAA standard plants
surv_SS_stand = surv_SS[(surv_SS$chemical =="standard"),]
stand=glm(result~light,data=surv_SS_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -0.02041    0.07056  -0.289    0.773  
### whole plant extract plants
surv_SS_crude = surv_SS[(surv_SS$chemical =="crude"),]
crude=glm(result~light,data=surv_SS_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# shade - open == 0 -0.09184    0.07086  -1.296    0.196


##########################
### /// X. texanum /// ###
##########################
surv2<-read.csv("survival3.csv") # 0 values for crude omited for model comparison (avoid overdispersion)
final2 = filter(surv2,week==5)
surv_XT = final[(final$species =="X.texanum"),] 
dotchart(surv_XT$result, xlab = "survival",
         ylab = "Order of the data")
### check distribution
Lm=lm(result~chemical*light,data=surv_XT)
plot(hist(resid(Lm))) #Plotting a histogram of the model residuals shows normally distributed data.

### --- GLM fit with logit link style negative binomial distributions
surv_XT$chemical<-as.factor(surv_XT$chemical)
surv_XT$light<-as.factor(surv_XT$light)
#Lm0=glm(result~1,data=surv_XT)
#Lm1=glm.nb(result~chemical*light,link=log,data=surv_XT)
#Lm2=glm(result~chemical*light,family = poisson,data=surv_XT)
Lm3=glm(result~chemical*light,family = binomial,data=surv_XT)
#AICc(Lm0,Lm1,Lm2,Lm3)
summary(Lm3)
anova(Lm3, test = "F")
#                Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                             587     485.87                      
# chemical        2   48.160       585     437.71 24.0800 3.485e-11 ***
# light           1   16.356       584     421.35 16.3563 5.248e-05 ***
# chemical:light  2    0.468       582     420.88  0.2341    0.7913      


### --- post-hoc tests for pairwise comparisons
### sun plants
surv_XT_sun = surv_XT[(surv_XT$light =="open"),] 
Lm1sun=glm(result~chemical,data=surv_XT_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      0.29592    0.05295   5.589   <1e-05 ***
# standard - crude == 0  0.27551    0.05295   5.203   <1e-05 ***
# standard - none == 0  -0.02041    0.05295  -0.385    0.921 
### shade plants
surv_XT_shade = surv_XT[(surv_XT$light =="shade"),] 
Lm1shade=glm(result~chemical,data=surv_XT_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      0.15306    0.03900   3.925   <0.001 ***
# standard - crude == 0  0.10204    0.03900   2.617   0.0241 *  
# standard - none == 0  -0.05102    0.03900  -1.308   0.3905 

### no chemical plants
surv_XT_no = surv_XT[(surv_XT$chemical =="none"),]
none=glm(result~light,data=surv_XT_no)
summary(none)
#                   Estimate Std. Error z value Pr(>|z|)  
# lightshade  -0.14286    0.05903   -2.42   0.0164 *  
### 2HPAA standard plants
surv_XT_stand = surv_XT[(surv_XT$chemical =="standard"),]
stand=glm(result~light,data=surv_XT_stand)
summary(stand)
#                   Estimate Std. Error z value Pr(>|z|) 
# lightshade  -0.17347    0.05479  -3.166   0.0018 ** 
### whole plant extract plants
surv_XT_crude = surv_XT[(surv_XT$chemical =="crude"),]
crude=glm(result~light,data=surv_XT_crude)
summary(crude)
# NA all plants died 



### ---------------------------------- ###
### /// Percent Survival over Time /// ###
### ---------------------------------- ###
### GOAL:
# 1.) Analyze with survival analyses to see if survival curves diverge over time by LIGHT level.
#   - log rank sum tests
# 2.) Analyze with survival analyses to see if survival curves diverge over time by CHEMICAL level.
#   - log rank sum tests

surv<-read.csv("survival.csv")
surv=surv[,-c(4)]
head(surv)

### use melted dataframe with rows of individual plants that survived or died (= 1 or 0)
surv.m <- melt(surv, id = c("species","chemical","shade","week"))
surv.m=surv.m[,-c(5)]
head(surv.m,50)
names(surv.m)<-c("species","chemical","light","week","result")
write.csv(surv.m,file="survival2.csv")

### --- 1.) First Analysis is across chemical level effects across light levels --- ###
# significant differences and summary stats will be shown in a Table

# subset for chemical treatments 
surv.cont<-filter(surv.m, chemical=='none')
surv.stand<-filter(surv.m, chemical=='standard')
surv.crude<-filter(surv.m, chemical=='crude')

### - CONTROL treatments
### M .maximus 
Mmax.cont<-filter(surv.cont, species=='M.maximus')
Mmax.cont.fit <- survdiff(Surv(week, result) ~ light, data = Mmax.cont)
# Chisq= 9.3  on 1 degrees of freedom, p= 0.002 
nrow(Mmax.cont)
### C. fasciculata 
Cfas.cont<-filter(surv.cont, species=='C.fasciculata')
Cfas.cont.fit <- survdiff(Surv(week, result) ~ light, data = Cfas.cont)
# Chisq= 2.7  on 1 degrees of freedom, p= 0.1
### S. scheelei 
Ssch.cont<-filter(surv.cont, species=='S.scheelei')
Ssch.cont.fit <- survdiff(Surv(week, result) ~ light, data = Ssch.cont)
# Chisq= 0  on 1 degrees of freedom, p= 1 
### X. texanum
Xtex.cont<-filter(surv.cont, species=='X.texanum')
Xtex.cont.fit <- survdiff(Surv(week, result) ~ light, data = Xtex.cont)
# Chisq= 2.8  on 1 degrees of freedom, p= 0.1 

### - STANDARD extract treatments
### M .maximus 
Mmax.stand<-filter(surv.stand, species=='M.maximus')
Mmax.stand.fit <- survdiff(Surv(week, result) ~ light, data = Mmax.stand)
# Chisq= 11.1  on 1 degrees of freedom, p= 8e-04 
# C. fasciculata 
Cfas.stand<-filter(surv.stand, species=='C.fasciculata')
Cfas.stand.fit <- survdiff(Surv(week, result) ~ light, data = Cfas.stand)
# Chisq= 71.6  on 1 degrees of freedom, p= <2e-16
### S. scheelei 
Ssch.stand<-filter(surv.stand, species=='S.scheelei')
Ssch.stand.fit <- survdiff(Surv(week, result) ~ light, data = Ssch.stand)
# Chisq= 0.1  on 1 degrees of freedom, p= 0.7 
### X. texanum
Xtex.stand<-filter(surv.stand, species=='X.texanum')
Xtex.stand.fit <- survdiff(Surv(week, result) ~ light, data = Xtex.stand)
# Chisq= 12.8  on 1 degrees of freedom, p= 3e-04

### - Whole-plant extract treatment
### M .maximus 
Mmax.crude<-filter(surv.crude, species=='M.maximus')
Mmax.crude.fit <- survdiff(Surv(week, result) ~ light, data = Mmax.crude)
# Chisq= 9.1  on 1 degrees of freedom, p= 0.003
# C. fasciculata 
Cfas.crude<-filter(surv.crude, species=='C.fasciculata')
Cfas.crude.fit <- survdiff(Surv(week, result) ~ light, data = Cfas.crude)
# Chisq= 0.4  on 1 degrees of freedom, p= 0.5 
### S. scheelei 
Ssch.crude<-filter(surv.crude, species=='S.scheelei')
Ssch.crude.fit <- survdiff(Surv(week, result) ~ light, data = Ssch.crude)
# Chisq= 13.8  on 1 degrees of freedom, p= 2e-04 
### X. texanum
Xtex.crude<-filter(surv.crude, species=='X.texanum')
Xtex.crude.fit <- survdiff(Surv(week, result) ~ light, data = Xtex.crude)
#

### --- 2.) Second Analysis is among chemical level effects within light level --- ###
# significant differences will be shown on the Figure with compact letter display

# subset for light treatments 
surv.shade<-filter(surv.m, light=='shade')
surv.sun<-filter(surv.m, light=='open')

### --- subset and analyze for species in SUN
### M. maximus pairwise comparisons (SUN)
Mmax.sun<-filter(surv.sun, species=='M.maximus')
Mm.sun.con.sta<-filter(Mmax.sun, chemical=='none' | chemical=='standard')
Mm.sun.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Mm.sun.con.sta)
# Chisq= 75.6  on 1 degrees of freedom, p= <2e-16
Mm.sun.con.cru<-filter(Mmax.sun, chemical=='none' | chemical=='crude')
Mm.sun.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Mm.sun.con.cru)
# Chisq= 86.8  on 1 degrees of freedom, p= <2e-16 
Mm.sun.sta.cru<-filter(Mmax.sun, chemical=='standard' | chemical=='crude')
Mm.sun.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Mm.sun.sta.cru)
# Chisq= 0.1  on 1 degrees of freedom, p= 0.7 

### C. fasciculata pairwise comparisons (SUN)
Cfas.sun<-filter(surv.sun, species=='C.fasciculata')
Cf.sun.con.sta<-filter(Cfas.sun, chemical=='none' | chemical=='standard')
Cf.sun.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Cf.sun.con.sta)
# Chisq= 48.7  on 1 degrees of freedom, p= 3e-12 
Cf.sun.con.cru<-filter(Cfas.sun, chemical=='none' | chemical=='crude')
Cf.sun.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Cf.sun.con.cru)
# Chisq= 6.6  on 1 degrees of freedom, p= 0.01
Cf.sun.sta.cru<-filter(Cfas.sun, chemical=='standard' | chemical=='crude')
Cf.sun.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Cf.sun.sta.cru)
# Chisq= 91.2  on 1 degrees of freedom, p= <2e-16

### S.scheelei pairwise comparisons (SUN)
Ssch.sun<-filter(surv.sun, species=='S.scheelei')
Ss.sun.con.sta<-filter(Ssch.sun, chemical=='none' | chemical=='standard')
Ss.sun.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Ss.sun.con.sta)
#  Chisq= 0.1  on 1 degrees of freedom, p= 0.7 
Ss.sun.con.cru<-filter(Ssch.sun, chemical=='none' | chemical=='crude')
Ss.sun.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Ss.sun.con.cru)
# Chisq= 1.4  on 1 degrees of freedom, p= 0.2 
Ss.sun.sta.cru<-filter(Ssch.sun, chemical=='standard' | chemical=='crude')
Ss.sun.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Ss.sun.sta.cru)
# Chisq= 0.7  on 1 degrees of freedom, p= 0.4 

### X. texanum pairwise comparisons (SUN)
Xtex.sun<-filter(surv.sun, species=='X.texanum')
Xt.sun.con.sta<-filter(Xtex.sun, chemical=='none' | chemical=='standard')
Xt.sun.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Xt.sun.con.sta)
# Chisq= 0  on 1 degrees of freedom, p= 0.9 
Xt.sun.con.cru<-filter(Xtex.sun, chemical=='none' | chemical=='crude')
Xr.sun.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Xt.sun.con.cru)
# Chisq= 82.5  on 1 degrees of freedom, p= <2e-16 
Xt.sun.sta.cru<-filter(Xtex.sun, chemical=='standard' | chemical=='crude')
Xr.sun.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Xt.sun.sta.cru)
# Chisq= 80  on 1 degrees of freedom, p= <2e-16 


### --- subset and analyze for species in SHADE
### M. maximus pairwise comparisons (SHADE)
Mmax.shade<-filter(surv.shade, species=='M.maximus')
Mm.shade.con.sta<-filter(Mmax.shade, chemical=='none' | chemical=='standard')
Mm.shade.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Mm.shade.con.sta)
# Chisq= 76.3  on 1 degrees of freedom, p= <2e-16 
Mm.shade.con.cru<-filter(Mmax.shade, chemical=='none' | chemical=='crude')
Mm.shade.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Mm.shade.con.cru)
# Chisq= 74.9  on 1 degrees of freedom, p= <2e-16 
Mm.shade.sta.cru<-filter(Mmax.shade, chemical=='standard' | chemical=='crude')
Mm.shade.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Mm.shade.sta.cru)
# Chisq= 0  on 1 degrees of freedom, p= 1 

### C. fasciculata pairwise comparisons (SHADE)
Cfas.shade<-filter(surv.shade, species=='C.fasciculata')
Cf.shade.con.sta<-filter(Cfas.shade, chemical=='none' | chemical=='standard')
Cf.shade.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Cf.shade.con.sta)
# Chisq= 0  on 1 degrees of freedom, p= 1 
Cf.shade.con.cru<-filter(Cfas.shade, chemical=='none' | chemical=='crude')
Cf.shade.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Cf.shade.con.cru)
# Chisq= 2.3  on 1 degrees of freedom, p= 0.1 
Cf.shade.sta.cru<-filter(Cfas.shade, chemical=='standard' | chemical=='crude')
Cf.shade.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Cf.shade.sta.cru)
# Chisq= 2.3  on 1 degrees of freedom, p= 0.1 

### S.scheelei pairwise comparisons (SHADE)
Ssch.shade<-filter(surv.shade, species=='S.scheelei')
Ss.shade.con.sta<-filter(Ssch.shade, chemical=='none' | chemical=='standard')
Ss.shade.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Ss.shade.con.sta)
#  Chisq= 0.5  on 1 degrees of freedom, p= 0.5 
Ss.shade.con.cru<-filter(Ssch.shade, chemical=='none' | chemical=='crude')
Ss.shade.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Ss.shade.con.cru)
#  Chisq= 7.1  on 1 degrees of freedom, p= 0.008 
Ss.shade.sta.cru<-filter(Ssch.shade, chemical=='standard' | chemical=='crude')
Ss.shade.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Ss.shade.sta.cru)
# Chisq= 11.3  on 1 degrees of freedom, p= 8e-04

### X. texanum pairwise comparisons (SHADE)
Xtex.shade<-filter(surv.shade, species=='X.texanum')
Xt.shade.con.sta<-filter(Xtex.shade, chemical=='none' | chemical=='standard')
Xt.shade.con.sta.fit <- survdiff(Surv(week, result) ~ chemical, data = Xt.shade.con.sta)
#  Chisq= 4.5  on 1 degrees of freedom, p= 0.03 
Xt.shade.con.cru<-filter(Xtex.shade, chemical=='none' | chemical=='crude')
Xr.shade.con.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Xt.shade.con.cru)
# Chisq= 54.8  on 1 degrees of freedom, p= 1e-13 
Xt.shade.sta.cru<-filter(Xtex.shade, chemical=='standard' | chemical=='crude')
Xr.shade.sta.cru.fit <- survdiff(Surv(week, result) ~ chemical, data = Xt.shade.sta.cru)
# Chisq= 33.2  on 1 degrees of freedom, p= 8e-09 



###########################################################################
###########################################################################
### --- Relative Shoot Height Across Common Garden Treatment Levels --- ###
###########################################################################
###########################################################################
### GOALS:
# 1.) Analyze effects of light, chemical and interaction on shoot height with gaussian GLS.
# 2.) Analyze pairwise height differences across shade treatments with post-hoc tests for each species.

### subset just the final measurements
grow<-read.csv("growth.csv")
#head(grow,10)
grow.final = filter(grow,week==4)
grow2 = grow.final[,-c(4:7)]
### melt the dataframe so that each growth value is a row
grow.m <- melt(grow2, id = c("species","chemical","shade"))
grow.m=na.omit(grow.m[,-c(4)])
grow.m$value<-as.numeric(grow.m$value)

grow.m$log.value=log(grow.m$value +1)
names(grow.m)<-c("species","chemical","light","height","log.height")
grow.m<-na.omit(grow.m)
write.csv(grow.m, file="growth3.csv")


##########################
### /// M. maximus /// ###
##########################
grow<-read.csv("growth3.csv")
head(grow,10)

grow_MM = grow[(grow$species =="M.maximus"),] #retaining only M.maximus
dotchart(grow_MM$height, xlab = "growth",
         ylab = "Order of the data")
### check distribution
Lm=lm(height~chemical*light,data=grow_MM)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals shows normally distributed data.

### check homogeneity of variance
boxplot(grow_MM$height~chemical, data = grow_MM) # Low heterogeneity of variance
boxplot(grow_MM$height~light, data = grow_MM) # Some heterogeneity of variance

### --- GLS fit with gaussian distribution
grow_MM$chemical<-as.factor(grow_MM$chemical)
grow_MM$light<-as.factor(grow_MM$light)

Lm1=gls(height~chemical*light,data=grow_MM,weights=varIdent(form=~1|as.factor(light)))
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                  numDF   F-value p-value 
#   (Intercept)        1 4185.462  <.0001
# chemical           2    2.0543  0.1301
# light              1 1140.3692  <.0001
# chemical:light     2    5.4799  0.0046

### --- post-hoc tests for pairwise comparisons
### sun plants
MM_sun = grow_MM[(grow_MM$light =="open"),] 
Lm1sun=glm(height~chemical,data=MM_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0        7.250      1.101   6.587  < 1e-04 ***
# standard - crude == 0    4.865      1.256   3.875 0.000309 ***
# standard - none == 0    -2.385      1.163  -2.051 0.099743 .
length(which(MM_sun$chemical =='none')) # 80
length(which(MM_sun$chemical =='standard')) # 46 
length(which(MM_sun$chemical =='crude')) # 55

### shade plants
MM_shade = grow_MM[(grow_MM$light =="shade"),] 
Lm1shade=glm(height~chemical,data=MM_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      2.29876    1.12363   2.046   0.0396 .
# standard - crude == 0  0.03824    1.25018   0.031   0.9995  
# standard - none == 0  -2.26053    0.93835  -2.409   0.0414 *
length(which(MM_shade$chemical =='none')) # 57
length(which(MM_shade$chemical =='standard')) # 28
length(which(MM_shade$chemical =='crude')) # 17

### no chemical plants
MM_no = grow_MM[(grow_MM$chemical =="none"),]
none=glm(height~light,data=MM_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade  -22.2945     0.7747  -28.78   <2e-16 ***

### 2HPAA standard plants
MM_stand = grow_MM[(grow_MM$chemical =="standard"),]
stand=glm(height~light,data=MM_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -22.170      1.392  -15.92   <2e-16 ***

### whole plant extract plants
MM_crude = grow_MM[(grow_MM$chemical =="crude"),]
crude=glm(height~light,data=MM_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -17.343      1.973  -8.789   <2e-16 ***


##############################
### /// C. fasciculata /// ###
##############################
grow_Cf = grow[(grow$species =="C.fasciculata"),] 
dotchart(grow_Cf$height, xlab = "growth",
         ylab = "Order of the data") # no outliers
### check distribution
Lm=lm(height~chemical*light,data=grow_Cf)
plot(hist(resid(Lm))) #  normally distributed data.

### check homogeneity of variance
boxplot(grow_Cf$height~chemical, data = grow_Cf) # Low heterogeneity of variance
boxplot(grow_Cf$height~light, data = grow_Cf) # Some heterogeneity of variance

### --- GLS fit with gaussian distribution
grow_Cf$chemical<-as.factor(grow_Cf$chemical)
grow_Cf$light<-as.factor(grow_Cf$light)

Lm1=gls(height~chemical*light,data=grow_Cf,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                Df Deviance Resid. Df Resid. Dev  Pr(>Chi) 
# chemical           2   25.6707  <.0001
# light              1   62.0453  <.0001
# chemical:light     2   14.9577  <.0001

### --- post-hoc tests for pairwise comparisons
### sun plants
CF_sun = grow_Cf[(grow_Cf$light =="open"),] 
CF_sun$chemical<-as.factor(CF_sun$chemical)
Lm1sun=glm(height~chemical,data=CF_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      -0.3038     0.8310  -0.366    0.929    
# standard - crude == 0   4.0011     0.7543   5.304   <1e-05 ***
# standard - none == 0    4.3048     0.6841   6.293   <1e-05 ***
length(which(CF_sun$chemical =='none')) # 47
length(which(CF_sun$chemical =='standard')) # 80 
length(which(CF_sun$chemical =='crude')) # 35

### shade plants
CF_shade = grow_Cf[(grow_Cf$light =="shade"),] 
Lm1shade=glm(height~chemical,data=CF_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0       3.6855     0.7700   4.787  < 1e-04 ***
# standard - crude == 0   3.1274     0.7971   3.923 0.000265 ***
# standard - none == 0   -0.5582     0.6087  -0.917 0.626776 
length(which(CF_shade$chemical =='none')) # 34
length(which(CF_shade$chemical =='standard')) # 27
length(which(CF_shade$chemical =='crude')) # 13

### no chemical plants
CF_no = grow_Cf[(grow_Cf$chemical =="none"),]
none=glm(height~light,data=CF_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade   -0.2806     0.7264  -0.386    0.699

### 2HPAA standard plants
CF_stand = grow_Cf[(grow_Cf$chemical =="standard"),]
stand=glm(height~light,data=CF_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -5.1436     0.7895  -6.515 2.57e-09 ***

### whole plant extract plants
CF_crude = grow_Cf[(grow_Cf$chemical =="crude"),]
crude=glm(height~light,data=CF_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -4.2699     1.0218  -4.179  0.00013 ***


###########################
### /// S. scheelei /// ###
###########################
grow_Ss = grow[(grow$species =="S.scheelei"),] 
dotchart(grow_Ss$height, xlab = "growth",
         ylab = "Order of the data") # no outliers
### check distribution
Lm=lm(height~chemical*light,data=grow_Ss)
plot(hist(resid(Lm))) #  normally distributed data.

boxplot(grow_Ss$height~chemical, data = grow_Ss) # Low heterogeneity of variance
boxplot(grow_Ss$height~light, data = grow_Ss) # Some heterogeneity of variance

### --- GLS fit with gaussian distribution
grow_Ss$chemical<-as.factor(grow_Ss$chemical)
grow_Ss$light<-as.factor(grow_Ss$light)

Lm1=gls(height~chemical*light,data=grow_Ss,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) 
#                Df Deviance Resid. Df Resid. Dev  Pr(>Chi) 
# chemical           2   17.6931  <.0001
# light              1   74.0791  <.0001
# chemical:light     2    1.6573  0.1924

### --- post-hoc tests for pairwise comparisons
### sun plants
SS_sun = grow_Ss[(grow_Ss$light =="open"),] 
SS_sun$chemical<-as.factor(SS_sun$chemical)
Lm1sun=glm(height~chemical,data=SS_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0       0.2103     0.5821   0.361   0.9306  
# standard - crude == 0   1.2299     0.5655   2.175   0.0754 .
# standard - none == 0    1.0196     0.5502   1.853   0.1524 
length(which(SS_sun$chemical =='none')) # 52
length(which(SS_sun$chemical =='standard')) # 59
length(which(SS_sun$chemical =='crude')) # 47

### shade plants
SS_shade = grow_Ss[(grow_Ss$light =="shade"),] 
Lm1shade=glm(height~chemical,data=SS_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0       1.3504     0.3926   3.440  0.00174 ** 
# standard - crude == 0   2.3158     0.3852   6.011  < 0.001 ***
# standard - none == 0    0.9654     0.3528   2.737  0.01700 * 
length(which(SS_shade$chemical =='none')) # 52
length(which(SS_shade$chemical =='standard')) # 57 
length(which(SS_shade$chemical =='crude')) # 38

### no chemical plants
SS_no = grow_Ss[(grow_Ss$chemical =="none"),]
none=glm(height~light,data=SS_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade   -2.0288     0.5193  -3.907 0.000168 ***

### 2HPAA standard plants
SS_stand = grow_Ss[(grow_Ss$chemical =="standard"),]
stand=glm(height~light,data=SS_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -2.0831     0.4395   -4.74 6.22e-06 ***

### whole plant extract plants
SS_crude = grow_Ss[(grow_Ss$chemical =="crude"),]
crude=glm(height~light,data=SS_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -3.1690     0.4968  -6.379 9.61e-09 ***


##########################
### /// X. texanum /// ###
##########################
grow_Xt = grow[(grow$species =="X.texanum"),] 
dotchart(grow_Xt$height, xlab = "growth",
         ylab = "Order of the data") # no outliers
### check distribution
Lm=lm(height~chemical*light,data=grow_Xt)
plot(hist(resid(Lm))) #  normally distributed data.

boxplot(grow_Xt$height~chemical, data = grow_Xt) # Some heterogeneity of variance
boxplot(grow_Xt$height~light, data = grow_Xt) # low heterogeneity of variance

### --- GLS fit with gaussian distribution
grow_Xt$chemical<-as.factor(grow_Xt$chemical)
grow_Xt$light<-as.factor(grow_Xt$light)

# Lm1=gls(height~chemical*light,data=grow_Ss,weights=varIdent(form=~1|as.factor(light)))
Lm1=gls(height~chemical*light,data=grow_Xt,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) 
#                Df Deviance Resid. Df Resid. Dev  Pr(>Chi) 
# chemical           2  66.83074  <.0001
# light              1   0.03707  0.8477
# chemical:light     2   0.06698  0.9353

### --- post-hoc tests for pairwise comparisons
### sun plants
Xt_sun = grow_Xt[(grow_Xt$light =="open"),] 
Lm1sun=glm(height~chemical,data=Xt_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0       0.6621     0.1775   3.730 0.000555 ***
# standard - crude == 0   1.1111     0.1797   6.183  < 1e-04 ***
# standard - none == 0    0.4490     0.1492   3.009 0.007355 ** 
length(which(Xt_sun$chemical =='none')) # 29
length(which(Xt_sun$chemical =='standard')) # 27   
length(which(Xt_sun$chemical =='crude')) # 1

### shade plants
Xt_shade = grow_Xt[(grow_Xt$light =="shade"),] 
Lm1shade=glm(height~chemical,data=Xt_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0       0.6533     0.1113   5.869   <1e-04 ***
# standard - crude == 0   1.1700     0.1245   9.401   <1e-04 ***
# standard - none == 0    0.5167     0.1245   4.152   <1e-04 ***
length(which(Xt_shade$chemical =='none')) # 15
length(which(Xt_shade$chemical =='standard')) # 10
length(which(Xt_shade$chemical =='crude')) # 1

### no chemical plants
XT_no = grow_Xt[(grow_Xt$chemical =="none"),]
none=glm(height~light,data=XT_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade  -0.0087356 0.11346337 -0.076991   0.939

### 2HPAA standard plants
XT_stand = grow_Xt[(grow_Xt$chemical =="standard"),]
stand=glm(height~light,data=XT_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   0.05889    0.27466   0.214     0.83

### whole plant extract plants
XT_crude = grow_Xt[(grow_Xt$chemical =="crude"),]
crude=glm(height~light,data=XT_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# NONE SURVIVED THE EXPERIMENT



####################################################################
####################################################################
### --- Biomass Gained Across Common Garden Treatment Levels --- ###
####################################################################
####################################################################
#biomass<-read.csv("biomass.csv")
biomass2<-read.csv("biomass2.csv") 

biomass2$mass<-as.numeric(biomass2$mass)
biomass2$species<-as.factor(biomass2$species)
biomass2$chemical<-as.factor(biomass2$chemical)
biomass2$light<-as.factor(biomass2$light)

total.mass=filter(biomass2, tissue=="total")

### -------------------------- ###
### ///Final TOTAL Biomass /// ###
### -------------------------- ###
### TOTAL BIOMASS GOALS:
# 1.) Analyze effects of light, chemical and their interaction on total biomass with gaussian GLS.
# 2.) Analyze pairwise biomass differences across shade treatments with post-hoc tests for each species.

##########################
### /// M. maximus /// ###
############s##############
mass_MM = total.mass[(total.mass$species =="M.maximus"),] # retaining only M.maximus
dotchart(mass_MM$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=mass_MM)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals shows normally distributed data.

### check homogeneity of variance
boxplot(mass_MM$mass~chemical, data = mass_MM) # Low heterogeneity of variance
boxplot(mass_MM$mass~light, data = mass_MM) # Very large heterogeneity of variance

### --- GLS fit with gaussian distribution
mass_MM<-na.omit(mass_MM)
Lm1=gls(mass~chemical*light,data=mass_MM,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                     numDF   F-value p-value 
# chemical           2    1.2245  0.2956
# light              1 1023.8719  <.0001
# chemical:light     2   26.1456  <.0001

### --- post-hoc tests for pairwise comparisons
### sun plants
MM_sun = mass_MM[(mass_MM$light =="open"),] 
Lm1sun=glm(mass~chemical,data=MM_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0        78.99      12.74   6.200   <1e-04 ***
# standard - crude == 0    96.87      14.57   6.651   <1e-04 ***
# standard - none == 0     17.88      13.48   1.327    0.379 
length(which(MM_sun$chemical =='none')) # 79
length(which(MM_sun$chemical =='standard')) # 45   
length(which(MM_sun$chemical =='crude')) # 54

### shade plants
MM_shade = mass_MM[(mass_MM$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=MM_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0       2.3652     1.9675   1.202    0.446
# standard - crude == 0  -0.6654     2.1633  -0.308    0.948
# standard - none == 0   -3.0306     1.5202  -1.994    0.111
length(which(MM_shade$chemical =='none')) # 54
length(which(MM_shade$chemical =='standard')) #  26  
length(which(MM_shade$chemical =='crude')) # 13

### no chemical plants
MM_no = mass_MM[(mass_MM$chemical =="none"),]
none=glm(mass~light,data=MM_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade  -192.252      9.511  -20.21   <2e-16 ***

### 2HPAA standard plants
MM_stand = mass_MM[(mass_MM$chemical =="standard"),]
stand=glm(mass~light,data=MM_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -213.16      18.21  -11.71   <2e-16 ***

### whole plant extract plants
MM_crude = mass_MM[(mass_MM$chemical =="crude"),]
crude=glm(mass~light,data=MM_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade  -115.62      15.22  -7.596 1.52e-10 ***


##############################
### /// C. fasciculata /// ###
##############################
mass_CF = total.mass[(total.mass$species =="C.fasciculata"),] 
dotchart(mass_CF$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=mass_CF)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals shows normally distributed data.

### check homogeneity of variance
boxplot(mass_CF$mass~chemical, data = mass_CF) # Some heterogeneity of variance
boxplot(mass_CF$mass~light, data = mass_CF) # Very large heterogeneity of variance

### --- GLS fit with gaussian distribution
mass_CF<-na.omit(mass_CF)
Lm1=gls(mass~chemical*light,data=mass_CF,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                     numDF   F-value p-value  
# chemical           2   2.2328  0.1097
# light              1 354.6666  <.0001
# chemical:light     2  16.5870  <.0001

### --- post-hoc tests for pairwise comparisons
### sun plants
CF_sun = mass_CF[(mass_CF$light =="open"),] 
Lm1sun=glm(mass~chemical,data=CF_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0        6.289      9.974   0.631    0.802    
# standard - crude == 0   45.033      8.987   5.011   <1e-04 ***
# standard - none == 0    38.745      8.401   4.612   <1e-04 ***
length(which(CF_sun$chemical =='none')) # 43
length(which(CF_sun$chemical =='standard')) # 74  
length(which(CF_sun$chemical =='crude')) # 35

### shade plants
CF_shade = mass_CF[(mass_CF$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=CF_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      -0.1909     0.9555  -0.200    0.978
# standard - crude == 0   0.5515     0.9993   0.552    0.844
# standard - none == 0    0.7424     0.7362   1.008    0.568
length(which(CF_shade$chemical =='none')) # 33
length(which(CF_shade$chemical =='standard')) # 24  
length(which(CF_shade$chemical =='crude')) # 11

### no chemical plants
CF_no = mass_CF[(mass_CF$chemical =="none"),]
none=glm(mass~light,data=CF_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade  -50.284      7.409  -6.786 1.15e-11 ***

### 2HPAA standard plants
CF_stand = mass_CF[(mass_CF$chemical =="standard"),]
stand=glm(mass~light,data=CF_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -88.286      9.774  -9.033 1.78e-14 ***

### whole plant extract plants
CF_crude = mass_CF[(mass_CF$chemical =="crude"),]
crude=glm(mass~light,data=CF_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# shade - open == 0   -43.80      10.98  -3.989 0.000248 ***

##########################
### /// S. setaria /// ###
##########################
mass_SS = total.mass[(total.mass$species =="S.scheelei"),] 
dotchart(mass_SS$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=mass_SS)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals shows normally distributed data.

### check homogeneity of variance
boxplot(mass_SS$mass~chemical, data = mass_SS) # Low heterogeneity of variance
boxplot(mass_SS$mass~light, data = mass_SS) # Very large heterogeneity of variance
var.test(mass_SS$mass~light, data = mass_SS) # F = 351.24, num df = 141, denom df = 119, p-value < 2.2e-16

### --- GLS fit with gaussian distribution
mass_SS<-na.omit(mass_SS)
Lm1=gls(mass~chemical*light,data=mass_SS,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) 
#                     numDF   F-value p-value  
# chemical           2    8.9496  0.0002
# light              1  239.9054  <.0001
# chemical:light     2    1.0588  0.3484

### --- post-hoc tests for pairwise comparisons
### sun plants
SS_sun = mass_SS[(mass_SS$light =="open"),] 
Lm1sun=glm(mass~chemical,data=SS_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0        3.567      4.629   0.771    0.721
# standard - crude == 0    7.612      4.589   1.659    0.221
# standard - none == 0     4.045      4.322   0.936    0.617
length(which(SS_sun$chemical =='none')) # 50
length(which(SS_sun$chemical =='standard')) #   52
length(which(SS_sun$chemical =='crude')) # 40

### shade plants
SS_shade = mass_SS[(mass_SS$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=SS_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0       0.8795     0.2610   3.370 0.002146 ** 
# standard - crude == 0   1.0379     0.2564   4.048 0.000147 ***
# standard - none == 0    0.1583     0.2315   0.684 0.772465 
length(which(SS_shade$chemical =='none')) # 43
length(which(SS_shade$chemical =='standard')) # 47   
length(which(SS_shade$chemical =='crude')) # 30

### no chemical plants
SS_no = mass_SS[(mass_SS$chemical =="none"),]
none=glm(mass~light,data=SS_no)
summary(none)
#                   Estimate Std. Error t value Pr(>|z|)  
# lightshade   -27.740      3.771  -7.355 8.09e-11 ***

### 2HPAA standard plants
SS_stand = mass_SS[(mass_SS$chemical =="standard"),]
stand=glm(mass~light,data=SS_stand)
summary(stand)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -31.627      3.255  -9.717 5.49e-16 ***

### whole plant extract plants
SS_crude = mass_SS[(mass_SS$chemical =="crude"),]
crude=glm(mass~light,data=SS_crude)
summary(crude)
#                   Estimate Std. Error t value Pr(>|z|) 
# lightshade   -25.053      3.077  -8.141 1.22e-11 ***


##########################
### /// X. texanum /// ###
##########################
mass_XT = total.mass[(total.mass$species =="X.texanum"),] 
mass_XT$chemical<-as.factor(mass_XT$chemical)
mass_XT$light<-as.factor(mass_XT$light)

dotchart(mass_XT$mass, xlab = "biomass",
         ylab = "Order of the data") # one huge outlier at 16.4 mg, I checked photos of largest plants and it's legitimate
### check distribution
Lm=lm(mass~chemical*light,data=mass_XT)
plot(hist(resid(Lm))) 

### check homogeneity of variance
boxplot(mass_XT$mass~chemical, data = mass_XT) # Low heterogeneity of variance 
boxplot(mass_XT$mass~light, data = mass_XT) # Low heterogeneity of variance

### --- GLS fit with gaussian distribution
mass_XT<-na.omit(mass_XT)
Lm1=gls(mass~chemical*light,data=mass_XT,weights=varIdent(form=~1|as.factor(light)))
summary(Lm1) 
anova(Lm1) 
#                     numDF   F-value p-value  
# chemical           2   5.35977  0.0073
# light              1  15.93748  0.0002
# chemical:light     2   0.44812  0.6410

### --- post-hoc tests for pairwise comparisons
### sun plants
XT_sun$chemical<-as.factor(XT_sun$chemical)
XT_sun = mass_XT[(mass_XT$light =="open"),] 
XT_sun = XT_sun[-c(2,3,5,15),]
Lm1sun=glm(mass~chemical,data=XT_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0       2.9476     0.7619   3.869 0.000293 ***
# standard - crude == 0   2.5333     0.7863   3.222 0.003324 ** 
# standard - none == 0   -0.4143     0.3726  -1.112 0.492139  
length(which(XT_sun$chemical =='none')) # 21
length(which(XT_sun$chemical =='standard')) #  12
length(which(XT_sun$chemical =='crude')) #  1

### shade plants
XT_shade = mass_XT[(mass_XT$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=XT_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      1.56000    0.59020   2.643   0.0208 *
# standard - crude == 0  1.64444    0.61291   2.683   0.0186 *
# standard - none == 0   0.08444    0.33058   0.255   0.9634 
length(which(XT_shade$chemical =='none')) # 15
length(which(XT_shade$chemical =='standard')) # 9
length(which(XT_shade$chemical =='crude')) # 1

### no chemical plants
XT_no = mass_XT[(mass_XT$chemical =="none"),]
none=glm(mass~light,data=XT_no)
summary(none)
#                   Estimate Std. Error z value Pr(>|z|)  
# lightshade   -1.9991     0.8060  -2.480   0.0181 *

### 2HPAA standard plants
XT_stand = mass_XT[(mass_XT$chemical =="standard"),]
stand=glm(mass~light,data=XT_stand)
summary(stand)
#                   Estimate Std. Error z value Pr(>|z|) 
# lightshade   -2.0289     0.8861  -2.290    0.032 * 

### whole plant extract plants
XT_crude = mass_XT[(mass_XT$chemical =="crude"),]
crude=glm(mass~light,data=XT_crude)
summary(crude)
#                   Estimate Std. Error z value Pr(>|z|) 
# # NONE SURVIVED THE EXPERIMENT IN EITHER LIGHT TREATMENT



### ----------------------------- ###
### ///Final Root-Shoot Ratio /// ###
### ----------------------------- ###
### ROOT:SHOOT RATIO ANALYSIS GOALS:
# 1.) Analyze effects of light, chemical and interaction on root:shoot ratio with gaussian GLS.
# 2.) Analyze pairwise root:shoot ratio differences across shade treatments with post-hoc tests for each species.
### ROOT:SHOOT Basics
# R:S is an index that is comparable across species
# R:S > 1 means plant invested more in shoot mass
# R:S < 1 means plant invested more in root mass
# R:S = 1 means plant invested equally in root & shoot mass

biomass2<-read.csv("biomass2.csv") 
rs=filter(biomass2, tissue=="root.shoot")

nrow(rs) # [1] 2352
rs<-na.omit(rs)
nrow(rs) # [1] 818
rs$mass<-as.numeric(rs$mass)
rs$log.mass<-as.numeric(rs$log.mass)
rs$species<-as.factor(rs$species)
rs$chemical<-as.factor(rs$chemical)
rs$light<-as.factor(rs$light)

####################################
### /// ALL Species Together /// ###
####################################
dotchart(rs$mass, xlab = "R:S",
         ylab = "Order of the data")

### check data 
Lm=lm(mass~light,data=rs)
plot(hist(resid(Lm))) # slight right skew
boxplot(rs$mass~light, data = rs) # Very low heterogeneity of variance

### is root:shoot lower in shade overall?
Lm1=gls(mass~light,data=rs)
summary(Lm1) 
anova(Lm1)
#             numDF   F-value p-value
# (Intercept)     1 1546.3963  <.0001
# light           1    9.4142  0.0022
### higher R:S in sun means that sun plants have more shoot mass than root mass
### OPPOSITE = lower R:S in shade means shade plants have relatively more shoot mass

##########################
### /// M. maximus /// ###
##########################
rs_MM = rs[(rs$species =="M.maximus"),] 
dotchart(rs_MM$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=rs_MM)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals follow normal distribution with 2 outliers with high R:S

### check homogeneity of variance
boxplot(rs_MM$mass~chemical, data = rs_MM) # Some heterogeneity of variance across levels
boxplot(rs_MM$mass~light, data = rs_MM) # Very low heterogeneity of variance

### --- GLS fit with gaussian distribution
rs_MM<-na.omit(rs_MM)
Lm1=gls(mass~chemical*light,data=rs_MM,weights=varIdent(form=~1|as.factor(chemical)))
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                 numDF  F-value p-value
# chemical           2  35.1693  <.0001
# light              1   6.5802  0.0109
# chemical:light     2   1.8124  0.1653

### --- post-hoc tests for pairwise comparisons
### sun plants
MM_sun = rs_MM[(rs_MM$light =="open"),] 
Lm1sun=glm(mass~chemical,data=MM_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0     -0.33072    0.06120  -5.404   <1e-04 ***
# standard - crude == 0  0.09850    0.06996   1.408    0.336    
# standard - none == 0   0.42922    0.06473   6.631   <1e-04 ***

### shade plants
MM_shade = rs_MM[(rs_MM$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=MM_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0     -0.11839    0.10140  -1.168 0.467222    
# standard - crude == 0  0.17428    0.11149   1.563 0.257013    
# standard - none == 0   0.29267    0.07835   3.736 0.000557 ***

### no chemical plants
MM_no = rs_MM[(rs_MM$chemical =="none"),]
none=glm(mass~light,data=MM_no)
summary(glht(none, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|)  
# shade - open == 0 -0.05976    0.04126  -1.448    0.148

### 2HPAA standard plants
MM_stand = rs_MM[(rs_MM$chemical =="standard"),]
stand=glm(mass~light,data=MM_stand[-29,])
summary(glht(stand, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0  -0.2253     0.1059  -2.128   0.0334 *
  
### whole plant extract plants
MM_crude = rs_MM[(rs_MM$chemical =="crude"),]
crude=glm(mass~light,data=MM_crude)
summary(glht(crude, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0  -0.2721     0.1257  -2.164   0.0305 *


##############################
### /// C. fasciculata /// ###
##############################
rs_CF = rs[(rs$species =="C.fasciculata"),] 
dotchart(rs_CF$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=rs_CF)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals follow normal distribution with 2 outliers with high R:S

### check homogeneity of variance
boxplot(rs_CF$mass~chemical, data = rs_CF) # Low heterogeneity of variance across levels
boxplot(rs_CF$mass~light, data = rs_CF) # Some heterogeneity of variance

### --- GLS fit with gaussian distribution
rs_CF<-na.omit(rs_CF)
Lm1=gls(mass~chemical*light,data=rs_CF,weights=varIdent(form=~1|as.factor(light)))
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                 numDF  F-value p-value
# chemical           2   5.7295  0.0038
# light              1  77.7649  <.0001
# chemical:light     2   4.8077  0.0091

### --- post-hoc tests for pairwise comparisons
### sun plants
CF_sun = rs_CF[(rs_CF$light =="open"),] 
Lm1sun=glm(mass~chemical,data=CF_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      0.03852    0.07433   0.518  0.86178   
# standard - crude == 0  0.23092    0.06698   3.447  0.00169 **
# standard - none == 0   0.19240    0.06261   3.073  0.00619 **

### shade plants
CF_shade = rs_CF[(rs_CF$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=CF_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0      0.06635    0.05324   1.246    0.422
# standard - crude == 0  0.04373    0.05568   0.785    0.709
# standard - none == 0  -0.02262    0.04103  -0.551    0.844

### no chemical plants
CF_no = rs_CF[(rs_CF$chemical =="none"),]
none=glm(mass~light,data=CF_no)
summary(glht(none, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|)  
# shade - open == 0 -0.18807    0.04974  -3.781 0.000156 ***

### 2HPAA standard plants
CF_stand = rs_CF[(rs_CF$chemical =="standard"),]
stand=glm(mass~light,data=CF_stand)
summary(glht(stand, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0 -0.40309    0.08398    -4.8 1.59e-06 ***

### whole plant extract plants
CF_crude = rs_CF[(rs_CF$chemical =="crude"),]
crude=glm(mass~light,data=CF_crude)
summary(glht(crude, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0 -0.21590    0.06828  -3.162  0.00157 **


##########################
### /// S. setaria /// ###
##########################
rs_SS = rs[(rs$species =="S.scheelei"),] 
dotchart(rs_SS$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=rs_SS)
plot(hist(resid(Lm))) # Plotting a histogram of the model residuals follow normal distribution with 2 outliers with high R:S

### check homogeneity of variance
boxplot(rs_SS$mass~chemical, data = rs_SS) # Low heterogeneity of variance across levels
boxplot(rs_SS$mass~light, data = rs_SS) # Very low heterogeneity of variance

### --- GLS fit with gaussian distribution
rs_SS<-na.omit(rs_SS)
Lm1=gls(mass~chemical*light,data=rs_SS)
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                 numDF  F-value p-value
# chemical           2   5.0832  0.0068
# light              1   0.5920  0.4423
# chemical:light     2   8.9041  0.0002

### --- post-hoc tests for pairwise comparisons
### sun plants
SS_sun = rs_SS[(rs_SS$light =="open"),] 
Lm1sun=glm(mass~chemical,data=SS_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0     -0.03684    0.10767  -0.342    0.937
# standard - crude == 0  0.10390    0.10674   0.973    0.594
# standard - none == 0   0.14074    0.10053   1.400    0.341

### shade plants
SS_shade = rs_SS[(rs_SS$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=SS_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0     -0.52417    0.12310  -4.258   <1e-04 ***
# standard - crude == 0 -0.56635    0.12093  -4.683   <1e-04 ***
# standard - none == 0  -0.04218    0.10920  -0.386    0.921    

### no chemical plants
SS_no = rs_SS[(rs_SS$chemical =="none"),]
none=glm(mass~light,data=SS_no)
summary(glht(none, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|)  
# shade - open == 0 -0.10799    0.08872  -1.217    0.224

### 2HPAA standard plants
SS_stand = rs_SS[(rs_SS$chemical =="standard"),]
stand=glm(mass~light,data=SS_stand)
summary(glht(stand, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0 -0.29091    0.09442  -3.081  0.00206 **

### whole plant extract plants
SS_crude = rs_SS[(rs_SS$chemical =="crude"),]
crude=glm(mass~light,data=SS_crude)
summary(glht(crude, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0   0.3793     0.1583   2.396   0.0166 *

##########################
### /// X. texanum /// ###
##########################
rs_XT = rs[(rs$species =="X.texanum"),] 
dotchart(rs_XT$mass, xlab = "biomass",
         ylab = "Order of the data")
### check distribution
Lm=lm(mass~chemical*light,data=rs_XT)
plot(hist(resid(Lm))) # one outlier with lage ratio is real, I double checked photos

### check homogeneity of variance
boxplot(rs_XT$mass~chemical, data = rs_XT) # Some heterogeneity of variance across levels
boxplot(rs_XT$mass~light, data = rs_XT) # Very low heterogeneity of variance

### --- GLS fit with gaussian distribution
rs_XT<-na.omit(rs_XT)
Lm1=gls(mass~chemical*light,data=rs_XT)
summary(Lm1) 
anova(Lm1) # the interaction is significant so keep both terms. Report p and F
#                 numDF  F-value p-value
# chemical           2   4.55281  0.0145
# light              1   0.50055  0.4820
# chemical:light     2   1.20832  0.3060

### --- post-hoc tests for pairwise comparisons
### sun plants
XT_sun = rs_XT[(rs_XT$light =="open"),] 
Lm1sun=glm(mass~chemical,data=XT_sun)
summary(glht(Lm1sun, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|)    
# none - crude == 0      0.654881   0.273708   2.393   0.0401 *
# standard - crude == 0  0.647318   0.278979   2.320   0.0484 *
# standard - none == 0  -0.007563   0.124094  -0.061   0.9978 

### shade plants
XT_shade = rs_XT[(rs_XT$light =="shade"),] 
Lm1shade=glm(mass~chemical,data=XT_shade)
summary(glht(Lm1shade, mcp(chemical="Tukey")))
#                       Estimate Std. Error z value Pr(>|z|) 
# none - crude == 0       0.5932     0.4555   1.302   0.3823  
# standard - crude == 0   0.9812     0.4730   2.074   0.0899 .
# standard - none == 0    0.3880     0.2551   1.521   0.2707 

### no chemical plants
XT_no = rs_XT[(rs_XT$chemical =="none"),]
none=glm(mass~light,data=XT_no)
summary(glht(none, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|)  
# shade - open == 0 -0.06169    0.17558  -0.351    0.725

### 2HPAA standard plants
XT_stand = rs_XT[(rs_XT$chemical =="standard"),]
stand=glm(mass~light,data=XT_stand)
summary(glht(stand, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# shade - open == 0   0.3338     0.1739    1.92   0.0549 .

### whole plant extract plants
XT_crude = rs_XT[(rs_XT$chemical =="crude"),]
crude=glm(mass~light,data=XT_crude)
summary(glht(crude, mcp(light="Tukey")))
#                   Estimate Std. Error z value Pr(>|z|) 
# NONE SURVIVED THE EXPERIMENT



### ----------------------------------------------------- ###
### /// Germination of seeds with different solutions /// ###
### ----------------------------------------------------- ###

germ<-read.csv('solvent.germination.csv')
germ=germ[,-8]
germ$percent=germ$percent*100
final<-filter(germ,week==4)
final

##########################
### /// M. maximus /// ###
##########################
MMgerm<-filter(final,species=="M.maximus")

MM_solv <- data.frame(
  "germinated" = c(80,0),"not" = c(16,96),
  row.names = c("control","solvent"))
test1<-chisq.test(MM_solv)
test1$p.value # 6.24146e-31

MM_meth <- data.frame(
  "germinated" = c(26,0),"not" = c(16,96),
  row.names = c("methanol","solvent"))
test2<-chisq.test(MM_meth)
test2$p.value # 8.753713e-17

MM_solut<- data.frame(
  "germinated" = c(83,62),"not" = c(17,38),
  row.names = c("control","methanol"))
test2<-chisq.test(MM_solut)
test2$p.value #  0.001538984
# LOWER GERMINATION SUCCESS IN METHANOL

##############################
### /// C. fasciculata /// ###
##############################
CFgerm<-filter(final,species=="C.fasciculata")

CF_solv <- data.frame(
  "germinated" = c(45,1),"not" = c(51,95),
  row.names = c("control","solvent"))
test1<-chisq.test(CF_solv)
test1$p.value # 3.581838e-13

CF_meth <- data.frame(
  "germinated" = c(36,1),"not" = c(6,96),
  row.names = c("methanol","solvent"))
test2<-chisq.test(CF_meth)
test2$p.value #  2.86211e-24

CF_solut<- data.frame(
  "germinated" = c(45,36),"not" = c(51,6),
  row.names = c("control","methanol"))
test2<-chisq.test(CF_solut)
test2$p.value #  4.584111e-05
# HIGHER GERMINATION SUCCESS IN METHANOL

###########################
### /// S. scheelei /// ###
###########################
SSgerm<-filter(final,species=="S.scheelei")

SS_solv <- data.frame(
  "germinated" = c(53,0),"not" = c(43,96),
  row.names = c("control","solvent"))
test1<-chisq.test(SS_solv)
test1$p.value # 4.668209e-17

SS_meth <- data.frame(
  "germinated" = c(33,0),"not" = c(9,96),
  row.names = c("methanol","solvent"))
test2<-chisq.test(SS_meth)
test2$p.value #  2.039914e-22

SS_solut<- data.frame(
  "germinated" = c(53,33),"not" = c(43,9),
  row.names = c("control","methanol"))
test2<-chisq.test(SS_solut)
test2$p.value #  [1] 0.01572909
# HIGHER GERMINATION SUCCESS IN METHANOL

##########################
### /// X. texanum /// ###
##########################
XTgerm<-filter(final,species=="X.texanum")

XT_solv <- data.frame(
  "germinated" = c(22,0),"not" = c(74,96),
  row.names = c("control","solvent"))
test1<-chisq.test(XT_solv)
test1$p.value # 1.954179e-06

XT_meth <- data.frame(
  "germinated" = c(16,0),"not" = c(26,96),
  row.names = c("methanol","solvent"))
test2<-chisq.test(XT_meth)
test2$p.value #  8.105687e-10

XT_solut<- data.frame(
  "germinated" = c(22,16),"not" = c(74,26),
  row.names = c("control","methanol"))
test2<-chisq.test(XT_solut)
test2$p.value #  [1] 0.1031811
# NO DIFFERENCE 


############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################