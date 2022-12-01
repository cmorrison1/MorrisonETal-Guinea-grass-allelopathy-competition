##################################################################################
##################################################################################
##################################################################################
# ------------------------------------------------------------------------------ #
# /////// Guinea grass Allelopathy-Competition Model Comparisons (G039)# /////// #
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


library(ggplot2)
library(nlme) ### Mixed effects modelpackage
library(MASS) #for glmmPQL
library(dplyr)
library(tidyr)

getwd()
setwd("~/Desktop/guinea.grass.chemistry/allelopathy/analyses")



### ----------------------------------------------- ###
### /// Final Percent Survival Model Comparison /// ###
### ----------------------------------------------- ###
surv<-read.csv("survival2.csv")
final = filter(surv,week==5)


##########################
### /// M. maximus /// ###
##########################
surv_MM = final[(final$species =="M.maximus"),] #retaining only M.maximus

### - GLM fit with logit link style negative binomial distributions
mfull=glm.nb(result~chemical*light, link=log,data=surv_MM)
summary(mfull) 
# Coefficients:
#                             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                  -0.5776     0.1348  -4.284 1.84e-05 ***
# chemicalnone                  0.3747     0.1752   2.139   0.0324 *  
# chemicalstandard             -0.1787     0.1998  -0.894   0.3711    
# lightshade                   -1.1741     0.2775  -4.231 2.33e-05 ***
# chemicalnone:lightshade       0.8351     0.3272   2.552   0.0107 *  
# chemicalstandard:lightshade   0.6777     0.3667   1.848   0.0646 .  

# (Dispersion parameter for Negative Binomial(12712.68) family taken to be 1)
# Null deviance: 413.89  on 587  degrees of freedom
# Residual deviance: 357.08  on 582  degrees of freedom
# AIC: 937.1


### - model only testing effect of allelochemicals
mchem <- glm.nb(result~chemical, link=log,data=surv_MM)
summary(mchem)
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -1.0014     0.1178  -8.497  < 2e-16 ***
# chemicalnone       0.6433     0.1456   4.419 9.89e-06 ***
#  chemicalstandard   0.0274     0.1655   0.166    0.869    

anova(mchem, mfull)
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# 1         chemical 15949.65       585       -952.5060                                   
# 2 chemical * light 12712.68       582       -923.0979 1 vs 2     3 29.40813 1.838057e-06
#### light is statistically significant predictor of seedling recruitment

### - model only testing effect of light
mlight <- glm.nb(result~light, link=log,data=surv_MM)
summary(mlight)
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.48508    0.07433  -6.526 6.76e-11 ***
# lightshade  -0.57352    0.12381  -4.632 3.62e-06 ***   

anova(mlight, mfull)
# Model    theta Resid. df    2 x log-lik.   Test    df  LR stat.      Pr(Chi)
# 1            light 16008.78       586       -957.5649                                   
# 2 chemical * light 12712.68       582       -923.0979 1 vs 2     4 34.46709 5.976383e-07

### --- Checking model assumption --- ###
mpois <- glm(result~chemical*light, data=surv_MM, family = "poisson")
pchisq(2 * (logLik(mfull) - logLik(mpois)), df = 1, lower.tail = FALSE)
# 'log Lik.' 1.715863e-119 (df=6)

### --- get the confidence intervals for the coefficients by profiling the likelihood function --- ###
(est <- cbind(Estimate = coef(mfull), confint(mfull)))
#                               Estimate       2.5 %     97.5 %
# (Intercept)                 -0.5776343 -0.85409781 -0.3244931
# chemicalnone                 0.3746934  0.03416961  0.7225251
# chemicalstandard            -0.1786918 -0.57385607  0.2118931
# lightshade                  -1.1741198 -1.74789565 -0.6527657
# chemicalnone:lightshade      0.8351445  0.20899952  1.4973973
# chemicalstandard:lightshade  0.6776830 -0.03285341  1.4106837

### --- looking at incident rate ratios rather than coefficients --- ###
# To do this, we can exponentiate our model coefficients.
exp(est)
#                              Estimate     2.5 %    97.5 %
# (Intercept)                 0.5612245 0.4256671 0.7228937
# chemicalnone                1.4545455 1.0347601 2.0596275
# chemicalstandard            0.8363636 0.5633489 1.2360158
# lightshade                  0.3090909 0.1741400 0.5206040
# chemicalnone:lightshade     2.3051471 1.2324444 4.4700396
# chemicalstandard:lightshade 1.9693095 0.9676804 4.0987568


#---------------------------------------------------------------------------------------------------- #
###############################
### /// C. fasciculata /// ###
##############################
surv_CF = final[(final$species =="C.fasciculata"),] 

### - GLM fit with logit link style negative binomial distributions
mfull=glm.nb(result~chemical*light,link='log',data=surv_CF)
summary(mfull) 
#                             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 -1.02962    0.16903  -6.091 1.12e-09 ***
# chemicalnone                 0.29480    0.22327   1.320  0.18671    
# chemicalstandard             0.82668    0.20267   4.079 4.52e-05 ***
# lightshade                  -0.99040    0.32480  -3.049  0.00229 ** 
# chemicalnone:lightshade      0.66661    0.39521   1.687  0.09165 .  
# chemicalstandard:lightshade -0.09579    0.39375  -0.243  0.80779    

### - model only testing effect of allelochemicals
mchem <- glm.nb(result~chemical, link=log,data=surv_CF)
summary(mchem)
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -1.4069     0.1443  -9.747  < 2e-16 ***
# chemicalnone       0.5232     0.1822   2.873  0.00407 ** 
# chemicalstandard   0.8016     0.1737   4.614 3.94e-06 ***   

anova(mchem, mfull)
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# 1         chemical 11595.958       585       -879.7576                                   
# 2 chemical * light  9554.061       582       -839.7485 1 vs 2     3 40.00917 1.060751e-08

### - model only testing effect of light
mlight <- glm.nb(result~light, link=log,data=surv_CF)
summary(mlight)
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.59598    0.07857  -7.585 3.31e-14 ***
# lightshade  -0.78353    0.14031  -5.584 2.35e-08 ***

anova(mlight, mfull)
# Response: result
#              Model     theta Resid. df    2 x log-lik.   Test    df  LR stat.      Pr(Chi)
# 1            light 11099.196       586       -869.2765                                   
# 2 chemical * light  9554.061       582       -839.7485 1 vs 2     4 29.52807 6.105607e-06


# ---------------------------------------------------------------------------------------------------- #
###########################
### /// S. scheelei /// ###
###########################
surv_SS = final[(final$species =="S.scheelei"),] 

### - GLM fit with logit link style negative binomial distributions
mfull=glm.nb(result~chemical*light,link=log,data=surv_SS)
summary(mfull) 
#                            Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                  -0.7348     0.1459  -5.038 4.71e-07 ***
# chemicalnone                  0.1011     0.2013   0.502    0.615    
# chemicalstandard              0.2274     0.1955   1.163    0.245    
# lightshade                   -0.2126     0.2182  -0.974    0.330    
# chemicalnone:lightshade       0.2316     0.2927   0.791    0.429    
# chemicalstandard:lightshade   0.1781     0.2865   0.622    0.534  
anova(mfull, test = "F")
#                Df Deviance Resid. Df Resid. Dev      F  Pr(>F)  
# NULL                             587     399.72                 
# chemical        2   4.9318       585     394.78 2.4659 0.08493 .
# light           1   0.3268       584     394.46 0.3268 0.56752  
# chemical:light  2   0.6719       582     393.78 0.3359 0.71467 


# ---------------------------------------------------------------------------------------------------- #
##########################
### /// X. texanum /// ###
##########################
surv2<-read.csv("survival3.csv") # 0 values for whole-plant extract added to database for model comparison so coefficients can be calculated 
final2 = filter(surv2,week==5)
#surv<-read.csv("survival2.csv")
#final = filter(surv,week==5)

surv_XT = final2[(final2$species =="X.texanum"),] 

### - GLM fit with logit link style negative binomial distributions
mfull=glm(result~chemical*light,family = binomial,data=surv_XT)
summary(mfull)
#                               Estimate Std. Error z value Pr(>|z|)
# (Intercept)                 -3.45526    0.58640  -5.892 3.81e-09 ***
# chemicalnone                 2.58845    0.62677   4.130 3.63e-05 ***
# chemicalstandard             2.48842    0.62848   3.959 7.51e-05 ***
# lightshade                  -1.11945    1.16368  -0.962    0.336    
# chemicalnone:lightshade      0.27547    1.21731   0.226    0.821    
# chemicalstandard:lightshade -0.08846    1.23152  -0.072    0.943  
anova(mfull, test = "F")
#                Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                             587     485.87                      
# chemical        2   48.160       585     437.71 24.0800 3.485e-11 ***
# light           1   16.356       584     421.35 16.3563 5.248e-05 ***
# chemical:light  2    0.468       582     420.88  0.2341    0.7913      



### ----------------------------------------- ###
### /// Relative Height Model Comparisons /// ###
### ----------------------------------------- ###
grow<-read.csv("growth3.csv")


##########################
### /// M. maximus /// ###
##########################
grow_MM = grow[(grow$species =="M.maximus"),] #retaining only M.maximus

### - GLM fit with logit link style negative binomial distributions
mfull=gls(height~chemical*light,data=grow_MM,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  24.454545 0.8493805  28.791036  0.0000
# chemicalnone                  7.250455 1.1033777   6.571145  0.0000
# chemicalstandard              4.865020 1.2585894   3.865454  0.0001
# lightshade                  -17.342781 1.2983813 -13.357233  0.0000
# chemicalnone:lightshade      -4.951693 1.5714328  -3.151069  0.0018
# chemicalstandard:lightshade  -4.826784 1.7702797  -2.726566  0.0068

### - model only testing effect of allelochemicals
mchem <- gls(height~chemical,weights=varIdent(form=~1|as.factor(light)),data=grow_MM, method = "ML")
summary(mchem)
# Coefficients:
#                 Value Std.Error   t-value p-value
# (Intercept)      23.966261 0.8466353 28.307655   0.000
# chemicalnone      6.343085 1.0918420  5.809527   0.000
# chemicalstandard  4.156751 1.2454143  3.337645   0.001

anova(mchem,mfull)
#     Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# mchem     1  5 2099.500 2117.727 -1044.7501                        
# mfull     2  8 1764.571 1793.735  -874.2855 1 vs 2 340.9292  <.0001

### - model only testing effect of light
mlight <- gls(height~light, weights=varIdent(form=~1|as.factor(light)),data=grow_MM, method = "ML")
summary(mlight)
# Coefficients:
#                 Value Std.Error   t-value p-value
# (Intercept)  28.89558 0.5187635  55.70088       0
# lightshade  -20.48872 0.6634396 -30.88257       0

anova(mlight, mfull)
#     Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# mlight     1  4 1804.084 1818.666 -898.0423                       
# mfull      2  8 1764.571 1793.735 -874.2855 1 vs 2 47.5136  <.0001


# ---------------------------------------------------------------------------------------------------- #
##############################
### /// C. fasciculata /// ###
##############################
grow_Cf = grow[(grow$species =="C.fasciculata"),] 

### - GLM fit with logit link style negative binomial distributions
mfull=gls(height~chemical*light,data=grow_Cf,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  9.431429 0.6313767 14.937878  0.0000
# chemicalnone                -0.303769 0.8339626 -0.364248  0.7160
# chemicalstandard             4.001071 0.7569941  5.285472  0.0000
# lightshade                  -4.269890 0.9060190 -4.712804  0.0000
# chemicalnone:lightshade      3.989289 1.1310039  3.527211  0.0005
# chemicalstandard:lightshade -0.873721 1.0947930 -0.798070  0.4257

### - model only testing effect of allelochemicals
mchem <- gls(height~chemical, weights=varIdent(form=~1|as.factor(light)),data=grow_Cf, method = "ML")
summary(mchem)
#                    Value Std.Error   t-value p-value
#(Intercept)      8.003024 0.5433226 14.729785  0.0000
# chemicalnone     0.985819 0.6797904  1.450181  0.1484
# chemicalstandard 3.816681 0.6545381  5.831106  0.0000

anova(mchem, mfull)
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# mchem     1  5 1301.556 1318.876 -645.7782                       
# mfull     2  8 1232.639 1260.350 -608.3197 1 vs 2  74.917  <.0001

### - model only testing effect of light
mlight <- gls(height~light, weights=varIdent(form=~1|as.factor(light)),data=grow_Cf, method = "ML")
summary(mlight)
#                Value Std.Error  t-value p-value
# (Intercept) 11.31914 0.3344419 33.84485       0
# lightshade  -3.32319 0.4571191 -7.26986       0   

anova(mlight, mfull)
# Response: result
#      Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# mlight     1  4 1290.983 1304.839 -641.4916                        
# mfull      2  8 1232.639 1260.350 -608.3197 1 vs 2 66.34382  <.0001


# ---------------------------------------------------------------------------------------------------- #
###########################
### /// S. scheelei /// ###
###########################
grow_Ss = grow[(grow$species =="S.scheelei"),] 

### - GLM fit with logit link style negative binomial distributions
mfull=gls(height~chemical*light,data=grow_Ss,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  7.353191 0.4220463 17.422714  0.0000
# chemicalnone                 0.210270 0.5823392  0.361078  0.7183
# chemicalstandard             1.229859 0.5657011  2.174045  0.0305
# lightshade                  -3.168981 0.5168173 -6.131724  0.0000
# chemicalnone:lightshade      1.140135 0.7022233  1.623607  0.1055
# chemicalstandard:lightshade  1.085930 0.6843338  1.586843  0.1136
anova(mfull, test="F")
#                numDF   F-value p-value
# (Intercept)        1 2445.4030  <.0001
# chemical           2   17.7050  <.0001
# light              1   74.0608  <.0001
# chemical:light     2    1.6570  0.1925


# ---------------------------------------------------------------------------------------------------- #
##########################
### /// X. texanum /// ###
##########################
grow_Xt = grow[(grow$species =="X.texanum"),] 

### - GLM fit with logit link style negative binomial distributions
mfull=gls(height~chemical*light,data=grow_Xt,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  0.0000000 0.6328013  0.0000000  1.0000
# chemicalnone                 0.6620690 0.6436192  1.0286656  0.3069
# chemicalstandard             1.1111111 0.6444133  1.7242212  0.0887
# lightshade                   0.0000000 0.7368757  0.0000000  1.0000
# chemicalnone:lightshade     -0.0087356 0.7525272 -0.0116084  0.9908
# chemicalstandard:lightshade  0.0588889 0.7563542  0.0778589  0.9381
anova(mfull, test="F")
#                numDF   F-value p-value
# (Intercept)        1 230.66903  <.0001
# chemical           2  12.57087  <.0001
# light              1   0.03441  0.8533
# chemical:light     2   0.04405  0.9569


# ---------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------- #


### ----------------------------------------- ###
### /// Biomass Model Comparisons /// ###
### ----------------------------------------- ###
biomass2<-read.csv("biomass2.csv") 
total.mass=filter(biomass2, tissue=="total")


##########################
### /// M. maximus /// ###
##########################
mass_MM = total.mass[(total.mass$species =="M.maximus"),] 

### - GLM fit with logit link style negative binomial distributions
mass_MM<-na.omit(mass_MM)
mfull=gls(mass~chemical*light,data=mass_MM,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  119.46296  9.846361  12.132702       0
# chemicalnone                  78.99273 12.775798   6.182998       0
# chemicalstandard              96.87037 14.604513   6.632906       0
# lightshade                  -115.62450 10.001923 -11.560227       0
# chemicalnone:lightshade      -76.62749 12.924858  -5.928691       0
# chemicalstandard:lightshade  -97.53575 14.762224  -6.607118       0

### - model only testing effect of allelochemicals
mchem <- gls(mass~chemical, weights=varIdent(form=~1|as.factor(light)),data=mass_MM, method = "ML")
summary(mchem)
# Coefficients:
#                      Value Std.Error   t-value p-value
# chemicalnone      62.06818 12.628112  4.915080       0
# chemicalstandard  80.95427 14.426541  5.611482       0 

anova(mchem, mfull)
#     Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# mchem     1  5 3265.843 3283.854 -1627.921                        
# mfull     2  8 2646.634 2675.451 -1315.317 1 vs 2 625.2093  <.0001

### - model only testing effect of light
mlight <- gls(mass~light, weights=varIdent(form=~1|as.factor(light)),data=mass_MM, method = "ML")
summary(mlight)
#                 Value Std.Error   t-value p-value
# (Intercept)  179.0112  6.164898  29.03718       0
# lightshade  -173.9854  6.201001 -28.05764       0

anova(mlight, mfull)
#    Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# mlight     1  4 2691.466 2705.874 -1341.733                        
# mfull      2  8 2646.634 2675.451 -1315.317 1 vs 2 52.83205  <.0001


# ---------------------------------------------------------------------------------------------------- #
##############################
### /// C. fasciculata /// ###
##############################
mass_CF = total.mass[(total.mass$species =="C.fasciculata"),] 

### - GLM fit with logit link style negative binomial distributions
mass_CF<-na.omit(mass_CF)
mfull=gls(mass~chemical*light,data=mass_CF,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  51.38571  7.433850  6.912396  0.0000
# chemicalnone                  6.28870 10.012139  0.628108  0.5306
# chemicalstandard             45.03320  9.022174  4.991392  0.0000
# lightshade                  -43.80390  7.478968 -5.856944  0.0000
# chemicalnone:lightshade      -6.47961 10.056841 -0.644299  0.5201
# chemicalstandard:lightshade -44.48169  9.076390 -4.900813  0.0000

### - model only testing effect of allelochemicals
mchem <- gls(mass~chemical, weights=varIdent(form=~1|as.factor(light)),data=mass_CF, method = "ML")
summary(mchem)
#                    Value Std.Error   t-value p-value
# (Intercept)       7.729528 0.8151916  9.481854  0.0000
# chemicalnone     -0.269042 0.9415369 -0.285747  0.7753
# chemicalstandard  0.692326 0.9844535  0.703259  0.4826

anova(mchem, mfull)
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# mchem     1  5 2110.324 2127.293 -1050.1622                        
# mfull     2  8 1920.612 1947.761  -952.3059 1 vs 2 195.7127  <.0001

### - model only testing effect of light
mlight <- gls(mass~light, weights=varIdent(form=~1|as.factor(light)),data=mass_CF, method = "ML")
summary(mlight)
#              Estimate Std. Error z value Pr(>|z|)    
#                 Value Std.Error   t-value p-value
# (Intercept)  75.08882   3.92307  19.14032       0
# lightshade  -67.40499   3.93688 -17.12142       0

anova(mlight, mfull)
# Response: result
#    Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# mlight     1  4 1945.410 1958.985 -968.7051                        
# mfull      2  8 1920.612 1947.761 -952.3059 1 vs 2 32.79856  <.0001


# ---------------------------------------------------------------------------------------------------- #
###########################
### /// S. scheelei /// ###
###########################
mass_SS = total.mass[(total.mass$species =="S.scheelei"),] 

### - GLM fit with logit link style negative binomial distributions
mass_SS<-na.omit(mass_SS)
mfull=gls(mass~chemical*light,data=mass_SS,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  27.512500  3.453532  7.966483  0.0000
# chemicalnone                  3.567500  4.633399  0.769953  0.4420
# chemicalstandard              7.612500  4.593626  1.657187  0.0987
# lightshade                  -25.052500  3.459324 -7.242022  0.0000
# chemicalnone:lightshade      -2.687965  4.640728 -0.579212  0.5630
# chemicalstandard:lightshade  -6.574628  4.600761 -1.429030  0.1542
anova(mfull, test="F")
#                numDF   F-value p-value
# (Intercept)        1 1068.7434  <.0001
# chemical           2    8.9686  0.0002
# light              1  239.4734  <.0001
# chemical:light     2    1.0569  0.3490


# ---------------------------------------------------------------------------------------------------- #
##########################
### /// X. texanum /// ###
##########################
mass_XT = total.mass[(total.mass$species =="X.texanum"),] 

### - GLM fit with logit link style negative binomial distributions
mass_XT<-na.omit(mass_XT)
mfull=gls(mass~chemical*light,data=mass_XT,weights=varIdent(form=~1|as.factor(light)), method = "ML")
summary(mfull)
#                                  Value Std.Error    t-value p-value
# (Intercept)                  0.000000  2.009169  0.0000000  1.0000
# chemicalnone                 3.559091  2.098508  1.6960100  0.0952
# chemicalstandard             3.673333  2.138923  1.7173749  0.0912
# lightshade                   0.000000  2.082378  0.0000000  1.0000
# chemicalnone:lightshade     -1.999091  2.177892 -0.9179017  0.3624
# chemicalstandard:lightshade -2.028889  2.222858 -0.9127390  0.3651
anova(mfull, test="F")
#                numDF   F-value p-value
# (Intercept)        1 133.65236  <.0001
# chemical           2   5.43329  0.0068
# light              1  15.74328  0.0002
# chemical:light     2   0.44201  0.6449


############################################################################################################
############################################################################################################
############################################################################################################