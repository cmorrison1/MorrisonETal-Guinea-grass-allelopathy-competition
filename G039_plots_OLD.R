#######################################################################
#######################################################################
# ///// Guinea grass Experiment (G039) OLD PLOTS That Got Cut # ///// #
#######################################################################
#######################################################################


getwd()
setwd("~/Desktop/guinea.grass.chemistry/allelopathy/analyses")

library(DescTools)
library(ggplot2)
library(multcomp)
library(vegan)
library(lme4)
library(plotrix)
library(psych)
library(reshape2)
library(psych)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(survival)
library(GGally)
source("coldiss.R")
source("panelutils.R")

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
  

#########################################################################################
#### --- FIGURE 1: 2-Hydroxyphenylacetic Acid Concentration Rancho La Paloma, TX --- ####
#########################################################################################

### ---------------------------------------------------------- ###
### /// Panel B: 2HPAA tissue concentrations with Bargraph /// ###
### ---------------------------------------------------------- ###

### make summary data table with mean, sd, and se for the plot
stats <- data_summary(lp.chem, varname="log.pg.mg", 
                      groupnames=c("type"))
stats

# set factor level order 
stats$type = factor(stats$type,
                    c("soil","leaf", "litter","root"))
# plot it
bar<-ggplot(stats,aes(y=log.pg.mg, x=type)) + 
  geom_bar(stat="identity",
           col="black",
           fill="grey80",
           width=0.7,
           lwd=0.75) +
  geom_errorbar(aes(ymin=log.pg.mg-se, ymax=log.pg.mg+se),
                width=0.3,
                position=position_dodge(0.7),
                col='black') +
  ylab('') +
  xlab('') + 
  scale_x_discrete(labels=c("Soil","Leaves & culms","Leaf litter","Roots")) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        legend.position="right",
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face ="bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 50,hjust=1,size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18),
        plot.margin = margin(t = 10, r = 0, b = 0, l =  0, unit = "pt")
  )
bar















#######################################################################################
#### --- log10 transformed biomass so the graphs are easier to view and compare --- ###
#######################################################################################

biomass2=biomass[,-c(5:10)]
### melt the dataframe so that each root:shoot value is a row
mass.m <- melt(biomass2, id = c("species","chemical","shade","tissue"))
mass.m=mass.m[,-c(5)]
mass.m$log.value=log(mass.m$value +1)
names(mass.m)<-c("species","chemical","light","tissue","mass","log.mass")
write.csv(mass.m, file="biomass2.csv")

### calculate summary stats for new dataframe with log10 values
log.total=filter(mass.m, tissue=='total')
log.total<-log.total[,-4]
mass.stats <- data_summary(log.total, varname="log.mass", 
                           groupnames=c("species","chemical","light"))
mass.stats
write.csv(mass.stats, file="biomass.stats.csv")

### --- Plot the Log10 transformed data
mass.stats<-read.csv("biomass.stats.csv")

### Sun Plants 
log.sun=filter(mass.stats, light == 'open')
# colors for chemical treatments
mass.colors<-c("white","grey80","grey40")
log.sun$chemical <- factor(log.sun$chemical,
                           levels = c("none", "standard", "crude"),
                           labels=c("none","2HPAA","whole plant"))
log.sun$species <- factor(log.sun$species,
                          levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                          labels = c("M. maximus",
                                     "C. fasciculata", 
                                     "S. scheelei",
                                     "X. texanum"))
sun.mass.log<-ggplot(log.sun, aes(x=species, y=log.mass, fill= chemical)) +
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se), width=0.4,position = position_dodge(width=0.75)) +
  scale_fill_manual(values = mass.colors) + 
  scale_x_discrete(labels=c("Megathyrsus maximus", 
                            "Chamaecrista fasciculata",
                            "Setaria scheelei", 
                            "Xanthisma texanum")) +
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("A. Sun Treatment") +
  xlab("") + 
  ylab("") +
  ylim(0,6) +
  #geom_text(x=0.75, y=11, label="a",size=5) +
  #geom_text(x=1.0, y=11.75, label="a",size=5) +
  #geom_text(x=1.25, y=8.75, label="b",size=5) + 
  
  #geom_text(x=1.75, y=10.75, label="a",size=5) +
  #geom_text(x=2.0, y=9.5, label="a",size=5) +
  #geom_text(x=2.25, y=6.75, label="b",size=5) +
  
  #geom_text(x=2.73, y=7, label="a",size=5) +
  #geom_text(x=3.0, y=8, label="b",size=5) +
  #geom_text(x=3.25, y=6, label="c",size=5) +

#geom_text(x=3.75, y=2, label="a",size=5) +
#geom_text(x=4.0, y=2.5, label="b",size=5) +
#geom_text(x=4.25, y=1, label="c",size=5) +
labs(fill = "Chemical Treatment") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=3,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=14),
        axis.text.x = element_text(angle = 50,hjust=1,size=14),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = 20, r = 10, b = 10, l = -5, unit = "pt"))
sun.mass.log

### Shade plants
log.shade=filter(mass.stats, light == 'shade')
# colors for chemical treatments
mass.colors<-c("white","grey80","grey40")
log.shade$chemical <- factor(log.shade$chemical,
                             levels = c("none", "standard", "crude"),
                             labels=c("none","2HPAA","whole plant"))
log.shade$species <- factor(log.shade$species,
                            levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                            labels = c("M. maximus",
                                       "C. fasciculata", 
                                       "S. scheelei",
                                       "X. texanum"))
shade.mass.log<-ggplot(log.shade, aes(x=species, y=log.mass, fill= chemical)) +
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se), width=0.4,position = position_dodge(width=0.75)) +
  scale_fill_manual(values = mass.colors) + 
  scale_x_discrete(labels=c("Megathyrsus maximus", 
                            "Chamaecrista fasciculata",
                            "Setaria scheelei", 
                            "Xanthisma texanum")) +
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("B. Shade Treatment") +
  xlab("") + 
  ylab("") +
  ylim(0,6) +
  #geom_text(x=0.75, y=11, label="a",size=5) +
  #geom_text(x=1.0, y=11.75, label="a",size=5) +
  #geom_text(x=1.25, y=8.75, label="b",size=5) + 
  
  #geom_text(x=1.75, y=10.75, label="a",size=5) +
  #geom_text(x=2.0, y=9.5, label="a",size=5) +
  #geom_text(x=2.25, y=6.75, label="b",size=5) +
  
  #geom_text(x=2.73, y=7, label="a",size=5) +
  #geom_text(x=3.0, y=8, label="b",size=5) +
  #geom_text(x=3.25, y=6, label="c",size=5) +

#geom_text(x=3.75, y=2, label="a",size=5) +
#geom_text(x=4.0, y=2.5, label="b",size=5) +
#geom_text(x=4.25, y=1, label="c",size=5) +
labs(fill = "Chemical Treatment") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=3,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=14),
        axis.text.x = element_text(angle = 50,hjust=1,size=14),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = 20, r = 10, b = 10, l = -5, unit = "pt"))
shade.mass.log

### --- make sun/shade biomass plots into multipanel figure
log.biom=ggarrange(sun.mass.log,shade.mass.log,common.legend = TRUE, legend = "bottom")
# dev plot
dev.new(
  title = "log.biomass",
  width = 12,
  height = 10,
  noRStudioGD = TRUE
)
log.biom
annotate_figure(log.biom,
                left = text_grob(expression(paste('log'[10]," Total Biomass (mg)")), hjust=0.25,vjust=2, rot=90,size = 16),
                bottom = text_grob("Chemical Treatment", vjust=-5, size = 16))



### -------------------------------------------------- ###
### /// Panels C and D: Percent Survival over Time /// ###
### -------------------------------------------------- ###
### GOAL:
# 1.) Visualize survival trends over time with separate time series 
#     plots for plants grown on each treatment for each species.
#   - Scatterplot with lines connecting the data points
surv2<-read.csv("survival.csv")
surv2=surv2[,-c(4)]
head(surv2)

### melt the dataframe so that each row is an individual plant that survived or died (= 1 or 0)
surv.m <- melt(surv2, id = c("species","chemical","shade","week"))
surv.m=surv.m[,-c(5)]
head(surv.m,50)
nrow(surv.m) # 14112

### make summary data table with mean, sd, and se for the plots
surv.stats<-data_summary(surv.m, varname="value", 
                         groupnames=c("species","chemical","shade","week"))
surv.stats$value=surv.stats$value*100
surv.stats$sd=surv.stats$sd*100
surv.stats$se=surv.stats$se*100
head(surv.stats)
# subset summary table for shade and sun light levels
surv.shade<-filter(surv.stats, shade=='shade')
surv.sun<-filter(surv.stats, shade=='open')

### --- Plot line graphs of survival curves for species from the SUN treatment 
time.colors<-c("grey80","grey50","black")

# order factor levels and modify the names
surv.sun$chemical <- factor(surv.sun$chemical,
                            levels = c("none", "standard", "crude"),
                            labels=c("none","2HPAA","whole-plant"))
surv.sun$species <- factor(surv.sun$species,
                           levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                           labels = c("M. maximus",
                                      "C. fasciculata", 
                                      "S. scheelei",
                                      "X. texanum"))

# compact letter display text for pairwise differences in chemical differences
Mm_text_sun <- data.frame(label = c("a", "b", "b"),  
                          species = factor("M. maximus",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(6, 6, 6),
                          y = c(82, 56, 48))
Cf_text_sun <- data.frame(label = c("a", "b", "b"), 
                          species = factor("C. fasciculata",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(6, 6, 6),
                          y = c(82, 49, 36))
Ss_text_sun <- data.frame(label = c("a", "a", "a"), 
                          species = factor("S. scheelei",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(6, 6, 6),
                          y = c(53, 61, 48))
Xt_text_sun <- data.frame(label = c("a", "a", "b"), 
                          species = factor("X. texanum",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(6, 6, 6),
                          y = c(31.25, 26.5, 0))

#
sun.time<-ggplot(surv.sun, aes(x=week, y=value)) +
  geom_line(aes(color=chemical),position=position_dodge(0.3),lwd=1) +
  geom_point(aes(color=chemical),position=position_dodge(0.1),
             size=3, shape=21, fill="white") +
  #  geom_errorbar(aes(color=chemical,ymin=value-se, ymax=value+se), 
  #                width=0.9, position=position_dodge(0.1),lwd=1) + 
  facet_grid(~ species) +
  scale_color_manual(values=time.colors) +
  ylim(0,100) +
  ylab('Recruitment over time (%)') +
  xlab('') +
  scale_x_discrete(breaks=c("0","1","2","3","4","5"),
                   limits=c("0","1","2","3","4","5")) +
  geom_text(data = Mm_text_sun, mapping = aes(x = x, y = y,
                                              label = label),size=8,fontface="bold") +
  geom_text(data = Cf_text_sun, mapping = aes(x = x, y = y,
                                              label = label),size=8,fontface="bold") +
  geom_text(data = Ss_text_sun, mapping = aes(x = x, y = y,
                                              label = label),size=8,fontface="bold") +
  geom_text(data = Xt_text_sun, mapping = aes(x = x, y = y,
                                              label = label),size=8,fontface="bold") +
  ggtitle('') +
  labs(color='Chemical treatment') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color="black"),
        panel.background = element_blank(),
        strip.text = element_text(face = "italic",size=18),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"),panel.spacing = unit(2, "lines"),
        plot.title = element_text(vjust=2,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black",fill='white'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(hjust=1,size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18,face="bold",vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face="bold",vjust=1),
        plot.margin = margin(t = 10, r = 5, b = 0, l = 0, unit = "pt"))
sun.time

### --- Plot line graphs of survival curves for species from the SHADE treatment 
surv.shade$chemical <- factor(surv.shade$chemical,
                              levels = c("none", "standard", "crude"),
                              labels=c("none","2HPAA","whole-plant"))
surv.shade$species <- factor(surv.shade$species,
                             levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                             labels = c("M. maximus",
                                        "C. fasciculata", 
                                        "S. scheelei",
                                        "X. texanum"))

# compact letter display text for pairwise differences in chemical differences
Mm_text_shade <- data.frame(label = c("a", "b", "b"),  
                            species = factor("M. maximus",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(6, 6, 6),
                            y = c(58, 17.5, 28.75))
Cf_text_shade <- data.frame(label = c("a", "a", "a"), 
                            species = factor("C. fasciculata",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(6, 6, 6),
                            y = c(28.25, 35, 13.5))
Ss_text_shade <- data.frame(label = c("a", "a", "b"), 
                            species = factor("S. scheelei",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(6, 6, 6),
                            y = c(53.25, 58.25, 39))
Xt_text_shade <- data.frame(label = c("a", "a", "b"), 
                            species = factor("X. texanum",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(6, 6, 6),
                            y = c(15.75, 10.5, 0))
# 
shade.time<-ggplot(surv.shade, aes(x=week, y=value)) +
  geom_line(aes(color=chemical),position=position_dodge(0.3),lwd=1) +
  geom_point(aes(color=chemical),position=position_dodge(0.1),
             size=3, shape=21, fill="white") +
  geom_errorbar(aes(color=chemical,ymin=value-se, ymax=value+se), 
                width=0.9, position=position_dodge(0.1),lwd=1) + 
  facet_grid(~ species) +
  scale_color_manual(values=time.colors) +
  ylim(0,100) +
  ylab('') +
  xlab('') +
  geom_text(data = Mm_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=8,fontface="bold") +
  geom_text(data = Cf_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=8,fontface="bold") +
  geom_text(data = Ss_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=8,fontface="bold") +
  geom_text(data = Xt_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=8,fontface="bold") +
  scale_x_discrete(breaks=c("0","1","2","3","4","5"),
                   limits=c("0","1","2","3","4","5")) +
  ggtitle('') +
  labs(color='Chemical treatment') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color="black"),
        panel.background = element_blank(),
        strip.text = element_text(face = "italic",size=18),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(vjust=2,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black",fill='white'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(hjust=1,size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18,face="bold",vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face="bold",vjust=1),
        plot.margin = margin(t = 10, r = 5, b = 0, l = 0, unit = "pt"))
shade.time

###################################################################
###################################################################
### --- FIGURE 4: Relative Shoot Growth & BIOMASS Bar Graphs--- ###
###################################################################
###################################################################

### ---------------------------------------- ###
### /// Panels A and B - Relative Growth /// ###
### ---------------------------------------- ###
### GOAL:
# 1.) Visualize variation in final plant height  with dofged bargraphs for each species.

grow<-read.csv("growth2.csv") # this version has no 0s as place holders for X.tex on crude treatment
head(grow,10)
View(grow)

### subset just the final measurements
grow.final = filter(grow,week==4)
### subset just the SHADE treatment
grow.shade = filter(grow.final,shade=='shade')

### --- GG barplot of Final Growth  --- ###
grow.shade$chemical <- factor(grow.shade$chemical,
                              levels = c("none", "standard", "crude"),
                              labels = c("none","2HPAA","whole-plant"))
# order the species in order of descending survival
grow.shade$species <- factor(grow.shade$species,
                             levels = c("M.maximus", "C.fasciculata","S.scheelei","X.texanum"),
                             labels = c("M. maximus",
                                        "C. fasciculata", 
                                        "S. scheelei",
                                        "X. texanum"))
# colors for chemical treatments
grow.colors<-c("white","grey80","grey40")
# plot it
shade.grow<-ggplot(grow.shade, aes(x=species, y=average,fill=chemical)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=average-SE, ymax=average+SE), width=0.4,position = position_dodge(width=0.75)) +
  scale_fill_manual(values = grow.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("Shade treatment") +
  xlab("") + 
  ylab("") +
  ylim(0,33) +
  geom_text(x=0.75, y=12, label="a",size=5) +
  geom_text(x=1.0, y=12.75, label="a",size=5) +
  geom_text(x=1.25, y=9.75, label="b",size=5) + 
  geom_text(x=1.75, y=11.75, label="a",size=5) +
  geom_text(x=2.0, y=10.5, label="a",size=5) +
  geom_text(x=2.25, y=8.25, label="b",size=5) +
  geom_text(x=2.73, y=8.5, label="a",size=5) +
  geom_text(x=3.0, y=9, label="b",size=5) +
  geom_text(x=3.25, y=7, label="c",size=5) +
  geom_text(x=3.75, y=3.5, label="a",size=5) +
  geom_text(x=4.0, y=4, label="b",size=5) +
  geom_text(x=4.25, y=2, label="c",size=5) +
  labs(fill = "Chemical treatment") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=16),
        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
shade.grow


### subset just the SUN treatment
grow.sun = filter(grow.final,shade=='open')
# order the chemical levels in  order of descending survival
grow.sun$chemical <- factor(grow.sun$chemical,
                            levels = c("none", "standard", "crude"),
                            labels = c("none","2HPAA","whole-plant"))
# order the species in order of descending survival
grow.sun$species <- factor(grow.sun$species,
                           levels = c("M.maximus", "C.fasciculata","S.scheelei","X.texanum"),
                           labels = c("M. maximus",
                                      "C. fasciculata", 
                                      "S. scheelei",
                                      "X. texanum"))
# colors for chemical treatments
grow.colors<-c("white","grey80","grey40")
# plot it
sun.grow<-ggplot(grow.sun, aes(x=species, y=average,fill=chemical)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=average-SE, ymax=average+SE), width=0.4,position = position_dodge(width=0.75)) +
  scale_fill_manual(values = grow.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("Sun treatment") +
  xlab("") + 
  ylab("Shoot height (cm)") +
  ylim(0,33) +
  geom_text(x=0.75, y=34, label="a",size=5) +
  geom_text(x=1.0, y=32, label="a",size=5) +
  geom_text(x=1.25, y=27.75, label="b",size=5) + 
  geom_text(x=1.75, y=11.75, label="a",size=5) +
  geom_text(x=2.0, y=15.75, label="b",size=5) +
  geom_text(x=2.25, y=12, label="a",size=5) +
  geom_text(x=2.73, y=10, label="a",size=5) +
  geom_text(x=3.0, y=11, label="a",size=5) +
  geom_text(x=3.25, y=9.75, label="a",size=5) +
  geom_text(x=3.75, y=3, label="a",size=5) +
  geom_text(x=4.0, y=3.5, label="b",size=5) +
  geom_text(x=4.25, y=2, label="c",size=5) +
  labs(fill = "Chemical treatment") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=16),
        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
sun.grow


### -------------------------------------- ###
### /// Panels C and D - Total Biomass /// ###
### -------------------------------------- ###
biomass<-read.csv("biomass.csv")
head(biomass,10)
### GOAL:
# 1.) Visualize variation in final above:blow ground mass with dodged bargraphs for each species.

### subset out just below and above ground mass
total = filter(biomass,tissue =='total')
### subset just SHADE treatment
shade=filter(total, shade == 'shade')
# colors for chemical treatments
mass.colors<-c("white","grey80","grey40")

shade$chemical <- factor(shade$chemical,
                         levels = c("none", "standard", "crude"),
                         labels=c("none","2HPAA","whole-plant"))
shade$species <- factor(shade$species,
                        levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                        labels = c("M. maximus",
                                   "C. fasciculata", 
                                   "S. scheelei",
                                   "X. texanum"))

shade.mass<-ggplot(shade, aes(x=species, y=AVG, fill= chemical)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=AVG-SE, ymax=AVG+SE), width=0.4,position = position_dodge(width=0.75)) +
  scale_fill_manual(values = mass.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  ylim(0,250) +
  geom_text(x=0.75, y=22, label="a",size=5) +
  geom_text(x=1.0, y=19, label="a",size=5) +
  geom_text(x=1.25, y=20, label="a",size=5) +
  geom_text(x=1.75, y=21, label="a",size=5) +
  geom_text(x=2.0, y=23, label="a",size=5) +
  geom_text(x=2.25, y=21, label="a",size=5) +
  geom_text(x=2.73, y=16, label="a",size=5) +
  geom_text(x=3.0, y=17, label="a",size=5) +
  geom_text(x=3.25, y=15, label="b",size=5) +
  geom_text(x=3.75, y=14, label="a",size=5) +
  geom_text(x=4.0, y=14, label="a",size=5) +
  geom_text(x=4.25, y=12, label="b",size=5) +
  labs(fill = "Chemical treatment") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=16),
        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))
shade.mass


### subset just sun treatment
sun=filter(total, shade == 'open')
# colors for chemical treatments
mass.colors<-c("white","grey80","grey40")
sun$chemical <- factor(sun$chemical,
                       levels = c("none", "standard", "crude"),
                       labels=c("none","2HPAA","whole-plant"))
sun$species <- factor(sun$species,
                      levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                      labels = c("M. maximus",
                                 "C. fasciculata", 
                                 "S. scheelei",
                                 "X. texanum"))
#
sun.mass<-ggplot(sun, aes(x=species, y=AVG, fill= chemical)) +
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=AVG-SE, ymax=AVG+SE), width=0.4,position = position_dodge(width=0.75)) +
  scale_fill_manual(values = mass.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("") +
  xlab("") + 
  ylab("Total biomass (mg)") +
  ylim(0,250) +
  geom_text(x=0.75, y=221, label="a",size=5) +
  geom_text(x=1.0, y=244, label="a",size=5) +
  geom_text(x=1.25, y=142, label="b",size=5) + 
  geom_text(x=1.75, y=78, label="a",size=5) +
  geom_text(x=2.0, y=116, label="b",size=5) +
  geom_text(x=2.25, y=71, label="a",size=5) +
  geom_text(x=2.73, y=48, label="a",size=5) +
  geom_text(x=3.0, y=52, label="a",size=5) +
  geom_text(x=3.25, y=44, label="a",size=5) +
  geom_text(x=3.75, y=17, label="a",size=5) +
  geom_text(x=4.0, y=18, label="a",size=5) +
  geom_text(x=4.25, y=15, label="b",size=5) +
  labs(fill = "Chemical treatment") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=16),
        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))
sun.mass


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################