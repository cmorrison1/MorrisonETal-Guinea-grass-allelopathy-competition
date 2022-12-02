##############################################################################
##############################################################################
##############################################################################
# -------------------------------------------------------------------------- #
# ///// Guinea grass Allelopathy-Competition Experiment PLOTS (G039) # ///// #
# -------------------------------------------------------------------------- #
##############################################################################
##############################################################################
##############################################################################

# Colin Richard Morrison 
# PhD Candidate
# The University of Texas at Austin 
# Department of Integrative Biology 
# Graduate Program In Ecology, Evolution and Behavior
# crmorrison@utexas.edu

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
}

#############################################################################################
#############################################################################################
#### --- FIGURE 1: 2-Hydroxyphenylacetic Acid Concentration in the field - south, TX --- ####
#############################################################################################
#############################################################################################
### GOALS:
# 1.) Visualize 2-HPAA concentrations across space and within tissue types 


### ---------------------------------------------------------- ###
### /// Panel A: 2HPAA concentrations across sites sampled /// ###
### ---------------------------------------------------------- ###
dats<-read.csv("LP.2HPAA.data.csv")
log.pg.mg<-log(dats$pg.mg+1) # add a constant to deal with negative values
dats=cbind.data.frame(dats,log.pg.mg)
dats

# subset Df for each substrate
names(dats)
leaf<-dats[which(dats$type=='leaf'),names(dats) %in%
             c("site","longitude","latitude","type","pg.mg","log.pg.mg")]
litter<-dats[which(dats$type=='litter'),names(dats) %in%
               c("site","longitude","latitude","type","pg.mg","log.pg.mg")]
root<-dats[which(dats$type=='root'),names(dats) %in%
             c("site","longitude","latitude","type","pg.mg","log.pg.mg")]
soil<-dats[which(dats$type=='soil'),names(dats) %in%
             c("site","longitude","latitude","type","pg.mg","log.pg.mg")]


A=ggplot(soil, aes(x=longitude, y=latitude)) +
  geom_point(aes(fill= abs(log.pg.mg),
                 size=abs(log.pg.mg)),
             shape = 21,color="black") +
  ylim(27.09,27.29) +
  xlim(-98,-97.93) +
  ggtitle("Soil") +
  scale_alpha_continuous(limits=c(0,5), 
                         labels = c("0-1", "1", "2", "3","4","5")) +
  scale_size_continuous(limits=c(0,5),
                        labels = c("0-1", "1", "2", "3","4","5")) +
  scale_fill_distiller(limits=c(0,5),
                       labels = c("0-1", "1", "2", "3","4","5"),
                       type = "seq",
                       direction = 1,
                       palette = "Greys") +
  xlab("") +
  ylab("") +
  guides(fill= guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)")))), 
         size=guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)"))))) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size=18,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 20, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 50,size=18,hjust=1),
        plot.margin = margin(t = 10, r = 20, b = 0, l = 20, unit = "pt")
  )
A

B=ggplot(leaf, aes(x=longitude, y=latitude)) +
  geom_point(aes(fill= abs(log.pg.mg),
                 size=abs(log.pg.mg)),
             shape = 21,color="black") +
  ylim(27.09,27.29) +
  xlim(-98,-97.93) +
  xlab("") +
  ylab("") +
  ggtitle("Leaves & culms") +
  scale_alpha_continuous(limits=c(0,5), 
                         labels = c("0-1", "1", "2", "3","4","5")) +
  scale_size_continuous(limits=c(0,5),
                        labels = c("0-1", "1", "2", "3","4","5")) +
  scale_fill_distiller(limits=c(0,5),
                       labels = c("0-1", "1", "2", "3","4","5"),
                       type = "seq",
                       direction = 1,
                       palette = "Greys") +
  guides(fill= guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)")))), 
         size=guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)"))))) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size=18,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 20, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 50,size=18,hjust=1),
        plot.margin = margin(t = 10, r = 20, b = 0, l = 20, unit = "pt")
  )
B

C=ggplot(litter, aes(x=longitude, y=latitude)) +
  geom_point(aes(fill= log.pg.mg,
                 size=log.pg.mg),
             shape = 21,color="black") +
  ylim(27.09,27.29) +
  xlim(-98,-97.93) +
  xlab("") +
  ylab("") +
  ggtitle("Leaf litter") +
  scale_alpha_continuous(limits=c(0,5), 
                         labels = c("0-1", "1", "2", "3","4","5")) +
  scale_size_continuous(limits=c(0,5),
                        labels = c("0-1", "1", "2", "3","4","5")) +
  scale_fill_distiller(limits=c(0,5),
                       labels = c("0-1", "1", "2", "3","4","5"),
                       type = "seq",
                       direction = 1,
                       palette = "Greys") +
  guides(fill= guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)")))), 
         size=guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)"))))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size=18,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 20, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 50,size=18,hjust=1),
        plot.margin = margin(t = 10, r = 20, b = 0, l = 20, unit = "pt")
  )
C

D=ggplot(root, aes(x=longitude, y=latitude)) +
  geom_point(aes(fill= log.pg.mg,
                 size=log.pg.mg),
             shape = 21,color="black") +
  ylim(27.09,27.29) +
  xlim(-98,-97.93) +
  ggtitle("Roots") +
  xlab("") +
  ylab("") +
  scale_alpha_continuous(limits=c(0,5), 
                         labels = c("0-1", "1", "2", "3","4","5")) +
  scale_size_continuous(limits=c(0,5),
                        labels = c("0-1", "1", "2", "3","4","5")) +
  scale_fill_distiller(limits=c(0,5),
                       labels = c("0-1", "1", "2", "3","4","5"),
                       type = "seq",
                       direction = 1,
                       palette = "Greys") +
  guides(fill= guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)")))), 
         size=guide_legend(title=expression(bold(paste('log'[10]," 2-Hydroxyphenylacetic acid (pg/mg)"))))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size=18,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 20, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 50,size=18,hjust=1),
        plot.margin = margin(t = 10, r = 20, b = 0, l = 20, unit = "pt")
  )
D

### --- make spatial concentration plots into multipanel figure
concen=ggarrange(A,B,C,D, common.legend = TRUE, legend = "bottom",
                 labels = c("(a)"), font.label = list(size = 18))
# dev plot
dev.new(
  title = "2HPAA",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
concen
annotate_figure(concen,
                left = text_grob(expression(bold("Latitude")), hjust=-0.5,vjust=3.5, rot=90,size = 18),
                bottom = text_grob(expression(bold("Longitude")), vjust=-4.5, hjust=0.25,size = 18))

### ---------------------------------------------------------- ###
### /// Panel B: 2HPAA tissue concentrations with Bargraph /// ###
### ---------------------------------------------------------- ###

# set factor level order 
dats$type = factor(dats$type,
                    c("soil","leaf", "litter","root"))
# plot it
bar<-ggplot(dats, 
            aes(x=type, y=log.pg.mg)) + 
  geom_boxplot(width=0.5,outlier.shape = NA) +
  geom_point(position=position_jitter(0.15),
             pch=21,size=2.5,color="black",fill="black") + 
  ggtitle("") +
  labs(x="",y="") +
  #geom_text(x=0.75, y=-0.05, label="a",size=5) +
  #geom_text(x=1.0, y=-0.05, label="b",size=5) +
  #geom_text(x=1.25, y=-0.05, label="bc",size=5) + 
  #geom_text(x=1.75, y=-0.05, label="c",size=5) +
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
        plot.margin = margin(t = 10, r = 16, b = 0, l =  0, unit = "pt")
  )
bar

HPAA.bar=ggarrange(bar,labels = c("(b)"), font.label = list(size = 18))
dev.new(
  title = "2HPAA",
  width = 8,
  height = 10,
  noRStudioGD = TRUE
)
HPAA.bar
annotate_figure(HPAA.bar,
                left = text_grob(expression(paste(bold('log'[10]),bold(" 2-Hydroxyphenylacetic acid (pg/mg)"))), 
                                 hjust=0.25,vjust=2, rot=90,size = 18),
                bottom = text_grob(expression(bold("Guinea grass sample type")), vjust=-1, size = 18))



#####################################################################
#####################################################################
### --- FIGURE 2: Guinea grass qualitative chemical profiles  --- ###
#####################################################################
#####################################################################
GGmetab<-load("chemical.classification/Morrison_GuineaGrass_metab_20220210.RData")
head(metab)

### --- break down metabolome by superclass
metab$superclass<-as.factor(metab$superclass)
levels(metab$superclass)

### --- break down metabolome by class
metab$class<-as.factor(metab$class)
levels(metab$class)

### --------------------------------------------------- ###
### --- Panel A. Frequency of Soil Chemical Classes --- ###
### --------------------------------------------------- ###
# subset just soil samples
names(metab)
soilMET<-metab[,c(2,17:21)]
soilMET_comps<-soilMET[soilMET$X8560_4_soil.mzXML.Peak.area != 0, ]
nrow(soilMET_comps) # [1] 4
soil.comps=soilMET_comps$class
# make this into a dataframe of frequency of occurences
soil.class<- as.data.frame(table(unlist(soil.comps)))
# sort by decreasing occurence of compounds 
soil.class=soil.class[order(soil.class$Freq, decreasing = T),]
soil.class$Freq<-as.numeric(soil.class$Freq)
soil.class$percentage<-(soil.class$Freq/sum(soil.class$Freq))*100
# round percentages to one digit
soil.class$percentage<-round(soil.class$percentage, 0)
# change the column names
names(soil.class)=c("class","frequency","percentage")
soil.class

# order the classes by descending frequency
soil.class$class <- factor(soil.class$class, 
                           levels = rev(as.character(soil.class$class)))


# plot it with pie chart 
soil.pie<- ggplot(soil.class[1:3,], aes(x="", y=frequency, fill=class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  ggtitle("(a) Soil") +
  labs(x = NULL,y=NULL,fill="Chemical class") +
  scale_fill_viridis_d(drop = FALSE) +
  geom_text(aes(label = frequency, x=1.3),
            position = position_stack(vjust = 0.5),
            size=4) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position='none',
        legend.title = element_text(size=18,face ="bold"),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 20, face = "bold", 
                                  colour = "black",
                                  vjust=1,hjust=0.5),
        plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = "pt"))
soil.pie

### ------------------------------------------------- ###
### --- Panel B: Frequency of Root Chemical Classes --- ###
### ------------------------------------------------- ###
# subset just root samples
rootMET<-metab[,c(3,17:21)]
rootMET_comps<-rootMET[rootMET$X8560_2_root.mzXML.Peak.area != 0, ]
nrow(rootMET_comps) # [1] 32
root.comps=rootMET_comps$class
# make this into a dataframe of frequency of occurences
root.class<- as.data.frame(table(unlist(root.comps)))
# sort by decreasing occurence of compounds 
root.class=root.class[order(root.class$Freq, decreasing = T),]
root.class$Freq<-as.numeric(root.class$Freq)
root.class$percentage<-(root.class$Freq/sum(root.class$Freq))*100
# round percentages to one digit
root.class$percentage<-round(root.class$percentage, 0)
# change the column names
names(root.class)=c("class","frequency","percentage")
root.class

# order the classes by descending frequency
root.class$class <- factor(root.class$class, 
                           levels = rev(as.character(root.class$class)))
# plot it with pie chart 
root.pie<- ggplot(root.class[1:13,], aes(x="", y=frequency, fill=class))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  ggtitle("(b) Roots") +
  labs(x = NULL,y=NULL,fill="Chemical class") +
  scale_fill_viridis_d(drop = FALSE) +
  geom_text(aes(label = frequency, x=1.3),
            position = position_stack(vjust = 0.5),
            size=4) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position='none',
        legend.title = element_text(size=18,face ="bold"),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 20, face = "bold", 
                                  colour = "black",
                                  vjust=1,hjust=0.5),
        plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = "pt"))
root.pie

### ----------------------------------------------------- ###
### --- Panel C: Frequency of Litter Chemical Classes --- ###
### ----------------------------------------------------- ###
# subset just litter samples
litterMET<-metab[,c(4,17:21)]
litterMET_comps<-litterMET[litterMET$X8560_3_litter.mzXML.Peak.area != 0, ]
nrow(litterMET_comps) # 40
litter.comps=litterMET_comps$class
# make this into a dataframe of frequency of occurences
litter.class<- as.data.frame(table(unlist(litter.comps)))
# sort by decreasing occurence of compounds 
litter.class=litter.class[order(litter.class$Freq, decreasing = T),]
litter.class$Freq<-as.numeric(litter.class$Freq)
litter.class$percentage<-(litter.class$Freq/sum(litter.class$Freq))*100
# round percentages to one digit
litter.class$percentage<-round(litter.class$percentage, 0)
# change the column names
names(litter.class)=c("class","frequency","percentage")
litter.class

# order the classes by descending frequency
litter.class$class <- factor(litter.class$class, 
                             levels = rev(as.character(litter.class$class)))
# plot it with pie chart 
litter.pie<- ggplot(litter.class, aes(x="", y=frequency, fill=class))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  ggtitle("(c) Leaf litter") +
  labs(x = NULL,y=NULL,fill="Chemical class") +
  scale_fill_viridis_d() +
  geom_text(aes(label = frequency, x=1.3),
            position = position_stack(vjust = 0.5),
            size=4) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size=18,face ="bold"),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 20, face = "bold", 
                                  colour = "black",
                                  vjust=1,hjust=0.5),
        plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = "pt"))
litter.pie

### --------------------------------------------------------- ###
### --- Panel D: Frequency of Leaf/Culm  Chemical Classes --- ###
### --------------------------------------------------------- ###
# subset just leaf/culm samples
leafMET<-metab[,c(5,17:21)]
leafMET_comps<-leafMET[leafMET$X8560_1_leaf.mzXML.Peak.area != 0, ]
nrow(leafMET_comps) # 49
leaf.comps=leafMET_comps$class
# make this into a dataframe of frequency of occurences
leaf.class<- as.data.frame(table(unlist(leaf.comps)))
# sort by decreasing occurence of compounds 
leaf.class=leaf.class[order(leaf.class$Freq, decreasing = T),]
leaf.class$Freq<-as.numeric(leaf.class$Freq)
leaf.class$percentage<-(leaf.class$Freq/sum(leaf.class$Freq))*100
# round percentages to one digit
leaf.class$percentage<-round(leaf.class$percentage, 0)
# change the column names
names(leaf.class)=c("class","frequency","percentage")
leaf.class

# order the classes by descending frequency
leaf.class$class <- factor(leaf.class$class, 
                           levels = rev(as.character(leaf.class$class)))
# plot it with pie chart 
leaf.pie<- ggplot(leaf.class[1:16,], aes(x="", y=frequency, fill=class))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  ggtitle("(d) Leaves & culms") +
  labs(x = NULL,y=NULL,fill="Chemical class") +
  scale_fill_viridis_d(drop = FALSE) +
  geom_text(aes(label = frequency, x=1.3),
            position = position_stack(vjust = 0.5),
            size=4) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position='none',
        legend.title = element_text(size=18,face ="bold"),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 20, face = "bold", 
                                  colour = "black",
                                  vjust=1,hjust=0.5),
        plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = "pt"))
leaf.pie

### --- make chemical profile pie charts into multipanel figure
profile=ggarrange(soil.pie,root.pie,litter.pie,leaf.pie, 
                  common.legend = TRUE, legend = "right",
                  nrow=2,ncol=2)
# dev plot
dev.new(
  title = "profile",
  width = 12,
  height = 10
)
profile
annotate_figure(profile,
                bottom = text_grob(expression(bold("Metabolite frequency")), 
                                   vjust=-0.25, hjust= 1.25,size = 18))


####################################################################################
####################################################################################
### --- FIGURE 3: Seedling Survival Across Common Garden Treatment Levels  --- ###
####################################################################################
####################################################################################
surv<-read.csv("recruitment.csv")
surv$percent=surv$percent*100
head(surv)

### -------------------------------------------- ###
### /// Panels A & B: Final Percent Survival /// ###
### -------------------------------------------- ###
### GOAL:
# 1.) Visualize variation in survival with dodgedbargraphs for each species.

final = filter(surv,week==5)
nrow(final) # 24
head(final)

### --- Barplots of Final Percent Survival --- ###
### subset just SHADE treatment
final.shade = filter(final,shade=="shade")
# order the chemical levels in  order of descending survival
final.shade$chemical <- factor(final.shade$chemical,
                               levels = c("none", "standard", "crude"),
                               labels=c("none","2HPAA","whole-plant"))
# order the species in order of descending survival
final.shade$species <- factor(final.shade$species,
                               levels = c("M.maximus", "C.fasciculata","S.scheelei","X.texanum"),
                              labels=c("M. maximus",
                                       "C. fasciculata", 
                                       "S. scheelei",
                                       "X. texanum"))
### compact letter display objects for nice looking text
none_letter_shade <- data.frame(label = c("a", "a", "a","a"), 
                              chemical = factor("none",
                                                levels = c("none",
                                                           "2HPAA",
                                                           "whole-plant")),
                              x = c(0.75, 1.75, 2.75, 3.75),
                              y = c(65,42,60,21))
stand_letter_shade <- data.frame(label = c("b", "a", "a","a"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1, 2, 3, 4),
                               y = c(34,33,64,16))
crude_letter_shade <- data.frame(label = c("b", "b", "b","b"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1.25, 2.25, 3.25, 4.25),
                               y = c(23,19,45,5))

# colors for chemical treatments
surv.colors<-c("white","grey80","grey40")
# plot it
shade.surv<-ggplot(final.shade, aes(x=species, y=percent,fill=chemical)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  scale_fill_manual(values = surv.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("Shade treatment") +
  xlab("") + 
  ylab("") +
  labs(fill = "Chemical treatment") +
  ylim(0,100) + 
  geom_text(data=none_letter_shade, size=8,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_shade, size=8,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_shade, size=8,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
shade.surv

### subset just SUN treatment
final.sun = filter(final,shade=="open")
final.sun$chemical <- factor(final.sun$chemical,
                               levels = c("none", "standard", "crude"),
                             labels=c("none","2HPAA","whole-plant"))
final.sun$species <- factor(final.sun$species,
                              levels = c("M.maximus", "C.fasciculata","S.scheelei","X.texanum"),
                            labels=c("M. maximus",
                                     "C. fasciculata", 
                                     "S. scheelei",
                                     "X. texanum"))
### compact letter display objects for nice looking text
none_letter_sun <- data.frame(label = c("a", "a", "a","a"), 
                            chemical = factor("none",
                                             levels = c("none",
                                                        "2HPAA",
                                                        "whole-plant")),
                            x = c(0.75, 1.75, 2.75, 3.75),
                            y = c(89,55,59,35))
stand_letter_sun <- data.frame(label = c("b", "b", "a","a"), 
                              chemical = factor("none",
                                                levels = c("none",
                                                           "2HPAA",
                                                           "whole-plant")),
                              x = c(1, 2, 3, 4),
                              y = c(54,89,67,33))
crude_letter_sun <- data.frame(label = c("b", "a", "a","b"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1.25, 2.25, 3.25, 4.25),
                               y = c(63,42,54,5))

# plot it
sun.surv<-ggplot(final.sun, aes(x=species, y=percent,fill=chemical)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  scale_fill_manual(values = surv.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("Sun treatment") +
  xlab("") + 
  ylab("Final recruitment (%)") +
  labs(fill = "Chemical treatment") +
  ylim(0,100) + 
  geom_text(data=none_letter_sun, size=8,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_sun, size=8,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_sun, size=8,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
sun.surv


### -------------------------------------------------- ###
### /// Panels C and D: Percent Survival over Time /// ###
### -------------------------------------------------- ###
### GOAL:
# 1.) Visualize survival trends over time with separate time series 
#     plots for plants grown on each treatment for each species.
#   - Scatterplot with lines connecting the data points
surv<-read.csv("recruitment.csv")
surv$percent=surv$percent*100
head(surv)


# subset summary table for shade and sun light levels
surv.shade.t<-filter(surv, shade=='shade')
surv.sun.t<-filter(surv, shade=='open')

### --- Plot line graphs of survival curves for species from the SUN treatment 
time.colors<-c("grey80","grey50","black")

# order factor levels and modify the names
surv.sun.t$chemical <- factor(surv.sun.t$chemical,
                              levels = c("none", "standard", "crude"),
                              labels=c("none","2HPAA","whole-plant"))
surv.sun.t$species <- factor(surv.sun.t$species,
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
                          x = c(5.75, 5.75, 5.75),
                          y = c(83, 57, 48))
Cf_text_sun <- data.frame(label = c("a", "b", "b"), 
                          species = factor("C. fasciculata",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(5.75, 5.75, 5.75),
                          y = c(83, 50, 37))
Ss_text_sun <- data.frame(label = c("a", "a", "a"), 
                          species = factor("S. scheelei",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(5.75, 5.75, 5.75),
                          y = c(55, 62, 49))
Xt_text_sun <- data.frame(label = c("a", "a", "b"), 
                          species = factor("X. texanum",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(5.75, 5.75, 5.75),
                          y = c(32.5, 27.75, 0))

#
sun.time<-ggplot(surv.sun.t, aes(x=week, y=percent)) +
  geom_line(aes(color=chemical),position=position_dodge(0.3),lwd=1) +
  geom_point(aes(color=chemical),position=position_dodge(0.1),
             size=3, shape=21, fill="white") +
  facet_grid(~ species) +
  scale_color_manual(values=time.colors) +
  ylim(0,100) +
  ylab('Recruitment over time (%)') +
  xlab('') +
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
  scale_x_discrete(breaks=c("0","1","2","3","4","5"),
                   limits=c("0","1","2","3","4","5")) +
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
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18,face="bold",vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face="bold",vjust=1),
        plot.margin = margin(t = 10, r = 5, b = 0, l = 0, unit = "pt"))
sun.time

### --- Plot line graphs of survival curves for species from the SHADE treatment 
surv.shade.t$chemical <- factor(surv.shade.t$chemical,
                                levels = c("none", "standard", "crude"),
                                labels=c("none","2HPAA","whole-plant"))
surv.shade.t$species <- factor(surv.shade.t$species,
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
                            x = c(5.75, 5.75, 5.75),
                            y = c(59.25, 18, 28.75))
Cf_text_shade <- data.frame(label = c("a", "a", "a"), 
                            species = factor("C. fasciculata",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(5.75, 5.75, 5.75),
                            y = c(29, 35.5, 13.5))
Ss_text_shade <- data.frame(label = c("a", "a", "b"), 
                            species = factor("S. scheelei",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(5.75, 5.75, 5.75),
                            y = c(55, 60, 39))
Xt_text_shade <- data.frame(label = c("a", "a", "b"), 
                            species = factor("X. texanum",
                                             levels = c("M. maximus",
                                                        "C. fasciculata",
                                                        "S. scheelei",
                                                        "X. texanum")),
                            x = c(5.75, 5.75, 5.75),
                            y = c(15.75, 10.5, 0))
# 
shade.time<-ggplot(surv.shade.t, aes(x=week, y=percent)) +
  geom_line(aes(color=chemical),position=position_dodge(0.3),lwd=1) +
  geom_point(aes(color=chemical),position=position_dodge(0.1),
             size=3, shape=21, fill="white") +
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
  ggtitle('') +
  labs(color='Chemical treatment') +
  scale_x_discrete(breaks=c("0","1","2","3","4","5"),
                   limits=c("0","1","2","3","4","5")) +
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
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18,face="bold",vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face="bold",vjust=1),
        plot.margin = margin(t = 10, r = 5, b = 0, l = 0, unit = "pt"))
shade.time


### --- make survival over time  plots into multipanel figure
time=ggarrange(sun.surv,shade.surv,
               sun.time,shade.time,
               labels = c("(a)", "(b)",
                          "(c)","(d)"),
               font.label = list(size = 18),
               common.legend = TRUE, legend = "bottom",
               nrow=2,ncol=2)
# dev plot
dev.new(
  title = "time.survival",
  width = 20,
  height = 15,
  noRStudioGD = TRUE
)
time
annotate_figure(time,
                top = text_grob("Species", vjust=38,hjust=-0.2, size = 18,face='bold'),
                bottom = text_grob("Week", vjust=-4.25,hjust=-0.2, size = 18,face='bold'))



##########################################################
##########################################################
### --- FIGURE 4: Growth - Shoot Height & BIOMASS  --- ###
##########################################################
##########################################################

### ---------------------------------------- ###
### /// Panels A and B - Seedling Height /// ###
### ---------------------------------------- ###
### GOAL:
# 1.) Visualize variation in final plant height with boxplots for each species.

grow<-read.csv("growth2.csv") # this version has no 0s as place holders for X.tex on crude treatment
grow=grow[,-c(4,6,7)]
head(grow,10)

### melt the dataframe 
grow.m <- melt(grow, id = c("species","chemical","shade","week"))
grow.m=grow.m[,-c(5)]
head(grow.m,50)
nrow(grow.m) # 9016
### subset just the final measurements
grow.final = filter(grow.m,week==4)
grow.final$value<-as.numeric(grow.final$value)
nrow(grow.final) # 2352


### --- GG jittered boxplot of Final Growth  in SHADE--- ###
### subset just the SUN treatment
grow.sun = filter(grow.final,shade=='open')
# order the chemicals in order of increasing severity
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

### compact letter display objects for nice looking text
none_letter_sun_a <- data.frame(label = c("a", "a", "a","a"), 
                              chemical = factor("none",
                                                levels = c("none",
                                                           "2HPAA",
                                                           "whole-plant")),
                              x = c(0.75, 1.75, 2.75, 3.75),
                              y = c(42,25,22,7))
stand_letter_sun_a <- data.frame(label = c("a", "b", "a","b"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1, 2, 3, 4),
                               y = c(42,25,22,7))
crude_letter_sun_a <- data.frame(label = c("b", "a", "a","c"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1.25, 2.25, 3.25, 4.25),
                               y = c(42,25,22,7))
# plot it
sun.grow<-ggplot(grow.sun, 
                 aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = grow.colors) +
  geom_hline(yintercept = 0,lwd=0.3, color = 'grey50') +
  ggtitle("Sun treatment") +
  labs(x="",y="Shoot height (cm)",fill = "Chemical treatment") +
  ylim(0,43) +
  geom_text(data=none_letter_sun_a, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_sun_a, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_sun_a, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major.y = element_line(color="grey50"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 18, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = -3, r = 5, b = -20, l = -10, unit = "pt"))
sun.grow


### --- GG jittered boxplot of Final Growth  in SHADE--- ###
### subset just the SHADE treatment
grow.shade = filter(grow.final,shade=='shade')
# order the chemicals in order of increasing severity
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
### compact letter display objects for nice looking text
none_letter_shade_b <- data.frame(label = c("a", "a", "a","a"), 
                                chemical = factor("none",
                                                  levels = c("none",
                                                             "2HPAA",
                                                             "whole-plant")),
                                x = c(0.75, 1.75, 2.75, 3.75),
                                y = c(23.5,16.5,13,4))
stand_letter_shade_b <- data.frame(label = c("b", "a", "b","b"), 
                                 chemical = factor("none",
                                                   levels = c("none",
                                                              "2HPAA",
                                                              "whole-plant")),
                                 x = c(1, 2, 3, 4),
                                 y = c(23.5,16.5,13,4))
crude_letter_shade_b <- data.frame(label = c("b", "b", "c","c"), 
                                 chemical = factor("none",
                                                   levels = c("none",
                                                              "2HPAA",
                                                              "whole-plant")),
                                 x = c(1.25, 2.25, 3.25, 4.25),
                                 y = c(23.5,16.5,13,4))
### plot it
shade.grow<-ggplot(grow.shade, 
                   aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = grow.colors) +
  geom_hline(yintercept = 0,lwd=0.3, color = 'grey50') +
  ggtitle("Shade treatment") +
  labs(x="",y="",fill = "Chemical treatment") +
  ylim(0,25) +
  geom_text(data=none_letter_shade_b, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_shade_b, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_shade_b, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major.y = element_line(color="grey50"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 18, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = -3, r = 5, b = -20, l = -10, unit = "pt"))
shade.grow


### -------------------------------------- ###
### /// Panels C and D - Total Biomass /// ###
### -------------------------------------- ###
biomass<-read.csv("biomass.csv")
head(biomass,10)
### GOAL:
# 1.) Visualize variation in final total (above+blow)  biomass with boxplots for each species.

### subset out just below and above ground mass
total = filter(biomass,tissue =='total')
total=total[,-c(4:10)]
head(total,10)

### melt the dataframe so that each row is an individual plant that survived or died (= 1 or 0)
total.m <- melt(total, id = c("species","chemical","shade"))
total.m=total.m[,-c(4)]
head(total.m,50)
nrow(total.m) # 2352


### --- GG jittered boxplot of Final Biomass  in SUN--- ###
### subset just sun treatment
sun.total=filter(total.m, shade == 'open')
# colors for chemical treatments
mass.colors<-c("white","grey80","grey40")
sun.total$chemical <- factor(sun.total$chemical,
                            levels = c("none", "standard", "crude"),
                            labels=c("none","2HPAA","whole-plant"))
sun.total$species <- factor(sun.total$species,
                           levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                           labels = c("M. maximus",
                                 "C. fasciculata", 
                                 "S. scheelei",
                                 "X. texanum"))
### compact letter display objects for nice looking text
none_letter_sun_c <- data.frame(label = c("a", "a", "a","a"), 
                              chemical = factor("none",
                                                levels = c("none",
                                                           "2HPAA",
                                                           "whole-plant")),
                              x = c(0.75, 1.75, 2.75, 3.75),
                              y = c(450,240,130,45))
stand_letter_sun_c <- data.frame(label = c("a", "b", "a","a"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1, 2, 3, 4),
                               y = c(450,240,130,45))
crude_letter_sun_c <- data.frame(label = c("b", "a", "a","b"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1.25, 2.25, 3.25, 4.25),
                               y = c(450,240,130,45))
### plot it 
sun.mass<-ggplot(sun.total, 
                 aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = mass.colors) +
  geom_hline(yintercept = 0,lwd=0.3, color = 'grey50') +
  ggtitle("") +
  labs(x="",y="Total biomass (mg)",fill = "Chemical treatment") +
  ylim(0,450) +
  geom_text(data=none_letter_sun_c, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_sun_c, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_sun_c, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major.y = element_line(color="grey50"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 18, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = -3, r = 5, b = -20, l = -10, unit = "pt"))
sun.mass


### --- GG jittered boxplot of Final Biomass  in SHADE--- ###
### subset just SHADE treatment
shade.total=filter(total.m, shade == 'shade')
# colors for chemical treatments
mass.colors<-c("white","grey80","grey40")

shade.total$chemical <- factor(shade.total$chemical,
                              levels = c("none", "standard", "crude"),
                              labels=c("none","2HPAA","whole-plant"))
shade.total$species <- factor(shade.total$species,
                             levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                             labels = c("M. maximus",
                                 "C. fasciculata", 
                                 "S. scheelei",
                                 "X. texanum"))
### compact letter display objects for nice looking text
none_letter_shade_d <- data.frame(label = c("a", "a", "a","a"), 
                              chemical = factor("none",
                                                levels = c("none",
                                                           "2HPAA",
                                                           "whole-plant")),
                              x = c(0.75, 1.75, 2.75, 3.75),
                              y = c(36,19,9,7))
stand_letter_shade_d <- data.frame(label = c("a", "a", "a","a"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1, 2, 3, 4),
                               y = c(36,19,9,7))
crude_letter_shade_d <- data.frame(label = c("a", "a", "b","b"), 
                               chemical = factor("none",
                                                 levels = c("none",
                                                            "2HPAA",
                                                            "whole-plant")),
                               x = c(1.25, 2.25, 3.25, 4.25),
                               y = c(36,19,9,7))
### plot it 
shade.mass<-ggplot(shade.total, 
                   aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = mass.colors) +
  geom_hline(yintercept = 0,lwd=0.3, color = 'grey50') +
  ggtitle("") +
  labs(x="",y="",fill = "Chemical treatment") +
  ylim(0,38) +
  geom_text(data=none_letter_shade_d, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_shade_d, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_shade_d, size=7,fontface="bold",
            mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major.y = element_line(color="grey50"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 18, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = -3, r = 5, b = -20, l = -10, unit = "pt"))
shade.mass


### --- make sun/shade growth plots into multipanel figure
growth=ggarrange(sun.grow,shade.grow,sun.mass,shade.mass,
                 labels = c("(a)", "(b)","(c)","(d)"),
                 font.label = list(size = 18),
                 common.legend = TRUE, legend = "bottom")
# dev plot
dev.new(
  title = "growth",
  width = 12,
  height = 12,
  noRStudioGD = TRUE
)
growth
annotate_figure(growth,
                bottom = text_grob("Species", vjust=-4, size = 18,face='bold'))



### ----------------------------------------- ###
### /// FIGURE S4: Final Root-Shoot Ratio /// ###
### ----------------------------------------- ###
biomass<-read.csv("biomass.csv")

### GOAL:
# 1.) Visualize root-shoot ratios for each species by treatment with box-whisker plots. 
### subset out just below and above ground mass
ratio = filter(biomass,tissue=='root.shoot' )
ratio2 = ratio[,-c(4:10)]


### subset just SUN treatment
sun.ratio=filter(ratio2, shade == 'open')
sun.ratio2=sun.ratio[,-c(3)]
### melt the dataframe so that each root:shoot value is a row
sun.ratio.m <- melt(sun.ratio2, id = c("species","chemical"))
sun.ratio.m=na.omit(sun.ratio.m)
nrow(sun.ratio.m) # [1] 511

# colors for chemical treatments
ratio.colors<-c("white","grey80","grey40")
sun.ratio.m$chemical <- factor(sun.ratio.m$chemical,
                               levels = c("none", "standard", "crude"),
                               labels=c("none","2HPAA","whole-plant"))
sun.ratio.m$species <- factor(sun.ratio.m$species,
                              levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                              labels = c("M. maximus",
                                         "C. fasciculata", 
                                         "S. scheelei",
                                         "X. texanum"))
### compact letter display objects for nice looking text
none_letter_sun_a <- data.frame(label = c("a", "a", "a","a"), 
                                chemical = factor("none",
                                                  levels = c("none",
                                                             "2HPAA",
                                                             "whole-plant")),
                                x = c(0.75, 1.75, 2.75, 3.75),
                                y = c(2.6,2.5,2.55,1.7))
stand_letter_sun_a <- data.frame(label = c("b", "b", "a","a"), 
                                 chemical = factor("none",
                                                   levels = c("none",
                                                              "2HPAA",
                                                              "whole-plant")),
                                 x = c(1, 2, 3, 4),
                                 y = c(2.6,2.5,2.55,1.7))
crude_letter_sun_a <- data.frame(label = c("b", "a", "a","b"), 
                                 chemical = factor("none",
                                                   levels = c("none",
                                                              "2HPAA",
                                                              "whole-plant")),
                                 x = c(1.25, 2.25, 3.25, 4.25),
                                 y = c(2.6,2.5,2.55,1.7))
### plot it 
rs.sun<-ggplot(sun.ratio.m, 
               aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = grow.colors) +
  geom_hline(yintercept = 0,lwd=0.3, color = 'grey50') +
  ggtitle("Sun treatment") +
  labs(x="",y="Root/shoot ratio",fill = "Chemical treatment") +
  ylim(0,2.6) +
  geom_text(data=none_letter_sun_a, size=7,fontface="bold",
           mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=stand_letter_sun_a, size=7,fontface="bold",
           mapping = aes(x = x, y = y,label=label)) +
  geom_text(data=crude_letter_sun_a, size=7,fontface="bold",
           mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major.y = element_line(color="grey50"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 18, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = -3, r = 5, b = -20, l = -10, unit = "pt"))
rs.sun



### subset just SHADE treatment
shade.ratio=filter(ratio2, shade == 'shade')[,-3]
### melt the dataframe so that each root:shoot value is a row
shade.ratio.m <- melt(shade.ratio, id = c("species","chemical"))[,-3]
shade.ratio.m=na.omit(shade.ratio.m)
nrow(shade.ratio.m) # 307
  
# colors for chemical treatments
ratio.colors<-c("white","grey80","grey40")
shade.ratio.m$chemical <- factor(shade.ratio.m$chemical,
                       levels = c("none", "standard", "crude"),
                       labels=c("none","2HPAA","whole-plant"))
shade.ratio.m$species <- factor(shade.ratio.m$species,
                      levels = c("M.maximus","C.fasciculata", "S.scheelei","X.texanum"),
                      labels = c("M. maximus",
                                 "C. fasciculata", 
                                 "S. scheelei",
                                 "X. texanum"))

### compact letter display objects for nice looking text
none_letter_shade_b <- data.frame(label = c("a", "a", "a","a"), 
                                  chemical = factor("none",
                                                    levels = c("none",
                                                               "2HPAA",
                                                               "whole-plant")),
                                  x = c(0.75, 1.75, 2.75, 3.75),
                                  y = c(2.1,1.15,2.4,1.85))
stand_letter_shade_b <- data.frame(label = c("b", "a", "a","b"), 
                                   chemical = factor("none",
                                                     levels = c("none",
                                                                "2HPAA",
                                                                "whole-plant")),
                                   x = c(1, 2, 3, 4),
                                   y = c(2.1,1.15,2.4,1.85))
crude_letter_shade_b <- data.frame(label = c("b", "a", "b","c"), 
                                   chemical = factor("none",
                                                     levels = c("none",
                                                                "2HPAA",
                                                                "whole-plant")),
                                   x = c(1.25, 2.25, 3.25, 4.25),
                                   y = c(2.1,1.15,2.4,1.85))
### plot it 
rs.shade<-ggplot(shade.ratio.m, 
                   aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = grow.colors) +
  geom_hline(yintercept = 0,lwd=0.3, color = 'grey50') +
  ggtitle("Shade treatment") +
  labs(x="",y="",fill = "Chemical treatment") +
  ylim(0,2.5) +
 geom_text(data=none_letter_shade_b, size=7,fontface="bold",
           mapping = aes(x = x, y = y,label=label)) +
 geom_text(data=stand_letter_shade_b, size=7,fontface="bold",
           mapping = aes(x = x, y = y,label=label)) +
 geom_text(data=crude_letter_shade_b, size=7,fontface="bold",
           mapping = aes(x = x, y = y,label=label)) +
  theme(panel.grid.major.y = element_line(color="grey50"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 18, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=18),
        axis.text.x = element_text(angle = 50,hjust=1,size=18, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=18),
        plot.margin = margin(t = -3, r = 5, b = -20, l = -10, unit = "pt"))
rs.shade


### --- make sun/shade growth plots into multipanel figure
biom=ggarrange(rs.sun,rs.shade,
               nrow=1,ncol=2,
               labels = c("(a)", "(b)"),
               font.label = list(size = 18),
               common.legend = TRUE, legend = "bottom")
# dev plot
dev.new(
  title = "root-shoot",
  width = 20,
  height = 15,
  noRStudioGD = TRUE
)
biom
annotate_figure(biom,
                bottom = text_grob("Species", vjust=-5.5, size = 16,face="bold"))


############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################