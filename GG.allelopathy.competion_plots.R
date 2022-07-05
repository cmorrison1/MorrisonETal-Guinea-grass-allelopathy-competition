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

#########################################################################################
#########################################################################################
#### --- FIGURE 1: 2-Hydroxyphenylacetic Acid Concentration Rancho La Paloma, TX --- ####
#########################################################################################
#########################################################################################
### GOALS:
# 1.) Visualize 2-HPAA concentrations with bargraph 

### ------------------------------------------------------------------------- ###
### /// Plot tissue type 2HPAA concentration correlations per individual  /// ###
### ------------------------------------------------------------------------- ###
chem<-read.csv('LP_2HPAA_correlates.csv')
chem

### just site, type and concentration for community df
lp.chem<-chem[,c(3,6:7)]
### log transform concentration values 
lp.chem$log.pg.mg<-log(lp.chem$pg.mg + 1)
#'cast' the data a molten dataframe using dcast()
lp.chem.w <-dcast(
  lp.chem, 
  site ~ type,
  value.var="log.pg.mg",
  fun.aggregate=sum
)

row.names(lp.chem.w)=lp.chem.w$site
lp.chem.w=lp.chem.w[,-1]
lp.chem.w [lp.chem.w  == 0] <- NA
lp.chem.w 

### --- plot multiple correlations
# with psych package
mult<-cor(lp.chem.w[, c('leaf', 'litter', 'root', 'soil')],
          method=c('spearman'),use = "complete.obs")
mult
pairs.panels(mult)

dev.new(
  title = "correlations",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
par(col="grey50")
par(pch=".")
pairs(
  lp.chem.w,
  lower.panel = panel.smooth,
  upper.panel = panel.cor,
  diag.panel = panel.hist,
)
### with GGally and ggpairs()
lp.chem.w2$chemical <- factor(final.shade$chemical,
                               levels = c("none", "standard", "crude"),
                               labels=c("none","2HPAA","whole plant"))

dev.new(
  title = "correlations",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)

ggpairs(lp.chem.w,
        upper = list(
  continuous = wrap('cor', method = "spearman"))) + 
  labs(x='Log10 2HPAA concentration', y='Log10 2HPAA concentration') + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=12),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=12),
        plot.margin = margin(t = 10, r = 10, b = 0, l =  0, unit = "pt")
        )


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
        legend.title = element_text(size=16,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=16),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 16, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
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
        legend.title = element_text(size=16,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=16),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 16, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
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
        legend.title = element_text(size=16,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=16),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 16, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
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
        legend.title = element_text(size=16,face='bold'),
        legend.key=element_blank(),
        legend.text = element_text(size=16),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        plot.title = element_text(vjust=1,hjust=0.5,size = 16, face = "bold", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 10, l = 10),size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
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
                left = text_grob(expression(bold("Latitude")), hjust=-0.5,vjust=3.5, rot=90,size = 16),
                bottom = text_grob(expression(bold("Longitude")), vjust=-4.5, hjust=0.25,size = 16))

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
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face ="bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 70,hjust=1,size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 0, b = 0, l =  0, unit = "pt")
  )
bar

HPAA.bar=ggarrange(bar,labels = c("(b)"), font.label = list(size = 18))
dev.new(
  title = "2HPAA",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
HPAA.bar
annotate_figure(HPAA.bar,
                left = text_grob(expression(paste(bold('log'[10]),bold(" 2-Hydroxyphenylacetic acid (pg/mg)"))), 
                                 hjust=0.25,vjust=2, rot=90,size = 16),
                bottom = text_grob(expression(bold("Guinea grass sample type")), vjust=-1, size = 16))



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

### ------------------------------------------------- ###
### --- Panel A. Frequency of Soil Chemical Classes --- ###
### ------------------------------------------------- ###
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
        legend.title = element_text(size=16,face ="bold"),
        legend.text = element_text(size=11),
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
  width = 10,
  height = 10
)
profile
annotate_figure(profile,
                bottom = text_grob("Metabolite frequency", 
                                   vjust=-0.25, hjust= 1.0,size = 16))


####################################################################################
####################################################################################
### --- FIGURE 3: Seedling Survival Across Common Garden Treatment Levels  --- ###
####################################################################################
####################################################################################
surv<-read.csv("germination.csv")
head(surv)

### -------------------------------------------- ###
### /// Panels A & B: Final Percent Survival /// ###
### -------------------------------------------- ###
### GOAL:
# 1.) Visualize variation in survival with bargraphs for each species.
#   - Plot survival across treatments for each species with dodged bargraphs.

final = filter(surv,week==5)
nrow(final) # 24
head(final)
final$percent=final$percent*100

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
  geom_text(x=0.75, y=64, label="a",size=5) +
  geom_text(x=1.0, y=33, label="b",size=5) +
  geom_text(x=1.25, y=22, label="b",size=5) +
  geom_text(x=1.75, y=41, label="a",size=5) +
  geom_text(x=2.05, y=32, label="a",size=5) +
  geom_text(x=2.25, y=18, label="b",size=5) +
  geom_text(x=2.73, y=59, label="ab",size=5) +
  geom_text(x=3.0, y=63, label="a",size=5) +
  geom_text(x=3.25, y=44, label="b",size=5) +
  geom_text(x=3.75, y=20, label="a",size=5) +
  geom_text(x=4.0, y=15, label="a",size=5) +
  geom_text(x=4.25, y=4, label="b",size=5) +
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
  geom_text(x=0.75, y=88, label="a",size=5) +
  geom_text(x=1.0, y=53, label="b",size=5) +
  geom_text(x=1.25, y=62, label="b",size=5) + 
  geom_text(x=1.75, y=54, label="a",size=5) +
  geom_text(x=2.0, y=88, label="b",size=5) +
  geom_text(x=2.25, y=41, label="a",size=5) +
  geom_text(x=2.73, y=58, label="a",size=5) +
  geom_text(x=3.0, y=66, label="a",size=5) +
  geom_text(x=3.25, y=53, label="a",size=5) +
  geom_text(x=3.75, y=34, label="a",size=5) +
  geom_text(x=4.0, y=32, label="a",size=5) +
  geom_text(x=4.25, y=4, label="b",size=5) +
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
sun.surv


### -------------------------------------------------- ###
### /// Panels C and D: Percent Survival over Time /// ###
### -------------------------------------------------- ###
### GOAL:
# 1.) Visualize survival trends over time with separate time series 
#     plots for plants grown on each treatment for each species.
#   - Scatterplot with lines connecting the data points
surv<-read.csv("survival.csv")
surv=surv[,-c(4)]
head(surv)

### melt the dataframe so that each row is an individual plant that survived or died (= 1 or 0)
surv.m <- melt(surv, id = c("species","chemical","shade","week"))
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
                      x = c(5.5, 5.5, 5.5),
                      y = c(82, 56, 48))
Cf_text_sun <- data.frame(label = c("a", "b", "b"), 
                      species = factor("C. fasciculata",
                                       levels = c("M. maximus",
                                                  "C. fasciculata",
                                                  "S. scheelei",
                                                  "X. texanum")),
                      x = c(5.5, 5.5, 5.5),
                      y = c(82, 49, 36))
Ss_text_sun <- data.frame(label = c("a", "a", "a"), 
                      species = factor("S. scheelei",
                                       levels = c("M. maximus",
                                                  "C. fasciculata",
                                                  "S. scheelei",
                                                  "X. texanum")),
                      x = c(5.5, 5.5, 5.5),
                      y = c(53, 61, 48))
Xt_text_sun <- data.frame(label = c("a", "a", "b"), 
                      species = factor("X. texanum",
                                       levels = c("M. maximus",
                                                  "C. fasciculata",
                                                  "S. scheelei",
                                                  "X. texanum")),
                      x = c(5.5, 5.5, 5.5),
                      y = c(31.25, 26.5, 0))

#
sun.time<-ggplot(surv.sun, aes(x=week, y=value)) +
  geom_line(aes(color=chemical),position=position_dodge(0.3),lwd=1) +
  geom_point(aes(color=chemical),position=position_dodge(0.1),
             size=3, shape=21, fill="white") +
  geom_errorbar(aes(color=chemical,ymin=value-se, ymax=value+se), 
                width=0.9, position=position_dodge(0.1),lwd=1) + 
  facet_grid(~ species) +
  scale_color_manual(values=time.colors) +
  ylim(0,100) +
  ylab('Recruitment over time (%)') +
  xlab('') +
  scale_x_discrete(breaks=c("0","1","2","3","4","5"),
                   limits=c("0","1","2","3","4","5")) +
  geom_text(data = Mm_text_sun, mapping = aes(x = x, y = y,
                                          label = label),size=5) +
  geom_text(data = Cf_text_sun, mapping = aes(x = x, y = y,
                                          label = label),size=5) +
  geom_text(data = Ss_text_sun, mapping = aes(x = x, y = y,
                                          label = label),size=5) +
  geom_text(data = Xt_text_sun, mapping = aes(x = x, y = y,
                                          label = label),size=5) +
  ggtitle('') +
  labs(color='Chemical treatment') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color="black"),
        panel.background = element_blank(),
        strip.text = element_text(face = "italic",size=14),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"),panel.spacing = unit(2, "lines"),
        plot.title = element_text(vjust=2,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black",fill='white'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=16),
        axis.text.x = element_text(hjust=1,size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,face="bold",vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16,face="bold",vjust=1),
        plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))
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
                          x = c(5.5, 5.5, 5.5),
                          y = c(58, 17.5, 28.75))
Cf_text_shade <- data.frame(label = c("a", "a", "a"), 
                          species = factor("C. fasciculata",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(5.5, 5.5, 5.5),
                          y = c(28.25, 35, 13.5))
Ss_text_shade <- data.frame(label = c("a", "a", "b"), 
                          species = factor("S. scheelei",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(5.5, 5.5, 5.5),
                          y = c(53.25, 58.25, 39))
Xt_text_shade <- data.frame(label = c("a", "a", "b"), 
                          species = factor("X. texanum",
                                           levels = c("M. maximus",
                                                      "C. fasciculata",
                                                      "S. scheelei",
                                                      "X. texanum")),
                          x = c(5.5, 5.5, 5.5),
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
                                                label = label),size=5) +
  geom_text(data = Cf_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=5) +
  geom_text(data = Ss_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=5) +
  geom_text(data = Xt_text_shade, mapping = aes(x = x, y = y,
                                                label = label),size=5) +
  scale_x_discrete(breaks=c("0","1","2","3","4","5"),
                   limits=c("0","1","2","3","4","5")) +
  ggtitle('') +
  labs(color='Chemical treatment') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color="black"),
        panel.background = element_blank(),
        strip.text = element_text(face = "italic",size=14),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"
        ),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(vjust=2,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black",fill='white'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=16),
        axis.text.x = element_text(hjust=1,size=16),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,face="bold",vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16,face='bold'),
        plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))
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
                top = text_grob("Species", vjust=42.5,hjust=-0.2, size = 16,face='bold'),
                bottom = text_grob("Week", vjust=-4.5,hjust=-0.2, size = 16,face='bold'))



##########################################################
##########################################################
### --- FIGURE 4: Relative Shoot Growth & BIOMASS  --- ###
##########################################################
##########################################################

### ---------------------------------------- ###
### /// Panels A and B - Relative Growth /// ###
### ---------------------------------------- ###
### GOAL:
# 1.) Visualize variation in final plant height  with dofged bargraphs for each species.
grow<-read.csv("growth.csv")
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


### --- make sun/shade growth plots into multipanel figure
growth=ggarrange(sun.grow,shade.grow,sun.mass,shade.mass,
                 labels = c("(a)", "(b)","(c)","(d)"),
                 font.label = list(size = 18),
                 common.legend = TRUE, legend = "bottom")
# dev plot
dev.new(
  title = "growth",
  width = 12,
  height = 10,
  noRStudioGD = TRUE
)
growth
annotate_figure(growth,
                bottom = text_grob("Species", vjust=-5.5, size = 16,face='bold'))






### ----------------------------------------- ###
### /// FIGURE S4: Final Root-Shoot Ratio /// ###
### ----------------------------------------- ###
biomass<-read.csv("biomass.csv")

### GOAL:
# 1.) Visualize root-shoot ratios for each species by treatment with box-whisker plots. 
### subset out just below and above ground mass
ratio = filter(biomass,tissue=='root.shoot' )
ratio2 = ratio[,-c(4:10)]

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
### --- Stacked boxplots with multiple groups
rs.shade<-ggplot(shade.ratio.m, 
                 aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size=1.25) + 
  scale_fill_manual(values = ratio.colors) +
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("") +
  labs(x="",y="",fill = "Chemical treatment") +
  ylim(0,2.5) +
  geom_text(x=0.75, y=-0.05, label="a",size=5) +
  geom_text(x=1.0, y=-0.05, label="b",size=5) +
  geom_text(x=1.25, y=-0.05, label="b",size=5) + 
  geom_text(x=1.75, y=-0.05, label="a",size=5) +
  geom_text(x=2.0, y=-0.05, label="a",size=5) +
  geom_text(x=2.25, y=-0.05, label="a",size=5) +
  geom_text(x=2.73, y=-0.05, label="a",size=5) +
  geom_text(x=3.0, y=-0.05, label="a",size=5) +
  geom_text(x=3.25, y=-0.05, label="b",size=5) +
  geom_text(x=3.75, y=-0.05, label="a",size=5) +
  geom_text(x=4.0, y=-0.05, label="b",size=5) +
  geom_text(x=4.25, y=-0.05, label="c",size=5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=14),
        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
rs.shade

#
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

### --- Stacked barplots with multiple groups
rs.sun<-ggplot(sun.ratio.m, 
               aes(x=species, y=value,fill = chemical)) + 
  geom_boxplot(position=position_dodge(width=0.75),
               outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             pch=21,size = 1.25) + 
  scale_fill_manual(values = ratio.colors) + 
  geom_hline(yintercept = 0,lwd=0.75) +
  ggtitle("") +
  labs(x="",y="Root/Shoot ratio",fill = "Chemical treatment") +
  ylim(0,2.5) +
  geom_text(x=0.75, y=-0.05, label="a",size=5) +
  geom_text(x=1.0, y=-0.05, label="b",size=5) +
  geom_text(x=1.25, y=-0.05, label="b",size=5) + 
  geom_text(x=1.75, y=-0.05, label="a",size=5) +
  geom_text(x=2.0, y=-0.05, label="b",size=5) +
  geom_text(x=2.25, y=-0.05, label="a",size=5) +
  geom_text(x=2.73, y=-0.05, label="a",size=5) +
  geom_text(x=3.0, y=-0.05, label="a",size=5) +
  geom_text(x=3.25, y=-0.05, label="a",size=5) +
  geom_text(x=3.75, y=-0.05, label="a",size=5) +
  geom_text(x=4.0, y=-0.05, label="a",size=5) +
  geom_text(x=4.25, y=-0.05, label="b",size=5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(color="black"),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=14),
        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16, face='bold',vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
rs.sun

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





############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################