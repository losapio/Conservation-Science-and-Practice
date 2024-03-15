library(reshape2)
library(effects)
library(vegan)
library(BiodiversityR)
library(FD)
library(norm)
library(pdiv)
library(igraph)
library(reshape2)
library(rms)
library(lme4)
library(multcomp)
library(ggplot2)
library(bipartite)
library(nnet)
library(effects)
library(car)
library(akima)
library(plotrix)
library(plyr)
library(betapart)
library(nlsem)
library(car)
library(lme4)
library(emmeans)
library(lmerTest)
library(lavaan)
library(semPlot)
library(nlsem)
library(rnetcarto)

setwd("/Users/gianalberto/pubblicazioni/lesvos_plantdiv")

###################################
#### alpha diversity###### 

survey <- read.table("plant_survey_2018final.csv",	stringsAsFactors=TRUE, header=TRUE, dec=".", sep=",")
head(survey)
survey$enclosure <- relevel(survey$enclosure, ref="out")
survey$plot <- relevel(survey$plot, ref="open")

#
nsp<-ncol(survey)-5
for (i in 1:nsp){survey[,i+4]<-as.numeric(survey[,i+4])}

# species accumulation curve

sac<-specaccum(survey[,6:ncol(survey)])

setEPS()
postscript("sac.eps", widt=3, height=3)
plot(sac, ci.type="polygon", ci.col="grey60")
dev.off()

#
survey$alphadiv<-0
survey$sprich<-0
for(i in 1:nrow(survey)){
  survey$alphadiv[i]<-vegan::diversity(survey[i,6:ncol(survey)], index="shannon")
  survey$sprich[i]<-vegan::specnumber(survey[i,6:ncol(survey)])
}

### stats ###
mod.alphadiv <- lmer(alphadiv ~ enclosure*plot + (1|site/blockn),  data=survey)

Anova(mod.alphadiv)
summary(mod.alphadiv)

emmeans(mod.alphadiv, pairwise ~ enclosure)
(2.16-2.05)/2.05*100
(2.32-2.05)/2.05*100

emmeans(mod.alphadiv, pairwise ~ plot)
(2.32-2.00)/2.00*100
(2.20-2.00)/2.00*100

emmeans(mod.alphadiv, pairwise ~ enclosure:plot)
(2.23-1.66)/1.66*100

write.csv(summary(mod.alphadiv)$coefficients, "alphadiv.csv")

setEPS()
postscript("fig_mod.alphadiv.eps", widt=3.5, height=3.5)
plot(allEffects(mod.alphadiv))
dev.off()

###################################################################
################ beta diversity ######################

### calculating beta diversity

betadiv.in.open  <-vegdist(survey[survey$land=="in"&survey$plot=="open",6:ncol(survey)], binary=TRUE)
betadiv.in.sarco <-vegdist(survey[survey$land=="in"&survey$plot=="sarco",6:ncol(survey)], binary=TRUE)
betadiv.in.remo  <-vegdist(survey[survey$land=="in"&survey$plot=="remo",6:ncol(survey)], binary=TRUE)

betadiv.out.out.open  <-vegdist(survey[survey$land=="out"&survey$enclosure=="out"&survey$plot=="open",6:ncol(survey)], binary=TRUE)
betadiv.out.out.sarco <-vegdist(survey[survey$land=="out"&survey$enclosure=="out"&survey$plot=="sarco",6:ncol(survey)], binary=TRUE)
betadiv.out.out.remo  <-vegdist(survey[survey$land=="out"&survey$enclosure=="out"&survey$plot=="remo",6:ncol(survey)], binary=TRUE)

betadiv.out.fen.open  <-vegdist(survey[survey$land=="out"&survey$enclosure=="fen"&survey$plot=="open",6:ncol(survey)], binary=TRUE)
betadiv.out.fen.sarco <-vegdist(survey[survey$land=="out"&survey$enclosure=="fen"&survey$plot=="sarco",6:ncol(survey)], binary=TRUE)
betadiv.out.fen.remo  <-vegdist(survey[survey$land=="out"&survey$enclosure=="fen"&survey$plot=="remo",6:ncol(survey)], binary=TRUE)

betadiv.df<-data.frame(betadiv=c(as.numeric(betadiv.in.open),as.numeric(betadiv.in.sarco),as.numeric(betadiv.in.remo),as.numeric(betadiv.out.out.open),as.numeric(betadiv.out.out.sarco),as.numeric(betadiv.out.out.remo),as.numeric(betadiv.out.fen.open),as.numeric(betadiv.out.fen.sarco),as.numeric(betadiv.out.fen.remo)))

pairbd<-length(betadiv.in.sarco)

betadiv.df$land=factor(c(rep("in",pairbd*3),rep("out",pairbd*6)))
betadiv.df$enclosure=factor(c(rep("na",pairbd*3),rep("out",pairbd*3),rep("fen",pairbd*3)))
betadiv.df$plot=factor(rep(c(rep("open",pairbd),rep("sarco",pairbd),rep("remo",pairbd)),3))

betadiv.df$site<-rep(c(rep("papa",7),rep("papl",4),rep("papa",6),rep("papl",4),rep("papa",5),rep("papl",4),rep("papa",4),rep("papl",4),rep("papa",3),rep("papl",4),rep("papa",2),rep("papl",4),rep("papa",1),rep("papl",4),rep("papl",4),rep("plpl",3),rep("plpl",2),rep("plpl",1)),9)
betadiv.df$block<-rep(c(paste(rep(1,11),c(2,3,4,5,6,7,8,1,2,3,4)),paste(rep(2,10),c(3,4,5,6,7,8,1,2,3,4)),paste(rep(3,9),c(4,5,6,7,8,1,2,3,4)),paste(rep(4,8),c(5,6,7,8,1,2,3,4)),paste(rep(5,7),c(6,7,8,1,2,3,4)),paste(rep(6,6),c(7,8,1,2,3,4)),paste(rep(7,5),c(8,1,2,3,4)),paste(rep(8,4),c(1,2,3,4)),paste(rep(1,3),c(2,3,4)),paste(rep(2,2),c(3,4)),paste(3,4)),9)

##### stats

betadiv.df$enclosure <- relevel(betadiv.df$enclosure, ref="out")
betadiv.df$plot <- relevel(betadiv.df$plot, ref="open")

mod.betadiv <- lmer(betadiv ~ enclosure*plot + (1|site),  data= betadiv.df)

Anova(mod.betadiv)
summary(mod.betadiv)

emmeans(mod.betadiv, pairwise ~ plot)
(0.599-0.539)/0.539*100

emmeans(mod.betadiv, pairwise ~ enclosure:plot)
(0.637-0.526)/0.526*100

write.csv(summary(mod.betadiv)$coefficients, "betadiv.csv")

setEPS()
postscript("fig_ mod.betadiv.eps", widt=3.5, height=3.5)
plot(allEffects(mod.betadiv))
dev.off()

#################################################################################
##### merging with biomass
head(biomass.df)
head(survey)
survey$string <- paste(survey$site, survey$land, survey$enclosure, survey$plot, survey$blockn, sep=".")
biomass.df$alphad <-0
biomass.df$sprich <- 0

for(i in 1:nrow(biomass.df)){
	biomass.df$alphad[i] <- survey$alphadiv[which(survey$string==biomass.df$string[i])]
	biomass.df$sprich[i] <- survey$sprich[which(survey$string==biomass.df$string[i])]
}

head(biomass.df)

biomass.df$blockn[which(biomass.df$site=='plaka')] = biomass.df$blockn[which(biomass.df$site=='plaka')]+8
########
biomass.df$plot <- relevel(survey$plot, ref="open")

mod.bm <-lmer(bm ~ enclosure * plot * alphad + (1|blockn), data= biomass.df, na.action=na.exclude)

Anova(mod.bm)
summary(mod.bm)

confint(mod.bm)

emmeans(mod.bm, pairwise ~ enclosure)
(23.9-19.8)/19.8*100

emmeans(mod.bm, pairwise ~ plot)
(18.3-23.4)/18.3*100
(18.3-21.7)/18.3*100

emmeans(mod.bm, pairwise ~ plot:enclosure)

emtrends(mod.bm, ~ enclosure * plot | alphad, var= "alphad")


pdf("fig_mod.bm.all.pdf")
plot(allEffects(mod.bm))
dev.off()

setEPS()
postscript("fig_mod.bm.2.eps")
plot(Effect(c("plot","enclosure"), mod.bm))
dev.off()


plot(Effect(c("plot"), mod.bm))
plot(Effect(c("enclosure"), mod.bm))

write.csv(summary(mod.bm)$coefficients, "mod.bm.csv")
write.csv(emtrends(mod.bm, ~ plot * enclosure | alphad, var= "alphad"), "mod.bmadiv.csv")

##################################
############ figures #############

survey$plot <- relevel(survey$plot, ref="sarco")
mod.alphadiv.plot <- lmer(alphadiv ~ plot*enclosure + (1|site/blockn),  data=survey)

pdf('mod.alphadiv.pdf', 9, 5)
plot(allEffects(mod.alphadiv.plot, partial.residuals=TRUE),
     lattice = list(layout=c(3, 1)))
dev.off()

betadiv.df$plot <- relevel(betadiv.df$plot, ref="sarco")
mod.betadiv.plot <- lmer(betadiv ~ plot*enclosure + (1|site),  data= betadiv.df)

pdf('mod.betadiv.pdf', 9, 5)
plot(allEffects(mod.betadiv.plot, partial.residuals=TRUE),
     lattice = list(layout=c(3, 1)))
dev.off()


biomass.df$plot <- relevel(biomass.df$plot, ref="sarco")

pdf('mod.bm.pdf', 9, 5)
plot(Effect(c("plot", "enclosure"), mod.bm, partial.residuals=TRUE),
     lattice = list(layout=c(3, 1)))
dev.off()

pdf('mod.bmalpha.pdf', 9, 5)
ggplot(biomass.df, aes(x = alphad,
                    y = bm,
                    colour=plot)) +
  facet_wrap(vars(enclosure), ncol = 3) +
  geom_point(size = 3) +
  geom_smooth(method = lm, formula = y ~ x) +
  theme_classic()
dev.off()

