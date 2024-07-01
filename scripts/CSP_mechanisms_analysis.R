library(lme4)
library(lmerTest)
library(aod)
library(beeswarm)
library(car)
library(influence.ME)


## import data ##
#all data combined (control, 30 min, & 2 hour exposure)
datatotal <- read.csv("CSP_mechanisms_data.csv")

#all control crosses
control <- subset(datatotal, treatment == ("Control"))

#only 2 hr exposure
data2 <- subset(datatotal, time == "2hr")

#2 hr exposure divided into control (left in GIM) & treatment (left in female)
control2 <- subset(data2, treatment == ("Control"))
treatment2 <- subset(data2, treatment == ("Treatment"))

#2 hr exposure divided into heterospecific and conspecific crosses
het2 <- subset(data2, species == ("Heterospecific"))
con2 <- subset(data2, species == ("Conspecific"))

#only 30 min exposure
data30 <- subset(datatotal, time == "30min")

#testes dissections + 2 hour controls
testes2 <- read.csv("testes_2hr.csv")

#scored immediately after mating + 2 hour controls
aftermating2 <- read.csv("aftermating_2hr.csv")



## set pathway for where to save figures ##
path = "~/Desktop"


## functions to test for overdispersion & to use quasi-likelihood ##
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

quasi_table <- function(model,ctab=coef(summary(model)),
                        phi=overdisp_fun(model)["ratio"]) {
  qctab <- within(as.data.frame(ctab),
                  {   `Std. Error` <- `Std. Error`*sqrt(phi)
                  `z value` <- Estimate/`Std. Error`
                  `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                  })
  return(qctab)
}



## figure 2 ##

#code for figure 2: a) copulation duration, b) sperm transfer

femaletype = factor(datatotal$female_type, levels(as.factor(datatotal$female_type))[c(2,1)])
print(levels(femaletype))

male_species = factor(datatotal$species, levels(as.factor(datatotal$species))[c(2,1)])
print(levels(male_species))

femaletypec = factor(control$female_type, levels(as.factor(control$female_type))[c(2,1)])
print(levels(femaletypec))

male_speciesc = factor(control$species, levels(as.factor(control$species))[c(2,1)])
print(levels(male_speciesc))

labels1 <- c("", "", "","")


png(filename=paste(path,"/figure2.png",sep=""), width=14, height=7,
    units="in", res=300, pointsize=14)

par(mfrow=c(1,2), mar = c(2, 5.5, 2, 2), oma = c(3,1,0,0))
beeswarm(duration~femaletype + male_species, data = datatotal,
         cex = 1.1, pch = c(17,17,19,19), 
         labels = labels1, 
         xlab = "", ylab = "", 
         cex.axis = 1.4,
         col = c('magenta4', 'goldenrod1'))
bxplot(duration~femaletype + male_species, data = datatotal, add = T)
legend("topright", cex = 1.4, bty = "n", legend=c("Sympatric", "Allopatric"), 
       inset = 0.03 ,fill=c('magenta4', 'goldenrod1'))
mtext(text = "Male type", side = 1, line = 3.5, cex = 2)
mtext(text = "Copulation duration (min)", side = 2, line = 4, cex = 2)
mtext(text = "Heterospecific       Conspecific", side = 1, line = 1, cex = 1.6)
mtext(text = "A", side = 1, line = -23.4, cex = 2, adj = -0.15, xpd = NA)

beeswarm(sperm_total ~ femaletypec + male_speciesc, data = control,
         cex = 1.1, pch = c(17,17,19,19), 
         labels = labels1, ylim = c(10, 240),
         xlab = "", ylab = "", 
         cex.axis = 1.4, 
         col = c('magenta4', 'goldenrod1'))
bxplot(sperm_total ~ femaletypec + male_speciesc, data = control, add = T)
mtext(text = "Male type", side = 1, line = 3.5, cex = 2)
mtext(text = "Sperm count", side = 2, line = 4, cex = 2)
mtext(text = "Heterospecific       Conspecific", side = 1, line = 1, cex = 1.6)
mtext(text = "B", side = 1, line = -23.4, cex = 2, adj = -0.15, xpd = NA)

dev.off()


#statistical analysis for copulation duration (figure 2a)

femaletype = factor(datatotal$female_type, levels(as.factor(datatotal$female_type))[c(2,1)])
print(levels(femaletype))

malespecies = factor(datatotal$species, levels(as.factor(datatotal$species))[c(2,1)])
print(levels(malespecies))

duration <- lmer(duration ~ femaletype + (1 | female_line) + malespecies + femaletype:malespecies, data = datatotal)
summary(duration)


#statistical analysis for quantity of sperm transferred (figure 2b)

spermtransfer <- lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = control)
summary(spermtransfer)

spermtransfer2 <- lmer(sperm_total ~ female_type + (1 | female_line) + species, data = control)
summary(spermtransfer2)

spermcontrol <- subset(data2, treatment == ("Control"))
spermtreatment <- subset(data2, treatment == ("Treatment"))

t.test(spermcontrol$sperm_total, spermtreatment$sperm_total)




## figure 3 ##

#code for figure 3: eusperm viability, a) heterospecific, b) conspecific

femaletype2h = factor(het2$female_type, levels(as.factor(het2$female_type))[c(2,1)])
print(levels(femaletype2h))

femaletype2c = factor(con2$female_type, levels(as.factor(con2$female_type))[c(2,1)])
print(levels(femaletype2c))

labels1 <- c("", "", "","")

png(filename=paste(path,"/figure3.png",sep=""), width=14, height=7,
    units="in", res=300, pointsize=14)

par(mfrow=c(1,2), mar = c(2, 5.5, 2, 2), oma = c(3,1,0,0))
beeswarm(percent_live_eu~treatment*femaletype2h, data = het2,
         cex = 1.1, pch = c(2,17,2,17), 
         labels = labels1, 
         ylim = c(0.05, 0.95), 
         xlab = "", ylab = "", 
         cex.axis = 1.4, 
         col = c('magenta4', 'magenta4', 'goldenrod1', 'goldenrod1'))
bxplot(percent_live_eu~treatment*femaletype2h, data = het2, add = T)
legend(x = 2, y = 0.35, cex = 1.4, bty = "n", legend=c("No FRT exposure", "2 hr FRT exposure"), 
       fill=c('white', 'black'))
legend(x = 2, y = 0.2, cex = 1.4, bty = "n", legend=c("Sympatric", "Allopatric"), 
       fill=c('magenta4', 'goldenrod1'))
mtext(text = "Heterospecific matings", side = 1, line = 2, cex = 2)
mtext(text = "Proportion of live eusperm", side = 2, line = 4, cex = 2)
mtext(text = "A", side = 1, line = -23.2, cex = 2, adj = -0.15, xpd = NA)

beeswarm(percent_live_eu~treatment*femaletype2c, data = con2,
         cex = 1.1, pch = c(1,19,1,19), 
         labels = labels1, 
         ylim = c(0.05, 0.95), 
         xlab = "", ylab = "", 
         cex.axis = 1.4, 
         col = c('magenta4', 'magenta4', 'goldenrod1', 'goldenrod1'))
bxplot(percent_live_eu~treatment*femaletype2c, data = con2, add = T)
legend(x = 2, y = 0.35, cex = 1.4, bty = "n", legend=c("No FRT exposure", "2 hr FRT exposure"), 
       fill=c('white', 'black'))
legend(x = 2, y = 0.2, cex = 1.4, bty = "n", legend=c("Sympatric", "Allopatric"), 
       fill=c('magenta4', 'goldenrod1'))
mtext(text = "Conspecific matings", side = 1, line = 2, cex = 2)
mtext(text = "B", side = 1, line = -23.2, cex = 2, adj = -0.15, xpd = NA)

dev.off()


#statistical analysis for eusperm viability of 2 hour matings (figure 3)

eusperm2 <- data2[,27:28]
eusperm2 <- data.matrix(eusperm2)

eusperm <- glmer(eusperm2 ~ female_type + (1 | female_line) + species + treatment, data = data2, family = "binomial")
summary(eusperm)
overdisp_fun(eusperm)

euspermB<-betabin(formula=cbind(data2[,27],data2[,28])~female_type+species+treatment,random=~female_line,data=data2)
summary(euspermB)




## figure 4 ##

#code for figure 4: parasperm proportion & eusperm viability

para_treatment <- lm(percent_live_eu ~ percent_parasperm, data = treatment2)

png(filename=paste(path,"/figure4.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=14)

par(mfrow=c(2,1), mar = c(2, 2, 2.5, 2), oma = c(3,4,0,9))
plot(percent_live_eu ~ percent_parasperm, data = treatment2, xlab = "", ylab = "",
     cex = 1.1, cex.axis = 1.2, xlim = c(0, 0.9), ylim = c(0.43,0.95),
     col = c('magenta4', 'goldenrod1')[as.factor(treatment2$female_type)], 
     pch = c(19,17)[as.factor(treatment2$species)])
abline(para_treatment, col = "black")
legend(x = 0.95, y = 1, cex = 1.3, bty = "n", legend=c("Heterospecific", "Conspecific"), 
       pch =c(17,19), xpd = NA)
legend(x = 0.94, y = 0.85, cex = 1.3, bty = "n", legend=c("Sympatric", "Allopatric"), 
       fill=c('magenta4', 'goldenrod1'), xpd = NA)
mtext(text = "A. 2 hour FRT Exposure", side = 1, line = -12.7, cex = 1.4, adj = 0.01)

plot(percent_live_eu ~ percent_parasperm, data = control2, xlab = "", ylab = "",
     cex = 1.1, cex.axis = 1.2, xlim = c(0, 0.9), ylim = c(0.43,0.95),
     col = c('magenta4', 'goldenrod1')[as.factor(control2$female_type)], 
     pch = c(1,2)[as.factor(control2$species)])
mtext(text = "Parasperm proportion", side = 1, line = 3.5, cex = 1.8)
mtext(text = "Proportion of live eusperm", side = 2, line = 4, cex = 1.8, adj = -0.8)
mtext(text = "B. No FRT Exposure", side = 1, line = -12.7, cex = 1.4, adj = 0.01)

dev.off()


#statistical analysis for effect of parasperm proportion on eusperm viability (figure 4)

eusperm2 <- data2[,27:28]
eusperm2 <- data.matrix(eusperm2)

eupara2 <- glmer(eusperm2 ~ treatment:percent_parasperm + (1|female_line), data = data2, family = "binomial")
summary(eupara2)
overdisp_fun(eupara2)

euparaB2<-betabin(formula=cbind(data2[,27],data2[,28]) ~ treatment:percent_parasperm, random=~female_line, data = data2)
summary(euparaB2)




## figure S1 ##

#code for figure S1: testes controls

print(levels(testes2$male_line))
maleline = factor(testes2$male_line, levels(as.factor(testes2$male_line))[c(2,4,3,1,5)])
print(levels(maleline))

print(levels(testes2$control_type))
control_type = factor(testes2$control_type, levels(as.factor(testes2$control_type))[c(2,1)])
print(levels(control_type))

labels3 = c("","","","","","","","","","")

png(filename=paste(path,"/figureS1.png",sep=""), width=8, height=4,
    units="in", res=300, pointsize=7)

par(mfrow = c(1,1), mar = c(6, 6, 1, 2), oma = c(0,0,0,13.5))
beeswarm(testes2$prop_live_eu~control_type*maleline, 
         cex = 1.6, pch = c(17,2,19,1,19,1,19,1,19,1),
         labels = labels3, corral = "gutter",
         ylim = c(0.05, 0.95), 
         xlab = "", ylab = "", 
         cex.axis = 1.6, 
         col = c('magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1'))
bxplot(testes2$prop_live_eu~control_type*maleline, add = T)
legend(x = 11.1, y = 1, cex = 1.7, bty = "n", legend=c("Heterospecific", "Conspecific"), 
       pch =c(17,19), xpd = NA)
legend(x = 11, y = 0.85, cex = 1.7, bty = "n", legend=c("Testes","2 hour control"), 
       fill=c('black', 'white'), xpd = NA)
legend(x = 11, y = 0.7, cex = 1.7, bty = "n", legend=c("Sympatric", "Allopatric"), fill = c('magenta4', 'goldenrod1'), xpd = NA)
mtext(text = "Proportion of live eusperm", side = 2, line = 4, cex = 2.4)
mtext(text = "MSH         Sierra          MSH         Lamoille        Zion", side = 1, line = 1, cex = 2)
mtext(text = "Male line", side = 1, line = 4, cex = 2.4)
abline(v = 2.5, lty = 2, col = "gray")

dev.off()


#statistical analysis for testes vs. 2 hour controls (figure S1)

testes <- testes2[,10:11]
testes <- data.matrix(testes)

testesm <- glmer(testes ~ species + control_type + (1|male_line), data = testes2, family = "binomial")
summary(testesm)
overdisp_fun(testesm)

testesb<-betabin(formula=cbind(testes2[,10],testes2[,11]) ~ species + control_type, random=~male_line, data = testes2)
summary(testesb)




## figure S2 ##

#code for figure S2: after mating controls
print(levels(as.factor(aftermating2$cross_type)))
crosstype = factor(aftermating2$cross_type, levels(as.factor(aftermating2$cross_type))[c(4,3,2,1)])
print(levels(crosstype))

print(levels(aftermating2$control_type))
control_type2 = factor(aftermating2$control_type, levels(as.factor(aftermating2$control_type))[c(2,1)])
print(levels(control_type2))

labels4 = c("","","","","","","","")

png(filename=paste(path,"/figureS2.png",sep=""), width=7, height=4,
    units="in", res=300, pointsize=7)

par(mfrow = c(1,1), mar = c(6, 6, 1, 2), oma = c(0,0,0,13.5))
beeswarm(aftermating2$prop_live_eu ~ control_type2*crosstype,
         cex = 1.6, pch = c(17,2,19,1,17,2,19,1), 
         labels = labels4, 
         ylim = c(0.05, 0.95), 
         xlab = "", ylab = "", 
         cex.axis = 1.6, 
         col = c('magenta4', 'magenta4', 'magenta4', 'magenta4', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1'))
bxplot(aftermating2$prop_live_eu ~ control_type2*crosstype, add = T)
legend(x = 9.1, y = 1, cex = 1.7, bty = "n", legend=c("Heterospecific", "Conspecific"), 
       pch =c(17,19), xpd = NA)
legend(x = 9, y = 0.85, cex = 1.7, bty = "n", legend=c("0 hour control","2 hour control"), 
       fill=c('black', 'white'), xpd = NA)
legend(x = 9, y = 0.7, cex = 1.7, bty = "n", legend=c("Sympatric", "Allopatric"), fill = c('magenta4', 'goldenrod1'), xpd = NA)
mtext(text = "Proportion of live eusperm", side = 2, line = 4, cex = 2.4)
mtext(text = "MSH                                Lamoille", side = 1, line = 1, cex = 2)
mtext(text = "Female line", side = 1, line = 4, cex = 2.4)

dev.off()


#statistical analysis for immediately after mating vs. 2 hour controls (figure S2)

aftermating <- aftermating2[,16:17]
aftermating <- data.matrix(aftermating)

aftermatingm <- glmer(aftermating ~ species + control_type + (1|female_line), data = aftermating2, family = "binomial")
summary(aftermatingm)
overdisp_fun(aftermatingm)

aftermatingb<-betabin(formula=cbind(aftermating2[,16],aftermating2[,17]) ~ species + control_type, random=~female_line, data = aftermating2)
summary(aftermatingb)




## figure S3 ##

#code for figure S3: eusperm viability across lines (2 hour exposure)

print(levels(con2$female_line))
femalelineC = factor(con2$female_line, levels(as.factor(con2$female_line))[c(3,2,1,4)])
print(levels(femalelineC))

print(levels(het2$female_line))
femalelineH = factor(het2$female_line, levels(as.factor(het2$female_line))[c(3,2,1,4)])
print(levels(femalelineH))

labels1 <- c("", "", "", "", "", "", "", "")

png(filename=paste(path,"/figureS3.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=14)

par(mfrow=c(2,1), mar = c(2, 2, 2.5, 2), oma = c(3,4,0,11))
beeswarm(percent_live_eu~treatment*femalelineH, data = het2, 
         cex = 1.1, pch = c(2,17,2,17,2,17,2,17), 
         labels = labels1, ylim = c(0.15, 0.95), 
         xlab = "", ylab = "", cex.axis = 1.2, 
         col = c('magenta4', 'magenta4','magenta4', 'magenta4', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1'))
bxplot(percent_live_eu~treatment*femalelineH, data = het2, add=T)
legend(x = 9, y = 1, cex = 1.3, bty = "n", legend=c("No FRT exposure", "2 hr FRT exposure"), 
       fill=c('white', 'black'), xpd = NA)
legend(x = 9, y = 0.75, cex = 1.3, bty = "n", legend=c("Sympatric", "Allopatric"), 
       fill=c('magenta4', 'goldenrod1'), xpd = NA)
mtext(text = "A. Heterospecific", side = 1, line = -12.7, cex = 1.4, adj = -0.05)
mtext(text = "Sierra         MSH        Lamoille       Zion", side = 1, line = 0.9, cex = 1.4)

beeswarm(percent_live_eu~treatment*femalelineC, data = con2, 
         cex = 1.1, pch = c(1,19,1,19,1,19,1,19), 
         labels = labels1, ylim = c(0.15, 0.95), 
         xlab = "", ylab = "", cex.axis = 1.2, 
         col = c('magenta4', 'magenta4', 'magenta4', 'magenta4', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1'))
bxplot(percent_live_eu~treatment*femalelineC, data = con2, add=T)
mtext(text = "B. Conspecific", side = 1, line = -12.7, cex = 1.4, adj = -0.05)
mtext(text = "Sierra         MSH        Lamoille       Zion", side = 1, line = 0.9, cex = 1.4)
mtext(text = "Female line", side = 1, line = 3, cex = 1.8)
mtext(text = "Proportion of live eusperm", adj = -0.8, side = 2, line = 4, cex = 1.8)

dev.off()


#statistical analysis for eusperm viability: variation among female lines (figure S3)

eusperm2 <- data2[,27:28]
eusperm2 <- data.matrix(eusperm2)

eu <- glm(eusperm2 ~0 + female_line + species + treatment, data = data2, family = "binomial")

overdisp_fun(eu)
quasi_table(eu)

wald.test(b=coef(eu), Sigma=vcov(eu)*2.4935, Terms=1:4)




## figure S4 ##

#code for figure S4: all sperm viability (30 minute exposure)

print(levels(as.factor(data30$female_line)))
femaleline = factor(data30$female_line, levels(as.factor(data30$female_line))[c(6,5,3,4,1,2,7,8)])
print(levels(femaleline))

labels2 = c("","","","","","","","","","","","","","","","")

png(filename=paste(path,"/figureS4.png",sep=""), width=8, height=4,
    units="in", res=300, pointsize=7)

par(mfrow = c(1,1), mar = c(8, 6, 1, 2), oma = c(0,0,0,15))
beeswarm(percent_live ~ treatment*femaleline, data = data30,
         cex = 1.6, pch = c(2,17,2,17,2,17,2,17,2,17,2,17,2,17,2,17), 
         labels = labels2, corral = "gutter",
         ylim = c(0.05, 0.95), 
         xlab = "", ylab = "", 
         cex.axis = 1.6, 
         col = c('magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1'))
bxplot(percent_live ~ treatment*femaleline, data = data30, add = T)
legend(x = 17, y = 1, cex = 1.7, bty = "n", legend=c("No FRT exposure", "30 min FRT exposure"), 
       fill=c('white', 'black'), xpd = NA)
legend(x = 17, y = 0.85, cex = 1.7, bty = "n", legend=c("Sympatric", "Allopatric"),
       fill=c('magenta4', 'goldenrod1'), xpd = NA)
mtext(text = "Proportion of live sperm", side = 2, line = 4, cex = 2.4)
mtext(text = "Line 1    Line 2    Line 3    Line 4    Line 5    Line 6    Line 7    Line 8", side = 1, line = 1, cex = 1.7)
mtext(text = "Sierra               MSH             Lamoille             Zion", side = 1, line = 3.5, cex = 2)
mtext(text = "Female line", side = 1, line = 6, cex = 2.4)
abline(v = c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5), lty = 2, col = "gray")

dev.off()


#statistical analysis for total sperm viability: 30 minutes (figure S4)

viability30m <- glmer(formula=cbind(data30[,21],data30[,20]) ~ female_type + (1|female_line), data = data30, family = "binomial")
summary(viability30m)
overdisp_fun(viability30m)

viabilityb<-betabin(formula=cbind(data30[,21],data30[,20]) ~ female_type, random=~female_line, data = data30)
summary(viabilityb)




## figure S5 ##

#code for figure S5: parasperm proportion (spermblue data)

print(levels(as.factor(data30$female_line)))
femaleline = factor(data30$female_line, levels(as.factor(data30$female_line))[c(6,5,3,4,1,2,7,8)])
print(levels(femaleline))

labels2 = c("","","","","","","","","","","","","","","","")

png(filename=paste(path,"/figureS5.png",sep=""), width=8, height=4,
    units="in", res=300, pointsize=7)

par(mfrow = c(1,1), mar = c(8, 6, 1, 2), oma = c(0,0,0,15))
beeswarm(percent_parasperm ~ treatment*femaleline, data = data30,
         cex = 1.6, pch = c(2,17,2,17,2,17,2,17,2,17,2,17,2,17,2,17), 
         labels = labels2, corral = "gutter",
         ylim = c(0.05, 0.95), 
         xlab = "", ylab = "", 
         cex.axis = 1.6, 
         col = c('magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'magenta4', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1', 'goldenrod1'))
bxplot(percent_parasperm ~ treatment*femaleline, data = data30, add = T)
legend(x = 17, y = 1, cex = 1.7, bty = "n", legend=c("No FRT exposure", "30 min FRT exposure"), 
       fill=c('white', 'black'), xpd = NA)
legend(x = 17, y = 0.85, cex = 1.7, bty = "n", legend=c("Sympatric", "Allopatric"),
       fill=c('magenta4', 'goldenrod1'), xpd = NA)
mtext(text = "Proportion of parasperm", side = 2, line = 4, cex = 2.4)
mtext(text = "Line 1    Line 2    Line 3    Line 4    Line 5    Line 6    Line 7    Line 8", side = 1, line = 1, cex = 1.7)
mtext(text = "Sierra               MSH             Lamoille             Zion", side = 1, line = 3.5, cex = 2)
mtext(text = "Female line", side = 1, line = 6, cex = 2.4)
abline(v = c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5), lty = 2, col = "gray")

dev.off()


#statistical analysis for parasperm proportion across female lines (figure S5)

parasperm30 <- data30[,23:24]
parasperm30 <- data.matrix(parasperm30)

para <- glmer(parasperm30 ~ treatment + (1|female_line), data = data30, family = "binomial")
summary(para)
overdisp_fun(para)

parab <- betabin(formula=cbind(data30[,23],data30[,24]) ~ treatment, random=~female_line, data = data30)
summary(parab)

parabfemale <- betabin(formula=cbind(data30[,23],data30[,24]) ~ female_type, random=~female_line, data = data30)
summary(parabfemale)




## figure S6 ##

#code for figure S6: eusperm viability & parasperm proportion

euspermt <- treatment2[,27:28]
euspermt <- data.matrix(euspermt)

euspermc <- control2[,27:28]
euspermc <- data.matrix(euspermc)

eupara_treatment <- glmer(euspermt ~ female_type + (1 | female_line) + species + percent_parasperm, data = treatment2, family = "binomial")
eupara_control <- glmer(euspermc ~ female_type + (1 | female_line) + species + percent_parasperm, data = control2, family = "binomial")

png(filename=paste(path,"/figureS6.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=14)

par(mfrow=c(2,1), mar = c(2, 2, 2.5, 2), oma = c(3,4,0,9))
plot(x = treatment2$percent_parasperm, y = fitted(eupara_treatment), 
     xlab = "", ylab = "",
     xlim = c(0.05,0.85), ylim = c(0.6,0.9), 
     cex = 1.1, cex.axis = 1.2, 
     col = c('goldenrod1','magenta4')[as.factor(treatment2$female_type)], 
     pch = c(19,17)[as.factor(treatment2$species)])
mtext(text = "A. 2 hour FRT Exposure", side = 1, line = -12.7, cex = 1.4, adj = 0.01)
legend(x = 0.91, y = 0.93, cex = 1.3, bty = "n", legend=c("Heterospecific", "Conspecific"), 
       pch =c(17,19), xpd = NA)
legend(x = 0.9, y = 0.85, cex = 1.3, bty = "n", legend=c("Sympatric", "Allopatric"), 
       fill=c('magenta4', 'goldenrod1'), xpd = NA)

plot(x = control2$percent_parasperm, y = fitted(eupara_control), 
     xlab = "", ylab = "",
     xlim = c(0.05,0.85), ylim = c(0.6,0.9), 
     cex = 1.1, cex.axis = 1.2,
     col = c('goldenrod1','magenta4')[as.factor(treatment2$female_type)], 
     pch = c(1,2)[as.factor(treatment2$species)])
mtext(text = "Proportion of live eusperm (model fitted values)", adj = 0.12, side = 2, line = 4, cex = 1.8)
mtext(text = "Parasperm proportion", side = 1, line = 3.5, cex = 1.8)
mtext(text = "B. No FRT Exposure", side = 1, line = -12.7, cex = 1.4, adj = 0.01)

dev.off()


#statistical analysis for parasperm proportion vs. eusperm viability (figure S6; same analysis as figure 3)

eupara2 <- glmer(eusperm2 ~ treatment:percent_parasperm + (1|female_line), data = data2, family = "binomial")
summary(eupara2)
overdisp_fun(eupara2)

euparaB2<-betabin(formula=cbind(data2[,27],data2[,28]) ~ treatment:percent_parasperm, random=~female_line, data = data2)
summary(euparaB2)




## influential outlier analysis ##

#figure 2a: copulation duration

dataf1a<-data.frame(datatotal$duration,datatotal$female_line,femaletype,malespecies)
colnames(dataf1a)<-c("duration", "female_line", "femaletype", "malespecies")
durationsub<-lmer(duration ~ femaletype + (1 | female_line) + malespecies + femaletype:malespecies, data = dataf1a)

#effect of female line outliers

fl<-influence(durationsub,"female_line")

plot(fl, which="dfbetas", parameters=c(2), xlab="DFbetaS",ylab="Female Line")
2/sqrt(8)
dfbetas(fl, parameters=c(2))

plot(fl, which="cook", cutoff=4/8, sort=TRUE, xlab="Cook´s Distance", ylab="Female Line")
sigtest(fl,test=1.96,)

dsMSH68<-exclude.influence(durationsub,"female_line", "MSH68")

#effect of individual observation level outliers

flo<-influence(durationsub, obs=TRUE)

plot(flo,which="cook",cutoff=4/160,sort=TRUE, xlab="Cook's Distance", ylab="Individual Observation")
cooks.distance(flo,sort=T)

summary(durationsub)
sigtest(flo,test=1.96)$femaletypeAllopatric[c(79,158,29,50),]
sigtest(flo,test=1.96)$malespeciesConspecific[c(79,158,29,50),]

dataf1a79<-dataf1a[-79,]
dataf1a158<-dataf1a[-158,]
dataf1a29<-dataf1a[-29,]
dataf1a50<-dataf1a[-50,]
summary(durationsub)
summary(lmer(duration ~ femaletype + (1 | female_line) + malespecies + femaletype:malespecies, data = dataf1a79))
summary(lmer(duration ~ femaletype + (1 | female_line) + malespecies + femaletype:malespecies, data = dataf1a158))
summary(lmer(duration ~ femaletype + (1 | female_line) + malespecies + femaletype:malespecies, data = dataf1a29))
summary(lmer(duration ~ femaletype + (1 | female_line) + malespecies + femaletype:malespecies, data = dataf1a50))


#figure 2b: sperm transfer

dataf1b<-data.frame(control$sperm_total,control$female_line,control$female_type,control$species,control$duration)
colnames(dataf1b)<-c("sperm_total", "female_line", "female_type", "species", "duration")
spermtotalsub<-lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b)

#effect of female line outliers

fl2<-influence(spermtotalsub,"female_line")

plot(fl2, which="dfbetas", parameters=c(2), xlab="DFbetaS",ylab="Female Line")
2/sqrt(8)
dfbetas(fl2, parameters=c(2))

plot(fl2, which="cook", cutoff=4/8, sort=TRUE, xlab="Cook´s Distance", ylab="Female Line")
sigtest(fl2,test=1.96,)

dsMVMS1<-exclude.influence(spermtotalsub,"female_line", "MVMS1-1")

#effect of individual observation level outliers

flo2<-influence(spermtotalsub, obs=TRUE)

plot(flo2,which="cook",cutoff=4/160,sort=TRUE, xlab="Cook's Distance", ylab="Individual Observation")
cooks.distance(flo2,sort=T)

summary(spermtotalsub)
sigtest(flo2,test=1.96)$female_typeSympatric[c(72,26,25,56,19,49,39,6),]
sigtest(flo2,test=1.96)$speciesHeterospecific[c(72,26,25,56,19,49,39,6),]
sigtest(flo2,test=1.96)$duration[c(72,26,25,56,19,49,39,6),]

dataf1b72<-dataf1b[-72,]
dataf1b26<-dataf1b[-26,]
dataf1b25<-dataf1b[-25,]
dataf1b56<-dataf1b[-56,]
dataf1b19<-dataf1b[-19,]
dataf1b49<-dataf1b[-49,]
dataf1b39<-dataf1b[-39,]
dataf1b6<-dataf1b[-6,]
summary(spermtotalsub)
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b72))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b26))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b25))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b56))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b19))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b49))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b39))
summary(lmer(sperm_total ~ female_type + (1 | female_line) + species + duration, data = dataf1b6))


#figure 3: eusperm viability

eusub<-data.frame(data2[,27:28],as.factor(data2$female_type),data2$female_line,data2$species,data2$treatment)
colnames(eusub)<-c("Live","Dead","female_type","female_line","species","treatment")

#effect of female line

feu<-influence(meusub,"female_line")
plot(feu, which="cook", cutoff=4/4, sort=TRUE, xlab="Cook´s Distance", ylab="Female Line")
meusub<-glmer(cbind(Live,Dead)~ female_type + species + treatment +(1 | female_line) , data = eusub, family = "binomial")

#effect of individual observation level outliers

feuo<-influence(meusub,obs=TRUE)
plot(feuo, which="cook", cutoff=4/64, sort=TRUE, xlab="Cook´s Distance", ylab="Observation")
cooks.distance(feuo,sort=T)
eusub[c(50,22,6,61,26,33,20,24,48,29,51,42,38,1),]

eusub2<-eusub[-c(50,22,6,61,26,33,20,24,48,29,51,42,38,1),]
meusub2<-glmer(cbind(Live,Dead)~ female_type + species + treatment +(1 | female_line) , data = eusub2, family = "binomial")
summary(meusub2)


#figure 4: eusperm viability and parasperm proportion

euparasub<-data.frame(data2[,27:28],data2$female_line,data2$percent_parasperm,data2$treatment)
colnames(euparasub)<-c("Live","Dead","female_line","percent_parasperm","treatment")
meuparasub<-glmer(cbind(Live,Dead)~ treatment:percent_parasperm +(1 | female_line) , data = euparasub, family = "binomial")

#effect of female line

feu2<-influence(meuparasub,"female_line")
plot(feu2, which="cook", cutoff=4/4, sort=TRUE, xlab="Cook´s Distance", ylab="Female Line")

#effect of individual observation level outliers

feuo2<-influence(meuparasub,obs=TRUE)
plot(feuo2, which="cook", cutoff=4/64, sort=TRUE, xlab="Cook´s Distance", ylab="Observation")
cooks.distance(feuo2,sort=T)
euparasub[c(1,38,20,42,51,2,48,36,26,54,33,4),]

euparasub2<-euparasub[-c(1,38,20,42,51,2,48,36,26,54,33,4),]
meuparasub2<-glmer(cbind(Live,Dead)~ treatment:percent_parasperm +(1 | female_line) , data = euparasub2, family = "binomial")
summary(meuparasub2)



