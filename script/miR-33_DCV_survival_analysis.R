#Installs the dplyr program for data cleanup
install.packages("dplyr")
library(dplyr)

#Installs the Cox Mixed Effect Model and downloads it into the function library 
install.packages("coxme")
library(coxme)

#Installs the survminer program for ggplot survival curves plotting.
install.packages("survminer")
library(survminer)

#Renames the data files, input file name here 

DCV <- read.csv("930_DCV.csv")

#Filters out the PBS samples to stop them from confounding the analysis

DCV_noPBS <- filter(DCV, Virus == 1)


#Fits the coxme to the survival data. 
#Where the fixed effects Virus and the random effects 
#are those of the technical and biological replications
#Tests if the survival curves of miRNA KO mutants are 
#significantly different  to the wild-type. 
#Only uses the data of virus-infected flies. 

fitDCV <- coxme(Surv(Days, Survival)~ Mutant + (1|BiolRep/TechRep), 
                data = DCV_noPBS)
fitDCV

#test for proportional hazards assumption
fitCoxphDCV <- coxph(Surv(Days, Survival)~ Mutant+Virus, data = DCV_noPBS)
cox.zph(fitDCV)


#Effect of mutant and virus infection on mortality 
#(simpler model, no interactions between mutant and virus)
##fitmiRNA_Virus1 <- coxme(Surv(Days, Survival)~ Mutant+Virus + 
##                           (1|BiolRep/TechRep), data=FHV)
##fitmiRNA_Virus1

#complex model
##fitmiRNA_Virus2 <- coxme(Surv(Days, Survival)~ Mutant*Virus + 
##                          (1|BiolRep/TechRep), data=FHV)
##fitmiRNA_Virus2



#Fits survival curves based on the survival data
#Wild-type vs mutants
miRNAcurves <- survfit(Surv(Days, Survival) ~ Mutant, data=DCV_noPBS)
survplot <- ggsurvplot(miRNAcurves,
           data = DCV_noPBS,               #change to the virus under analysis
           line = c(1,1),
           size = 5,
           xlab = "Days",
           ylab = "Proportional survival",
           font.x = c(60),
           font.y = c(60),
           font.tickslab = c(40, "bold"),
           xlim = c(1,18.5),
           ylim = c(0,1.1),
           break.time.by = 1, 
           break.y.by = 0.2,
           ggtheme = theme_classic(),
           axes.offset = FALSE,
           legend.labs = 
             c("w1118", "miR-33-KO DCV")     #change the legends to the miRNA
)

survplot$plot <- survplot$plot + 
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1, 
    label = "HR = 0.08 \n p < 0.001",
    size = 30) +
  theme(axis.ticks.length = unit(0.5, "cm"))

survplot

#save the plot as miRNA_Virus in the folder of your choice
#use the export function to save as pdf
#pdf -> device size -> preview
#preview before saving to make sure it is correct.


#The code below is for legacy purposes only, use the ggsurvplot code above
#Effect of virus
#Title for the name of the graphs Graph_Survival_Sus_FHv_108_v01
curvesVirus <- survfit(Surv(Days, Survival) ~ Mutant, data=FHV_noPBS)
plot(curvesVirus, main = "Virus", xlab="Days postinfection", ylab="Proportion Alive", xlim=c(0,10), ylim=c(0.1,1), col=c('gray0','gray70'), lty=c(1,2))
legend( 5,0.5 , c("p=0.013"),bty ="n")



#Effect of Mutant and VIrus
#Save File Name Graph_Wol_V_108_v01
curvesmiRNA_Virus <- survfit(Surv(Days, Survival) ~ Virus + Mutant, data=FHV)
plot(curvesmiRNA_Virus, main="DCV", xlab="Days postinfection", ylab="Proportion Alive", xlim=c(0,10), ylim=c(0.1,1), col=c('gray0','gray70','gray0','gray70'), lty=c(1,1,2,2))
legend( 4,0.5, c("pl-Tet PBS","pl-wPanCI PBS", "pl-Tet DCV", "pl-wPanCI DCV"), col=c('gray0','gray70','gray0','gray70'), lty=c(1,1,2,2))
