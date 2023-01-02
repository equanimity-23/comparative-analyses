# Model 1 -- wild-type vs mutant in specified virus

## DCV
DCV <- read.csv("data/DCV_930.csv")

DCV_noPBS <- filter(DCV, Virus == 1)

fitDCV <- coxme(Surv(Days, Survival) ~ Mutant + (1|BiolRep/TechRep), DCV_noPBS)
fitDCV

cox.zph(fitDCV)

DCVcurves <- survfit(Surv(Days, Survival) ~ Mutant, data=DCV_noPBS)
DCVplot <- ggsurvplot(DCVcurves,
                       data = DCV_noPBS,               
                       line = c(1,1),
                       size = 2,
                       xlab = "Days",
                       ylab = "Proportional survival",
                       xlim = c(1,18.5),
                       ylim = c(0,1.1),
                       break.time.by = 1, 
                       break.y.by = 0.2,
                       ggtheme = theme_classic(),
                       axes.offset = FALSE,
                       legend.labs = 
                         c("w1118", "miR-33-KO DCV")     #change the legends to the miRNA
)
DCVplot

sink("output/miR-33_DCV_coxme.txt") # add miRNA name and virus | creates the txt file in WD
print(fitDCV) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**


## FHV
FHV <- read.csv("data/FHV_930.csv")

FHV_noPBS <- filter(FHV, Virus == 1)

fitFHV <- coxme(Surv(Days, Survival) ~ Mutant + (1|BiolRep/TechRep), FHV_noPBS)
fitFHV

cox.zph(fitFHV)

FHVcurves <- survfit(Surv(Days, Survival) ~ Mutant, data=FHV_noPBS)
FHVplot <- ggsurvplot(FHVcurves,
                       data = FHV_noPBS,               
                       line = c(1,1),
                       size = 2,
                       xlab = "Days",
                       ylab = "Proportional survival",
                       xlim = c(1,12),
                       ylim = c(0,1.1),
                       break.time.by = 1, 
                       break.y.by = 0.2,
                       ggtheme = theme_classic(),
                       axes.offset = FALSE,
                       legend.labs = 
                         c("w1118", "miR-33-KO FHV")     #change the legends to the miRNA
)
FHVplot

sink("output/miR-33_FHV_coxme.txt") # add miRNA name and virus | creates the txt file in WD
print(fitFHV) # prints coxme output to txt file
sink() # closes the connection **CRITICAL**
