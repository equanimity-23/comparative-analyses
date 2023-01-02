# Model 2 -- wild-type vs mutant, but with differences in the virus (comparative)
# studying the interactions between the miRNA knockout and virus if there is a difference between DCV and FHV

dat <- read.csv("data/Compare_930.csv")

dat$Genotype <- relevel(dat$Genotype, ref = "WT")

dat$Genotype <- factor(dat$Genotype, levels = c("WT", "m33-KO"))

fitCompare <- coxme(Surv(Days, Survival) ~ Genotype + Virus + (1|BiolRep/TechRep), dat)
fitCompare

fitInteraction <- coxme(Surv(Days, Survival) ~ Genotype + Virus + Genotype:Virus + (1|BiolRep/TechRep), dat)
fitInteraction

fitAlt <- coxph(Surv(Days, Survival) ~ Mutant + Virus + Mutant:Virus, dat)
fitAlt
