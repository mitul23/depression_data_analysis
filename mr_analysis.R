library(dplyr)
library(tidyr)
library(devtools)
library(TwoSampleMR)
library(metafor)
library(meta)
library(data.table)

# depression data only with significant snps with p <= 5e-8
MDD_sig <- fread("data/PGC_UKB_depression_genome-wide.fix.LDclumped.tsv", sep="\t",header=T)%>% 
  as_tibble() %>% rename(beta=LogOR, se = StdErrLogOR,
                         effect_allele = A1, other_allele = A2,
                         eaf = Freq, pval = pval.exposure)

# Harmonised data from the MRbase site (http://www.mrbase.org/)- choose exposure = MDD_sig, outcome_exposure = c(1025, 1026)
harmonised <- read.csv("data/mr_base_harmonised_data.csv") %>% as_tibble()
unique(harmonised$SNP) # removed duplicated SNPs, keep the one in pair from 1025
harmonised %>% filter(harmonised$outcome == "Multiple sclerosis || id:1024") %>% dim() # 20 snps from 1024, 7 from 1025

# perform MR using pacakge: [MRBase](http://www.mrbase.org/)
mr_results <- mr(harmonised)
write.csv(mr_results, file="data/mr_results.csv")

# Meta-analysis using meta package: for continuous outcome (https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/fixed.html)
# m <- metagen(beta.outcome,
#             se.outcome,
#             data=harmonised,
#             studlab=paste(outcome.deprecated),
#             comb.fixed = TRUE,
#             comb.random = FALSE,
#             prediction=TRUE,
#             sm="SMD")
# save 
#sink("results.txt")
#print(m)
#sink()
# forest plot of the above result
#forest(m)
# save the plot in a particular format -Journal of the American Medical Association
#pdf(file = 'forestplot.pdf') 
#forest.jama <- forest(m,
#                      layout = "JAMA",
#                      text.predict = "95% PI",
#                      col.predict = "black")
#dev.off() 

# alternate package - metafor
#m.fe <- rma(beta.outcome, (se.outcome)**2, data=harmonised, method = "FE") # requires - effect size, variance

# Meta-analysis for binary outcome
# exponential of the beta.outcome
harmoised_altered <- harmonised %>% mutate(exp.beta.outcome = exp(beta.outcome))

# consider IVW and MR egger cols of the mr_results to perform Meta-analysis. 
mr_results <- mr_results %>% mutate(exp.b = exp(b)) %>% rename(Study = outcome)
IVW <- mr_results[c(3,8),]
IVW_MRegger <- mr_results[c(3,8,1),]

m_IVW <- metagen(exp.b, 
                 se, 
                 studlab = id.outcome,
                 method.tau = "SJ",
                 sm = "exp.b",
                 pval = pval,
                 comb.fixed = T,
                 comb.random = F,
                 data = IVW) 
m_IVW
sink("m_IVW.txt")
print(m_IVW)
sink()
forest(m_IVW,
       studlab = T)
pdf("data/m_IVW_forestplot.pdf")
m_IVW_MRegger <- metagen(exp.b, 
                         se, 
                         studlab = id.outcome,
                         method.tau = "SJ",
                         sm = "exp.b",
                         pval = pval,
                         comb.fixed = T,
                         comb.random = F,
                         data = IVW_MRegger) 
sink("m_IVW_MRegger.txt")
print(m_IVW_MRegger)
sink()
forest(m_IVW_MRegger,
       studlab = T)
pdf("data/m_IVW_MRegger_forestplot.pdf")

m_harmoised_alt <- metagen(exp(beta.exposure), 
                           se.exposure, 
                           studlab = SNP,
                           method.tau = "SJ",
                           sm = "beta.exposure",
                           pval = pval.exposure,
                           comb.fixed = T,
                           comb.random = F,
                           data = harmoised_altered) 
forest(m_harmoised_alt, studlab = T)
