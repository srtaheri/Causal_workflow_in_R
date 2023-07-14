library(dagitty)
library(ggdag)
library(dosearch)
source("/Users/sarataheri/GitHub/OptimalAdjustmentSet/main_functions.R")

EColi_dagitty <- dagitty( "dag {
  chiX [latent]
  citX [latent]
  dsrA [latent]
  gcvB [latent]
  srIR [latent]
  appX [latent]
  yjjQ [latent]
  citD [latent]
  appY -> appA
  appY -> appB
  appY -> appX
  appY -> hyaA
  appY -> hyaB
  appY -> hyaF
  arcA -> aceE
  arcA -> appY
  arcA -> citX
  arcA -> cydD
  arcA -> dpiA
  arcA -> dpiB
  arcA -> fnr
  arcA -> gcvB
  arcA -> hyaA
  arcA -> hyaB
  arcA -> hyaF
  arcA -> mdh
  arcA -> rpoS
  btsR -> mdh
  chiX -> dpiA
  chiX -> dpiB
  citX -> dpiB
  cra -> cyoA
  crp -> aceE
  crp -> cirA
  crp -> citX
  crp -> cyoA
  crp -> dcuR
  crp -> dpiA
  crp -> dpiB
  crp -> exuT
  crp -> fis
  crp -> fur
  crp -> gadX
  crp -> mdh
  crp -> oxyR
  crp -> srIR
  cspA -> hns
  dcuR -> dpiA
  dcuR -> dpiB
  dpiA -> appY
  dpiA -> citC
  dpiA -> citD
  dpiA -> citX
  dpiA -> dpiB
  dpiA -> exuT
  dpiA -> mdh
  dsrA -> hns
  dsrA -> lrp
  dsrA -> rpoS
  fis -> cyoA
  fis -> gadX
  fis -> hns
  fis -> hyaA
  fis -> hyaB
  fis -> hyaF
  fnr -> aceE
  fnr -> amtB
  fnr -> aspC
  fnr -> citX
  fnr -> cydD
  fnr -> cyoA
  fnr -> dcuR
  fnr -> dpiA
  fnr -> dpiB
  fnr -> gadX
  fnr -> hcp
  fnr -> narL
  fur -> amtB
  fur -> aspC
  fur -> cirA
  fur -> cyoA
  fur -> fnr
  gadX -> amtB
  gadX -> hns
  gcvB -> lrp
  gcvB -> oxyR
  gcvB -> ydeO
  hns -> appY
  hns -> srIR
  hns -> ydeO
  hns -> yjjQ
  ihfA -> crp
  ihfA -> fnr
  ihfA -> ihfB
  ihfB -> fnr
  iscR -> hyaA
  iscR -> hyaB
  iscR -> hyaF
  lrp -> aspC
  lrp -> soxS
  modE -> narL
  narL -> citX
  narL -> cydD
  narL -> dcuR
  narL -> dpiA
  narL -> dpiB
  narL -> hcp
  narL -> hyaA
  narL -> hyaB
  narL -> hyaF
  narP -> hyaA
  narP -> hyaB
  narP -> hyaF
  oxyR -> fur
  oxyR -> hcp
  phoB -> cra
  rpoD -> aceE
  rpoD -> appY
  rpoD -> arcA
  rpoD -> cirA
  rpoD -> crp
  rpoD -> cydD
  rpoD -> dcuR
  rpoD -> dsrA
  rpoD -> fis
  rpoD -> fnr
  rpoD -> fur
  rpoD -> gcvB
  rpoD -> hns
  rpoD -> hyaA
  rpoD -> hyaB
  rpoD -> hyaF
  rpoD -> ihfB
  rpoD -> mdh
  rpoD -> narL
  rpoD -> oxyR
  rpoD -> phoB
  rpoD -> soxS
  rpoD -> srIR
  rpoD -> ydeO
  rpoD -> yjjQ
  rpoH -> cra
  rpoS -> aceE
  rpoS -> appY
  rpoS -> hyaA
  rpoS -> hyaB
  rpoS -> hyaF
  rpoS -> ihfA
  rpoS -> ihfB
  rpoS -> oxyR
  soxS -> fur
  srIR -> gutM
  ydeO -> hyaA
  ydeO -> hyaB
  ydeO -> hyaF  
}")

## read the data
gene_blattner = read.csv("/Users/sarataheri/GitHub/OptimalAdjustmentSet/EColi/real_data/gene-blattner.tab.csv", sep = "")
## read the data
whole_data <- read.csv("/Users/sarataheri/GitHub/OptimalAdjustmentSet/EColi/real_data/Expression_data.csv")
#whole_data <- read.csv("data/Ecoli_ternary_data.csv")
#Extract the related rows:
gene_blattner_select = gene_blattner[which(gene_blattner$Name %in% names(EColi_dagitty)),]
gene_blattner_select = gene_blattner_select[order(gene_blattner_select$Blattner),]
selected_rows = which(whole_data$log.TPM %in% gene_blattner_select$Blattner)

data <- whole_data[selected_rows,]
#names are covariates plus exposure and outcome
n_df = data.frame("names" = gene_blattner_select[which(gene_blattner_select$Blattner %in% data$log.TPM),]$Name)
data <- cbind(n_df,data)
#Create observational data
remove = c("fur__delfur_dpd__1",
           "fur__delfur_fe2__1",
           "fur__delfur_fe2__2",
           "acid__delgade_ph5__1",
           "acid__delgade_ph5__2",
           "oxidative__delsoxr_pq__1",
           "oxidative__delsoxr_pq__2",
           "oxidative__delsoxs_pq__1",
           "oxidative__delsoxs_pq__2",
           "ompr__bw_delompr_nacl__1",
           "ompr__bw_delompr_nacl__2",
           "crp__delcrp_fru__1",
           "crp__delcrp_fru__2",
           "crp__delcrp_fru__3",
           "crp__delcrp_glc__1",
           "crp__delcrp_glc__2",
           "crp__delcrp_glc__3",
           "crp__delcrp_glyc__1",
           "crp__delcrp_glyc__2",
           "crp__delcrp_glyc__3"
)
obs_data = data[ , !(names(data) %in% remove)]
obs_data = t(obs_data[,3:ncol(obs_data)])
colnames(obs_data) = data$names
rownames(obs_data) = seq(1:nrow(obs_data))
obs_data = as.data.frame(obs_data)
#write.csv(obs_data, "/Users/sarataheri/GitHub/OptimalAdjustmentSet/EColi/real_data/obs_data_EColi_medium.csv")
#test conditional independencies, add a bi-directed edge for each failed test
EColi_dagitty <- add_bi_directed_edges_for_failed_cond_indpd_tests(g = EColi_dagitty, data = obs_data, type = "cis")

# Find all the mediators (variables on the proper causal paths) between exposure (fur) and outcome (dpiA)
nodes_on_proper_causal_paths <- get_nodes_on_proper_paths(g = EColi_dagitty, exposure = "fur", outcome = "dpiA", current_path = character(0), all_nodes = character(0))
mediators <- nodes_on_proper_causal_paths[!nodes_on_proper_causal_paths == "fur"]
# (1)Find all the descendants of the mediators excluding the mediators
descendants_of_mediators <- get_descendants_of_mediators(g = EColi_dagitty, nodes = mediators)

#update the network and assign all the descendants of mediators to as well as all the variables that are not exposure, outcome, or part of the canonical adjustment set to latent
EColi_dagitty_some_nodes_set_as_latent <- assign_observable_nodes_as_latent(g = EColi_dagitty, nodes_to_assign_as_latent = descendants_of_mediators)

# Apply Evan's simplification rules
EColi_dagitty_some_nodes_set_as_latent <- convert_bi_dir_edges_to_latent_nodes(EColi_dagitty_some_nodes_set_as_latent)
EColi_simplified <- generate_simplified_graph(g = EColi_dagitty_some_nodes_set_as_latent)

#test conditional independencies again, add a bi-directed edge for each failed test
EColi_final <- add_bi_directed_edges_for_failed_cond_indpd_tests(g = EColi_simplified, data = obs_data, type = "cis")

canonical_adj1 <- adjustmentSets(EColi_final, "fur", "dpiA", type = "canonical") #Find the canonical adjustment set


#Estimate E(dpiA | do(fur=0)) with canonical adjustment set
lm <- lm(dpiA ~ fur + arcA + crp + ihfA + ihfB + lrp + modE + oxyR + rpoD + rpoS + soxS, obs_data)
do_calc_data <- obs_data
do_calc_data["fur"] = 0
summary(predict(lm, do_calc_data))

#Estimate E(dpiA | do(fur=0)) with canonical adjustment set
lm2 <- lm(dpiA ~ fur + arcA + crp + ihfA  + lrp + oxyR + rpoD + rpoS + soxS, obs_data)
do_calc_data <- obs_data
do_calc_data["fur"] = 0
summary(predict(lm2, do_calc_data))

#validate the results
#Create interventional data
intv_data <- as.numeric(data[which(data$names == "dpiA"), c("fur__delfur_dpd__1", "fur__delfur_fe2__1", "fur__delfur_fe2__2")])
mean(intv_data)

#further simplification
#(2) Find all the variables that are not exposure, outcome, or part of the canonical adjustment set
canonical_adj1 <- adjustmentSets(EColi_dagitty, "fur", "dpiA", type = "canonical") #Find the canonical adjustment set
#nodes that are not in canonical adjustment set
nodes_except_exposure_outcome <- names(EColi_dagitty)[!(names(EColi_dagitty) %in% c("fur","dpiA"))]
nodes_assign_to_latent <- nodes_except_exposure_outcome[!(nodes_except_exposure_outcome %in% canonical_adj1[[1]])]