library(bnlearn)
library(dagitty)
library(dplyr)
library(dosearch) #can check if query is identifiable from any combination of observational and interventional distributions
source("/Users/sarataheri/GitHub/OptimalAdjustmentSet/main_functions.R")

# observational data
# The data are continuous, as they represent the concentration of the molecules
# under investigation.
sachs <- read.table("Sachs_case_study/sachs.data.txt", header = TRUE)

# discretizing the observational data
dsachs <- discretize(sachs, method = "hartemink", breaks = 3, ibreaks = 60, idisc = "quantile")

# interventional data
isachs <- read.table("Sachs_case_study/sachs.interventional.txt", header = TRUE, colClasses = "factor")

isachs_noIntv = isachs[which(isachs$INT == 0),] #1800 rows
isachs_MekIntv = isachs[which(isachs$INT == 2),] #600 rows
isachs_PIP2Intv = isachs[which(isachs$INT == 4),] #600 rows
isachs_AktIntv = isachs[which(isachs$INT == 7),] #600 rows
isachs_PKAIntv = isachs[which(isachs$INT == 8),] #600 rows
isachs_PKCIntv = isachs[which(isachs$INT == 9),] #1200 rows

# Prior knowledge on the network
sachs.prior.dag <- dagitty('dag {
                              PKA [latent]
                              Plcg -> PKC
                              Plcg -> PIP2
                              Plcg -> PIP3
                              PIP3 -> PIP2
                              PIP3 -> Akt
                              PIP2 -> PKC
                              PKC -> Mek
                              PKC -> Raf
                              PKC -> PKA
                              PKC -> Jnk 
                              PKC -> P38
                              PKA -> Raf
                              PKA -> Mek
                              PKA -> Erk
                              PKA -> Akt
                              PKA -> Jnk
                              PKA -> P38
                              Raf -> Mek -> Erk -> Akt
                            }')

sachs.dag_validated_against_obs_data <- add_bi_directed_edges_for_failed_cond_indpd_tests(g = sachs.prior.dag, data = dsachs, type = "cis.chisq")
sachs.dag_validated <- add_bi_directed_edges_for_failed_cond_indpd_tests(g = sachs.dag_validated_against_obs_data, data = isachs_PKCIntv, type = "cis.chisq")

#assign all the variables that are not part of the estimand to latent
sachs_variables_not_in_estimant_as_latent <- assign_observable_nodes_as_latent(g = sachs.dag_validated, nodes_to_assign_as_latent = c("Akt", "Jnk", "P38"))
# Apply Evan's simplification rules
converted_g <- convert_bi_dir_edges_to_latent_nodes(sachs_variables_not_in_estimant_as_latent)
sachs_simplified <- generate_simplified_graph(g = converted_g)

# P(Erk|do(Raf, Mek)) is not identifiable from only observational data when PKA is latent. 
# Total run time: 12 mins
do_search_sachs_obs_data = "P(Plcg,PIP3,PIP2,PKC,Raf,Mek,Erk,Akt,Jnk,P38)"
start = Sys.time()
dosearch(data = do_search_sachs_obs_data,
         query = "P(Erk|do(Raf, Mek))",
         graph = sachs.dag_validated)
end = Sys.time()
end - start

# P(Erk|do(Raf, Mek)) is identifiable from both observational and interventional data
# only on Mek when PKA is latent. It is equal to P(Akt|do(Mek))
# Total run time: 20 seconds
# do_search_sachs_obs_intv_data <- "P(Plcg,PIP3,PIP2,PKC,Raf,Mek,Erk,Akt,Jnk,P38)
#                                   P(Plcg,PIP3,PIP2,PKC,Raf,    Erk,Akt,Jnk,P38|do(Mek))"
do_search_sachs_obs_intv_data <- "P(Plcg,PIP3,PIP2,PKC,Raf,    Erk,Akt,Jnk,P38|do(Mek))"

start = Sys.time()
dosearch(data = do_search_sachs_obs_intv_data,
         query = "P(Erk|do(Mek, PKC))",
         graph = sachs.dag_validated)
end = Sys.time()
end - start

# Estimating P(Erk = 1|do(Raf=1, Mek=1)) with the plug-in estimator
## PIP2 = 0, PIP3 = 0, PKC = 0, Plcg = 0
set.seed(1)
selected_rows = sample(nrow(isachs_MekIntv),300)
isachs_MekIntv_selected_rows <- isachs_MekIntv[selected_rows,]
isachs_MekIntv_select_held_out <- isachs_MekIntv[-selected_rows,]
P_erk.do.Raf1.Mek1 <- function(erk) {
  result = 0
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:3) {
        for (l in 1:3) {
          P_Mek.PKC.given.PIP2.PIP3.Plcg = nrow(filter(isachs_MekIntv_selected_rows, PIP2 == i, PIP3 == j, PKC == k,Plcg == l))/nrow(filter(isachs_MekIntv_selected_rows, PIP2 == i, PIP3 == j,Plcg == l))
          P_Mek.PIP2.given.PIP3.Plcg = nrow(filter(isachs_MekIntv_selected_rows, PIP2 == i, PIP3 == j,Plcg == l))/nrow(filter(isachs_MekIntv_selected_rows, PIP3 == j,Plcg == l))
          P_Mek.PIP3.given.Plcg = nrow(filter(isachs_MekIntv_selected_rows, PIP3 == j,Plcg == l))/nrow(filter(isachs_MekIntv_selected_rows,Plcg == l))
          P_Mek.Erk.given.PIP2.PIP3.PKC.Plcg = nrow(filter(isachs_MekIntv_selected_rows, Erk == erk, PIP2 == i, PIP3 == j, PKC == k,Plcg == l))/nrow(filter(isachs_MekIntv_selected_rows, PIP2 == i, PIP3 == j, PKC == k, Plcg == l))
          P_Mek.Plcg = nrow(filter(isachs_MekIntv_selected_rows, Plcg == l))/nrow(isachs_MekIntv_selected_rows)
          
          prod = P_Mek.PKC.given.PIP2.PIP3.Plcg * P_Mek.PIP2.given.PIP3.Plcg * P_Mek.PIP3.given.Plcg * P_Mek.Erk.given.PIP2.PIP3.PKC.Plcg * P_Mek.Plcg
          if (is.na(prod)) {
            next
          }
          result = result + prod
        }
      }
    }
  }
  return(result)
}

# Estimating P[Erk = 1 |do(Raf=1, Mek=1)] with the plug-in estimator using training data
P_erk.do.Raf1.Mek1(erk = 1)

# Estimating P(Erk = 1|do(Raf=1, Mek=1)) from P(Erk=1|do(Mek=1)) using heldout data
sum(as.numeric(isachs_MekIntv_select_held_out$Erk) == 1)/length(isachs_MekIntv_select_held_out$Erk)
