library(tidyverse)


#ecoli = read_csv("obs_data_ecoli_without_zeros.csv")


intervention = function(model,data,cause,do){
  newdata = data
  newdata[[cause]] = do
  
  Ey_dox = mean(predict(model,newdata=newdata))
  return(Ey_dox)
}



exhaustive_backdoor_lms = function(cause,effect,default_adj,potential_adj,data,do){
  
  default_vars = c(cause,default_adj)
  
  default_formula = paste(effect,
                          paste(default_vars, collapse = " + "),
                          sep = " ~ ")
  
  my_combi <- unlist(lapply(1:length(potential_adj),    # Get all combinations
                            combinat::combn, 
                            x = potential_adj,
                            simplify = FALSE), 
                     recursive = FALSE)
  
  my_combi_vars = lapply(my_combi,function(x) c(default_vars,x))
  my_combi_vars[[length(my_combi_vars)+1]] = default_vars #adding the default in 
  
  all_forms = sapply(my_combi_vars,function(x) paste(effect,
                                                     paste(x, collapse = " + "),
                                                     sep = " ~ "))
  
  models = lapply(all_forms,function(x) lm(x,data=data))
  
  intervention_result = sapply(models,function(x) intervention(x,data=data,cause=cause,do=do))
  
  intervention_result = tibble(formula=all_forms,Ey_dox=intervention_result)
  return(intervention_result)
  
}

#exhaustive_backdoor_lms("fur","dpiA",default_adj = c("crp","rpoD"),potential_adj = c("oxyR","narL"),data=ecoli,do=1) 
