library(bnlearn)
library(dagitty)
library(dplyr)
library(dosearch) #can check if query is identifiable from any combination of observational and interventional distributions
source("/Users/sarataheri/GitHub/OptimalAdjustmentSet/main_functions.R")
# Prior knowledge on the network
covid.prior.dag <- dagitty('dag {
                              SARS_COV2 -> ACE2
                              SARS_COV2 -> PRR
                              ACE2 -> Ang
                              Ang -> AGTR1
                              AGTR1 -> ADAM17
                              ADAM17 -> EGF
                              ADAM17 -> TNF
                              ADAM17 -> Sil6r
                              PRR -> NFKB
                              Gefi -> EGFR
                              EGF -> EGFR
                              TNF -> NFKB
                              Sil6r -> IL6STAT3
                              Toci -> Sil6r
                              NFKB -> IL6AMP
                              EGFR -> NFKB
                              IL6STAT3 -> IL6AMP
                              IL6AMP -> cytok
                              SARS_COV2 <-> Ang
                              ADAM17 <-> Sil6r
                              EGFR <-> IL6STAT3
                              EGFR <-> EGF
                              EGFR <-> TNF
                              PRR <-> NFKB
                            }')
covid_obs_data <- read.csv("Covid_case_study/Covid_obs_1.csv")
covid_obs_data <- covid_obs_data[,names(covid.prior.dag)]

K=100
rep_full_lm_N1000 <- readRDS("/Users/sarataheri/Desktop/Simplified_LVM_not_Shared/Covid_case_study/data/Covid_hmc_fit/full_lm_ndpu1000.RData")[1:K]
rep_full_aipw_N1000 <- unname(as.vector(read.csv("/Users/sarataheri/Desktop/Simplified_LVM_not_Shared/Covid_case_study/data/Covid_hmc_fit/full_aipw_ndpu1000.csv")))[[1]][1:K]
rep_full_ipw_N1000 <- unname(as.vector(read.csv("/Users/sarataheri/Desktop/Simplified_LVM_not_Shared/Covid_case_study/data/Covid_hmc_fit/full_ipw_ndpu1000.csv")))[[1]][1:K]
rep_full_gformula_N1000 <- unname(as.vector(read.csv("/Users/sarataheri/Desktop/Simplified_LVM_not_Shared/Covid_case_study/data/Covid_hmc_fit/full_gformula_ndpu1000.csv")))[[1]][1:K]

df_ananke <- data.frame("CE" = c(rep_full_aipw_N1000,
                                 rep_full_gformula_N1000,
                                 rep_full_ipw_N1000,
                                rep_full_lm_N1000
                                 ),
                        "model" = c(rep("AIPW approach", K),
                                    rep("gformnula approach", K),
                                    rep("IPW approach", K),
                                    rep("Linear regression", K)
                        )
)

gg2 <- df_ananke %>% ggplot(aes(x = as.factor(model), y = CE, fill=model)) +
  geom_boxplot(width = 0.4) +
  geom_abline(intercept=0.79, slope=0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_blank(),#element_text(size = 35,face="bold"),
        legend.key.size = unit(2, "cm"),
        axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=35,face="bold"),
        axis.title.x = element_text(size=35,face="bold", #vjust = -0.3
        ),
        axis.text.x = element_text(#angle = 90, vjust = 0.5, hjust=1,
          size = 35, face="bold"),
        legend.position = "none",#"top",
        plot.title = element_text(size = 35, hjust = 0.5, face = "bold")
  ) +
  scale_fill_brewer(palette="BuPu") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  xlab("Different estimators") +
  ylim(-1.5,3) + # for ananke
  #ylim(-0.5,2) + # for lm
  ylab(expression(bold("ATE"))) +
  scale_x_discrete(labels = c('AIPW', 'gformula', 'IPW', 'Linear \n regression'))

gg2
ggsave("/Users/sarataheri/GitHub/Simplified_LVM/img/Covid_boxplots.pdf", plot = gg2, width = 10, height = 10, dpi = 300, units = "in")
