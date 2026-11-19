#install.packages("C:/Users/user/Downloads/httk_2.3.1.tar.gz", repos = NULL, type = "source")

## Supplemental Material for "Concentrations, Composition Profiles, 
## and In Vitro–In Silico-Based Mixture Risk Assessment 
## of Bisphenol A and its Analogs in Plant-Based Foods"

## Load the httk package
## This study used httk version 2.3.1. 
## To reproduce the results, please ensure you are using this version.
## If httk 2.3.1 is not installed, download and install "httk_2.3.1.tar.gz" from CRAN.

library(httk)
packageVersion("httk")     ## Check httk version

## List CAS numbers and names of the 8 bisphenols (BPs) included in this study
BPs_info <- list(
  "BPs" = c("2,2-BPF ", "4,4-BPF", "BPA", "BPAF", 
            "BPB ", "BPC", "BPE", "BPS"),
  "CAS" = c("2467-02-9", "620-92-8", "80-05-7", "1478-61-1", 
            "77-40-7", "79-97-0", "2081-08-5", "80-09-1"))

## Check the availability of HTTK data for 8 BPs before loading Sipes et al. (2017) data
## Return result: FALSE
check_beforeload <- all(BPs_info$CAS %in% get_cheminfo())

## Load CLint and Fup QSPR predictions from Sipes et al. (2017)
chem.physical_and_invitro.data <- load_sipes2017(overwrite = FALSE)

## Check the availability of HTTK data for the 8 BPs after loading
## Return result: TRUE
check_afterload <- all(BPs_info$CAS %in% get_cheminfo())

## View chemical parameters for 8 BPs from chem.physical_and_invitro.data
library(dplyr)
params_BPs <- chem.physical_and_invitro.data %>% filter(CAS %in% BPs_info$CAS)


## Retrieve in vitro measured Caco-2 data for each BP
getCaco2 <- function(chem_cas, chem_name) {
  params <- get_caco2(chem.cas = chem_cas) 
  params_df <- as.data.frame(params)              
  params_df$Chem.CAS <- chem_cas
  return(params_df)
}

## Retrieve fraction of chemical absorbed from the gut
getfabsgut <- function(chem_cas) {
  params <- get_fabsgut(chem.cas = chem_cas) 
  params_df <- as.data.frame(params)              
  params_df$Chem.CAS <- chem_cas
  return(params_df)
}

## Merge Caco-2 and gut absorption parameters into a single dataset
combined_params_Caco2 <- do.call(rbind, lapply(BPs_info$CAS, getCaco2))
combined_params_fabsgut <- do.call(rbind, lapply(BPs_info$CAS, getfabsgut))
combined_Caco2_fabsgut <- full_join(
  combined_params_Caco2, combined_params_fabsgut, by = "Chem.CAS")

## Retrieve PBTK model parameters for each BP
getpbtkparams <- function(chem_cas) {
  params <- parameterize_pbtk(chem.cas = chem_cas) 
  params_df <- as.data.frame(params)              
  params_df$Chem.CAS <- chem_cas
  return(params_df)
}
combined_params_pbtk <- do.call(rbind, lapply(BPs_info$CAS, getpbtkparams))


## Monte Carlo simulations for Css: Loop over all 8 BPs
BPs.Css <- c()
clearance <- c()
BPs_Css_table <- as.data.frame(BPs_info)
for (this.cas in BPs_info$CAS) {
  set.seed(3112)
  BPs.Css.mc <- calc_mc_css(chem.cas = this.cas, output.units = "uM", 
                            model = "pbtk", species = "human",
                            daily.dose = 1, concentration = "plasma", 
                            invitrouv = TRUE,
                            samples = 10000, return.samples = TRUE,
                            httkpop.generate.arg.list = list(
                            method = "vi",
                            gendernum = NULL,
                            agelim_years = c(20,50))
  )
  
  BPs.name <- BPs_Css_table[BPs_Css_table$CAS == this.cas, "BPs"] 
  BPs.Css.samples <- as.data.frame(BPs.Css.mc)
  names(BPs.Css.samples) <- BPs.name
  BPs.Css <- c(BPs.Css, BPs.Css.samples)
  
  ## Compute summary statistics for Css
  percentile <- quantile(BPs.Css.mc, probs = c(0.025,0.25,0.5,0.75,0.975), na.rm = TRUE)
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"Mean"] <- mean(BPs.Css.mc)
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"SD"] <- sd(BPs.Css.mc)
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"P2.5"] <- percentile[1]
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"P25"] <- percentile[2]
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"P50"] <- percentile[3]
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"P75"] <- percentile[4]
  BPs_Css_table[BPs_Css_table$CAS == this.cas,"P97.5"] <- percentile[5]
  
  ## calculate the total plasma clearance
  clearance.calculation <- calc_total_clearance(chem.cas = this.cas, species = "Human")
  clearance.value <- as.data.frame(clearance.calculation)
  names(clearance.value) <- BPs.name
  clearance <- c(clearance, clearance.value)
}


## Save Css statistics to a CSV file
write.csv(BPs_Css_table, file = "BPs_Css_statistics.csv", row.names = FALSE)

## Save all Css samples to a CSV file
write.csv(BPs.Css, file = "BPs_Css_allsamples.csv", row.names = FALSE)

## Save clearance to a CSV file
write.csv(clearance, file = "BPs_Clearance.csv", row.names = FALSE)

## Save PBTK model parameters to a CSV file 
write.csv(combined_params_pbtk, file = "BPs_Pbtk_params.csv", row.names = FALSE)


