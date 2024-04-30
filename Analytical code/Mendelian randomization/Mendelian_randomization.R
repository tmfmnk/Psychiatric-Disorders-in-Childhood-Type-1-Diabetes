## ---------------------------
##
## Script name: Mendelian Randomization T1D - Psychiatric Disorders
##
## Purpose of script: Confirmatory Two-Sample MR for Observational Analysis using Czech Register Data
##                    
## Author: B Perry
##
## Version Date: 2023-03-30
##
## Notes: Amendments made based on comments from SB
##        
## ---------------------------

# load packages
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)
library(LDlinkR)
library(data.table)
library(export)
library(ieugwasr)
library(genetics.binaRies)

# set working directory
setwd("path/")


# create places to store results and graphs
dir.create("proxies1")
dir.create("scatterplots1")
dir.create("leaveoneout1")
dir.create("singlesnp1")
dir.create("Fstats1")
proxies<-NULL
harmonised<-NULL
MR_analysis_main<-NULL
bidirect_out<-NULL
bidirect_T1D<-NULL
bidirect_harmonised<-NULL
bidirect_results<-NULL



###### ###### ###### ###### ######
##### Define Cis Instruments ######
###### ###### ###### ###### ######

CTSH<-c("rs12592898", "rs60254670", "rs12148472", "rs34843303", "rs34593439","rs2289702")

GLIS3<-c("rs4380994", "rs3892354", "rs1574285", "rs10974435", "rs34494309", "rs57884925",
         "rs10758591", "rs7024686", "rs7041847", "rs7034200", "rs10814914","rs10116772",
         "rs10814915","rs6476839", "rs6476842","rs7020673","rs10974438", "rs10758593",
         "rs7867224", "rs10814916", "rs34706136", "rs1075859", "rs4339696", "rs10814917")

THEMIS<-c("rs67707912", "rs13204742", "rs6939", "rs9491889", "rs9491890", "rs9491891", 
          "rs147626184", "rs118097399", "rs9491892", "rs9482848" ,"rs9491893", "rs113297984",
          "rs72973797", "rs72973800", "rs761332", "rs9482849", "rs12111314", "rs11753289",
          "rs9482850", "rs9482851", "rs72975913", "rs72975916", "rs7738609", "rs138300818",
          "rs3901020", "rs4510698")

IL10<-c("rs3024505", "rs3024495", "rs3024493", "rs3122605")

IL2RA<-c("rs12722563", "rs12722558", "rs12722552", "rs12722522",
         "rs12722508", "rs7909519", "rs61839660", "rs12722496", 
         "rs12722495", "rs79092647", "rs41295049", "rs41295061",
         "rs41295065", "rs35285258", "rs11594656", "rs6602437",
         "rs34975410")

IKZF3<-c("rs2941522", "rs12946510", "rs72538185", "rs907091",
         "rs907092", "rs2952140", "rs2313430", "rs10445308",
         "rs12942330", "rs11658993", "rs2952144", "rs4795395",
         "rs9909593","rs71152606","rs9303277", "rs3816470", "rs4795397",
         "rs11557466", "rs11078925",  "rs34120102",  "rs11655198",  "rs11650661",
         "rs11655292", "rs12709365", "rs13380815", "rs11557467", "rs12936231",
         "rs11870965", "rs9903250", "rs9905959", "rs11658278", "rs10852935",
         "rs10852936", "rs9891174", "rs59716545", "rs36095411", "rs12939457",
         "rs367998020","rs34189114","rs35736272", "rs1054609", "rs9907088",
         "rs36038753","rs35569035", "rs9910826" ,"rs71355426", "rs9904624",
         "rs4795398", "rs12939565", "rs12939566", "rs71152620","rs12232497",
         "rs12232498", "rs12941333", "rs2872507", "rs9908132", "rs9901146",
         "rs12936409","rs12103884", "rs9906951", "rs12950209","rs12950743",
         "rs7359623",  "rs68122720",  "rs8067378", "rs12453507", "rs34170568",
         "rs11651596", "rs12949100", "rs8069176", "rs883770","rs11078926",
         "rs4795399","rs2305480", "rs2305479", "rs35196450", "rs56750287",
         "rs62067034","rs36000226",  "rs36084703", "rs11078927", "rs11078928",
         "rs12939832", "rs2290400", "rs1008723","rs4795400","rs869402", "rs1011082",
         "rs921650","rs921649",  "rs6503524", "rs7216389",  "rs7216558",  "rs150597688",
         "rs143385463",  "rs1031458", "rs1031460", "rs8065777","rs7219923",
         "rs7224129", "rs8074437","rs4065275", "rs12603332")
         
         
         
         

###### ###### ###### ###### ######
##### 1) load in exposure ######
###### ###### ###### ###### ######

T1D_GWAS<-read_exposure_data(
  filename = "MR_Analysis/Exposures/34012112-GCST90014023-EFO_0001359-Build38.f.tsv",
  sep = "\t",
  phenotype_col = "Type 1 Diabetes",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  eaf_col = "effect_allele_frequency",
  samplesize_col = NA) %>%
  mutate(beta.exposure = (beta.exposure - mean(beta.exposure))/sd(beta.exposure))



exposures<-list(
  
  # trans instrument - all genome-wide significant SNPs clumped for independence that are in LD ref panel
  "T1D_trans_dat"=T1D_GWAS %>%
    subset(pval.exposure<5e-8) %>% 
    clump_data(clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR") %>% mutate(exposure = "T1D_trans",
                                                                                                       id.exposure = "T1D_trans") %>% ld_reflookup(),
  # and the cis instruments
  "T1D_CTSH_dat"= T1D_GWAS %>% 
    filter(SNP %in% CTSH) %>%
    mutate(exposure = "T1D_CTSH",
           id.exposure = "T1D_CTSH") %>% ld_reflookup(SNP),
  
  "T1D_GLIS3_dat"= T1D_GWAS %>% 
    filter(SNP %in% GLIS3) %>%
    mutate(exposure = "T1D_GLIS3",
           id.exposure = "T1D_GLIS3") %>% ld_reflookup(),  
  
  "T1D_IKZF3_dat" = T1D_GWAS %>%
    filter(SNP %in% IKZF3) %>%
    mutate(exposure="T1D_IKZF3",
           id.exposure="T1D_IKZF3"),
  
  
  "T1D_IL10_dat" = T1D_GWAS %>%
    filter(SNP %in% IL10) %>%
    mutate(exposure="T1D_IL10",
           id.exposure="T1D_IL10"),
  
  "T1D_IL2RA_dat" = T1D_GWAS %>%
    filter(SNP %in% IL2RA) %>%
    mutate(exposure="T1D_IL2RA",
           id.exposure="T1D_IL2RA"),
  
  "T1D_THEMIS_dat" = T1D_GWAS %>%
    filter(SNP %in% THEMIS) %>%
    mutate(exposure="T1D_THEMIS",
           id.exposure="T1D_THEMIS")
  )



###### ###### ###### ###### ###### ###### ############ ############ ############ ###### ###### ####
##### 2) Define outcomes (pre-prepared using read_outcome_data() then saved as .txt)  ######
###### ###### ###### ###### ############ ############ ############ ############ ############ ###### 

outcomes<-c("AN", "AUD", "BPAD", "DEP",  "GAD", "SCZ")


###### ###### ###### ###### 
##### 3) Run Script  ######
###### ###### ###### ###### 


for (k in outcomes) {
  
  cat("\n", "loading in", k, "...")

  # read in outcome data
  temp_outcome <- fread(paste0("MR_Analysis/Outcomes/", k,".txt"), header =T, stringsAsFactors = F, data.table = F, fill =T)
  
  
  
  # harmonise data
  cat("\n","harmonising exposures and", k, "...")
  
  
  for (j in names(exposures)){
    temp_exp<-exposures[[j]]
    harmonised[[k]][[j]]<-harmonise_data(temp_exp, temp_outcome)
    rm(temp_exp)}
  
  
  
  # find and replace missing SNPs
  
  cat("\n","find and replace missing exposure SNPs for GWAS of ", k, "...")
  
  
  missing_SNP<-NULL
  for (j in names(exposures)){
    
    temp_all<-exposures[[j]]
    temp_got<-harmonised[[k]][[j]]
    missing_SNP[[k]][[j]]<-setdiff(temp_all$SNP, temp_got$SNP)
    if (length(missing_SNP[[k]][[j]])>0){
    missing_SNP[[k]][[j]]<-ld_reflookup(missing_SNP[[k]][[j]])}
    rm(temp_all, temp_got)}
  
  
  
  # create folder to store proxies
  setwd("path/proxies1")
  
  # store proxy SNPs
  
  for (j in names(exposures)){
    if (length(missing_SNP[[k]][[j]])>0){
      dir.create(paste0(j, "+", k, sep=""))
      setwd(paste0("path/proxies1/", j,"+" ,k, sep=""))
      LDproxy_batch(snp = missing_SNP[[k]][[j]],
                    pop ="CEU",
                    r2d ="r2",
                    token = "be25d7e95d26", 
                    genome_build = "grch37")
      setwd("path/proxies1")
    } 
  }

  
  # load in those stored proxy files, find matching proxy SNPs and keep the one with the lowest p-value
  
  for (j in names(exposures)){
    if (length(missing_SNP[[k]][[j]])>0){
      setwd(paste0("path/proxies1/", j, "+", k, sep=""))
      temp_SNP<-missing_SNP[[k]][[j]]
      for (i in temp_SNP){
        proxies[[k]][[j]][[i]]<-fread(paste0(i, "_grch37.txt"), head=T, stringsAsFactors = F, data.table = F, fill = T) %>%
          select(c(Coord, RegulomeDB))
        proxies[[k]][[j]][[i]]<-filter(temp_outcome, SNP %in% proxies[[k]][[j]][[i]][['Coord']])
        proxies[[k]][[j]][[i]]<-filter(T1D_GWAS, SNP %in% proxies[[k]][[j]][[i]][['SNP']]) %>%
          slice_min(pval.exposure) %>%
          slice(1) %>%
          mutate(exposure=paste0(exposures[[i]]["exposure"][1,], sep=""),
                 id.exposure=paste0(exposures[[i]]["id.exposure"][1,], sep="")) %>%
          mutate_at(c('beta.exposure', 'se.exposure', 'pval.exposure', 'eaf.exposure'), as.numeric)
      }
      rm(temp_SNP)
      setwd("path/proxies1")
    } 
  }
  
  #  Align alleles and use estimates from missing SNPs
  for (i in names(exposures)){
    if (length(missing_SNP[[k]][[i]])>0){
      temp_SNP<-filter(T1D_GWAS, SNP %in% missing_SNP[[k]][[i]])
      for (j in 1:nrow(temp_SNP)){
        temp_EA<-proxies[[k]][[i]][[j]]$effect_allele.exposure
        temp_OA<-proxies[[k]][[i]][[j]]$other_allele.exposure
        if (proxies[[k]][[i]][[j]]$beta.exposure>0 & temp_SNP[j,]$beta.exposure <0 |
            proxies[[k]][[i]][[j]]$beta.exposure<0 & temp_SNP[j,]$beta.exposure>0) {
          proxies[[k]][[i]][[j]]$effect_allele.exposure<-temp_OA
          proxies[[k]][[i]][[j]]$other_allele.exposure<-temp_EA}
        proxies[[k]][[i]][[j]]$beta.exposure<-temp_SNP[j,]$beta.exposure
        proxies[[k]][[i]][[j]]$se.exposure<-temp_SNP[j,]$se.exposure
        proxies[[k]][[i]][[j]]$pval.exposure<-temp_SNP[j,]$pval.exposure
        proxies[[k]][[i]][[j]]$exposure<-exposures[[i]]$exposure[1]
        proxies[[k]][[i]][[j]]$id.exposure<-exposures[[i]]$id.exposure[1]
        
        rm(temp_EA, temp_OA)}}
  }
  
  
  # harmonise proxy SNPs and add them to main harmonised dataframes
  
  cat("\n","harmonising proxy SNPs", "...", "\n")

  for (i in names(exposures)){
    if (length(missing_SNP[[k]][[i]])>0){
      for (j in 1:length(missing_SNP[[k]][[i]])){
        proxies[[k]][[i]][[j]]<-harmonise_data(proxies[[k]][[i]][[j]], temp_outcome)}
      proxies[[k]][[i]]<-as.data.frame(do.call(rbind, proxies[[k]][[i]])) 
      rownames(proxies[[k]][[i]])<-NULL
      harmonised[[k]][[i]]<-rbind(harmonised[[k]][[i]], proxies[[k]][[i]]) %>%
        mutate_at(c('beta.exposure', 'se.exposure', 'pval.exposure', 'eaf.exposure',
                    'beta.outcome', 'se.outcome', 'pval.outcome', 'samplesize.outcome'), as.numeric) %>%
        distinct(SNP, .keep_all = T)
    }}
  
  setwd("path/")

    
  # Run Normal MR for Trans instrument, correlation adjusted for cis instruments 
  
  setwd("path/")
  
  
  for (j in names(exposures)){
    temp_harm<-data.frame(harmonised[[k]][[j]])
    cat("starting analysis for", j, "on", k, "...", "\n", "main analysis...", "\n")

    MR_analysis_main[[k]][[j]][["main"]]<-mr(as.data.frame(temp_harm)) %>% generate_odds_ratios()
    MR_analysis_main[[k]][[j]][["heterogeneity"]]<-mr_heterogeneity(as.data.frame(temp_harm))
    MR_analysis_main[[k]][[j]][["pleiotropy"]]<-mr_pleiotropy_test(as.data.frame(temp_harm))
    cat("\n", "doing mrpresso [this can take some time]...", "\n")
    MR_analysis_main[[k]][[j]][["MRPRESSO"]]<-mr_presso(BetaOutcome = "beta.outcome",
                                                        BetaExposure = "beta.exposure",
                                                        SdOutcome = "se.outcome",
                                                        SdExposure = "se.exposure",
                                                        OUTLIERtest = TRUE,
                                                        DISTORTIONtest = TRUE,
                                                        NbDistribution = 5000,
                                                        data = as.data.frame(temp_harm))
    cat("\n", "doing plots...", "\n")

    MR_analysis_main[[k]][[j]][["scatter"]]<-mr_scatter_plot(MR_analysis_main[[k]][[j]][["main"]], as.data.frame(temp_harm)) %>% print()
    graph2vector( x = NULL, file = paste0("scatterplots1/", k, "-", names(exposures[j]), "_scatterplot.pdf"),
                  fun = NULL,type = "PDF", aspectr = NULL, width = 5, height = 5, scaling = 100, font = ifelse(Sys.info()["sysname"] == "Windows", "Arial", "Helvetica")[[1]],
                  bg = "white", colormodel = "rgb", cairo = TRUE, fallback_resolution = 600)
    dev.off()
    MR_analysis_main[[k]][[j]][["loo"]]<-mr_leaveoneout(as.data.frame(temp_harm)) %>% mr_forest_plot(exponentiate = T) %>% print()
    graph2vector( x = NULL, file = paste0("leaveoneout1/", k, "-", names(exposures[j]), "_leaveoneout.pdf"),
                  fun = NULL,type = "PDF", aspectr = NULL, width = 5, height = 5, scaling = 100, font = ifelse(Sys.info()["sysname"] == "Windows", "Arial", "Helvetica")[[1]],
                  bg = "white", colormodel = "rgb", cairo = TRUE, fallback_resolution = 600)
    dev.off()
    MR_analysis_main[[k]][[j]][["singlesnp"]]<-mr_singlesnp(as.data.frame(temp_harm)) %>% mr_forest_plot(exponentiate = T) %>% print()
    graph2vector( x = NULL, file = paste0("singlesnp1/", k, "-", names(exposures[j]), "_singlesnp.pdf"),
                  fun = NULL,type = "PDF", aspectr = NULL, width = 5, height = 5, scaling = 100, font = ifelse(Sys.info()["sysname"] == "Windows", "Arial", "Helvetica")[[1]],
                  bg = "white", colormodel = "rgb", cairo = TRUE, fallback_resolution = 600)
    dev.off()
    cat("\n", "finishing sensitivity analyses", "\n")

    cat("\n")
    MR_analysis_main[[k]][[j]][["instrumentstrength"]][["Fstat"]]<-data.frame("SNP"=exposures[[j]][["SNP"]],
                                                                              "Fstat"=exposures[[j]][["beta.exposure"]]^2/exposures[[j]][["se.exposure"]]^2)
    MR_analysis_main[[k]][[j]][["instrumentstrength"]][["median"]]<-median(MR_analysis_main[[k]][[j]][["instrumentstrength"]][["Fstat"]][["Fstat"]])
    MR_analysis_main[[k]][[j]][["I2GX"]]<-Isq(abs(exposures[[j]][["beta.exposure"]]), exposures[[j]][["se.exposure"]])
    if (j != "T1D_trans_dat"){
      cat("\n", "cis instrument detected... doing correlation adjusted IVW...", "\n")

      MR_temp<-dat_to_MRInput(temp_harm, get_correlations = T)
      MR_analysis_main[[k]][[j]][["corr_adjusted"]]<-MendelianRandomization::mr_ivw(MR_temp[[1]], correl = T)
      rm(MR_temp)}
    rm(temp_harm)
    write.csv(MR_analysis_main[[k]][[j]][["instrumentstrength"]][["Fstat"]], file=paste0("Fstats1/",k, "-", names(exposures[j]), "_Fstats.csv"),
              append=F, row.names = F)

  }
  
  setwd("path/")
  
  
  # Bidirectional MR 
  
  # make trans variants for outcomes (bidirectional analysis), harmonise and run MR
  
  cat("\n", "finding trans variants for bidirectional analysis of", k, "\n")

  bidirect_out<-temp_outcome %>%
    subset(pval.outcome<5e-8) %>%
    clump_data(clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR") %>% 
    mutate(exposure = paste0(k), id.exposure = paste0(k))
  names(bidirect_out)<-str_replace(names(bidirect_out), "outcome", "exposure")
  bidirect_T1D<-filter(T1D_GWAS, SNP %in% bidirect_out[["SNP"]])
  names(bidirect_T1D)<-str_replace(names(bidirect_T1D), "exposure", "outcome")
  temp_exp<-bidirect_out
  temp_out<-bidirect_T1D
  cat("\n", "harmonising data for bidirectional analysis of", k, "\n")
  bidirect_harmonised[[k]]<-harmonise_data(temp_exp, temp_out)
  rm(temp_exp, temp_out)
  cat("\n", "performing bidirectional analysis of", k, "\n")
  bidirect_results[[k]]<-mr(bidirect_harmonised[[k]]) %>% generate_odds_ratios()
  
  rm(temp_outcome)
  cat("\n", "done analysis of", k, "\n")

}

# Check warnings
# NB currently unable to suppress ggplot-related warnings from mr_forest_plot()
# because suppressWarnings() doesn't work :(
warnings()

# View Results
View(MR_analysis_main)
View(bidirect_results)

# manually adjust p-values for multiple testing
p.adjust(p=c(), method="holm", n=length(p))
