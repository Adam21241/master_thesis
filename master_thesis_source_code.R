# setwd("/home/adam/Dokumenty/studia/PRACA_MGR/TARGET/wyniki/kod ostateczny/")


### wczytanie bazy danych

t <- readRDS(gzcon(url("http://zbo.ipipan.waw.pl/files/data/OrigCecRNA_en.rds")))

# saveRDS(t, "baza_danych.RDS")
# t <- readRDS("baza_danych.RDS")

### podział na S i L (tu generuje się rysunek 1)

target <- subset(t, subtype_Survival..months. >= 0)
x <- target$subtype_Survival..months.
x <- x[complete.cases(x)]

library(ash)
library(Cairo)

cairo_pdf(file = "rysunek_1.pdf")
h <- hist(x, nclass = 300, plot=FALSE)
plot(h, main = "", xlab = "Przeżywalność (w miesiącach)",
     ylab = "Udział procentowy", xlim = c(0,250), ylim = c(0,0.1), freq=FALSE)
f <- ash1(bin1(x,nbin=300),5)
lines(f , type="l", col = "blue") 
lines(density(x, bw = "sj", kernel = "cosine"), col = "red")
legend(x = 100, y = 0.08, legend = c("Averaged Shifted Histogram", "Kernel Density Estimation"), col = c("blue", "red"), lty = 1)
dev.off()


### ANALIZA 1 (tylko ekspresje genów): MCFS

target <- t

to.remove <- c("barcode", "sample", "shortLetterCode", "definition", "classification_of_tumor", "last_known_disease_status", "updated_datetime.x",
               "primary_diagnosis", "tumor_stage", "age_at_diagnosis", "morphology", "days_to_death", "days_to_last_known_disease_status",
               "created_datetime.x", "state.x", "days_to_recurrence", "diagnosis_id", "tumor_grade", "treatments", "tissue_or_organ_of_origin",
               "days_to_birth", "progression_or_recurrence", "prior_malignancy", "site_of_resection_or_biopsy", "days_to_last_follow_up",
               "cigarettes_per_day", "weight", "updated_datetime.y", "alcohol_history", "alcohol_intensity", "bmi", "years_smoked", "created_datetime.y",
               "state.y", "exposure_id", "height", "updated_datetime", "created_datetime", "gender", "year_of_birth", "state", "race", "demographic_id",
               "ethnicity", "year_of_death", "bcr_patient_barcode", "dbgap_accession_number", "disease_type", "released", "state.1", "primary_site",
               "project_id", "name", "subtype_patient", "subtype_Tissue.source.site", "subtype_Study", "subtype_Whole.genome", "subtype_SNP6",
               "subtype_U133a", "subtype_HM450", "subtype_HM27", "subtype_RPPA", "subtype_Histology", "subtype_Grade", "subtype_Age..years.at.diagnosis.",
               "subtype_Gender", "subtype_Vital.status..1.dead.", "subtype_Karnofsky.Performance.Score",
               "subtype_Mutation.Count", "subtype_Percent.aneuploidy", "subtype_IDH.status", "subtype_X1p.19q.codeletion", "subtype_IDH.codel.subtype",
               "subtype_MGMT.promoter.status", "subtype_Chr.7.gain.Chr.10.loss", "subtype_Chr.19.20.co.gain", "subtype_TERT.promoter.status",
               "subtype_TERT.expression..log2.", "subtype_TERT.expression.status", "subtype_ATRX.status", "subtype_DAXX.status",
               "subtype_Telomere.Maintenance", "subtype_BRAF.V600E.status", "subtype_BRAF.KIAA1549.fusion", "subtype_ABSOLUTE.purity",
               "subtype_ABSOLUTE.ploidy", "subtype_ESTIMATE.stromal.score", "subtype_ESTIMATE.immune.score", "subtype_ESTIMATE.combined.score",
               "subtype_Original.Subtype", "subtype_Transcriptome.Subtype", "subtype_Pan.Glioma.DNA.Methylation.Cluster",
               "subtype_IDH.specific.DNA.Methylation.Cluster", "subtype_Supervised.DNA.Methylation.Cluster", "subtype_Random.Forest.Sturm.Cluster",
               "subtype_RPPA.cluster", "subtype_Telomere.length.estimate.in.blood.normal..Kb.", "subtype_Telomere.length.estimate.in.tumor..Kb.", "subtype_BCR",
               "subtype_Whole.exome", "therapeutic_agents", "submitter_id", "treatment_id", "days_to_treatment", "treatment_intent_type",
               "treatment_or_therapy", "subtype_RNAseq", "subtype_Pan.Glioma.RNA.Expression.Cluster", "subtype_IDH.specific.RNA.Expression.Cluster")

target <- target[ , !(names(target) %in% to.remove)]

target <- subset(target, (subtype_Survival..months. >= 1 & subtype_Survival..months. <= 13 & vital_status == "dead") | subtype_Survival..months. >= 39 ) 

target <- target[ , -c(1:2)]
target <- target[complete.cases(target), ]

dim(target)


library(rmcfs)

input <- fix.data(target)

sort(input$subtype_Survival..months.)
input$subtype_Survival..months.[1:105] <- "S"
input$subtype_Survival..months.[106:203] <- "L"
input$subtype_Survival..months. <- as.factor(input$subtype_Survival..months.)

set.seed(123) # przykładowo, przy takim ustawieniu ziarna (123), metoda sample losuje dobrze, a runif słabo

# mycv <- round(runif(nrow(input), 1, 3)) # przedziały nie są do końca równoliczne (np. 60, 60, 80)

s <- sample(c(1:3), nrow(input), prob = c(0.33, 0.34, 0.33), replace = TRUE)

saveRDS(input[s == 1, ], "input_1_1.RDS")  # konwencja: nr eksperymentu _ nr z s
saveRDS(input[s == 2, ], "input_1_2.RDS")  
saveRDS(input[s == 3, ], "input_1_3.RDS") 


#

results <- mcfs(subtype_Survival..months.~., input[s != 3,], featureFreq = 100, splitSetSize = 100, threadsNumber = 10)
saveRDS(results, "results_1_12.RDS")

results <- mcfs(subtype_Survival..months.~., input[s != 2,], featureFreq = 100, splitSetSize = 100, threadsNumber = 10)
saveRDS(results, "results_1_13.RDS")

results <- mcfs(subtype_Survival..months.~., input[s != 1,], featureFreq = 100, splitSetSize = 100, threadsNumber = 10)
saveRDS(results, "results_1_23.RDS")


### ANALIZA 2 (dodatkowo dane kliniczne): MCFS

target_2 <- t

target_2 <- target_2[colSums(is.na(target_2)) < 200] 

to.remove <- c("patient", "barcode", "sample", "classification_of_tumor", "last_known_disease_status", "tumor_stage", "updated_datetime.x", "diagnosis_id", 
               "tumor_grade", "state.x", "days_to_birth", "progression_or_recurrence", "prior_malignancy", "site_of_resection_or_biopsy", "updated_datetime.y", "state.y",
               "updated_datetime", "state", "shortLetterCode", "treatments", "primary_diagnosis", "bcr_patient_barcode", "released", "state.1", "primary_site", "project_id",
               "name", "subtype_patient", "subtype_Whole.exome", "updated_datetime.1", "submitter_id", "treatment_id", "state.2", "exposure_id", "race", 
               "demographic_id", "ethnicity", "definition", "subtype_Study", "subtype_RNAseq", "subtype_Chr.19.20.co.gain", "subtype_DAXX.status", "subtype_BRAF.V600E.status",
               "subtype_BRAF.KIAA1549.fusion", "subtype_U133a"
)

target_2 <- target_2[ , !(names(target_2) %in% to.remove)]
target_2 <- subset(target_2, (subtype_Survival..months. >= 1 & subtype_Survival..months. <= 13 & vital_status == "dead") | subtype_Survival..months. >= 39 ) 

target_2$disease_type <- as.character(target_2$disease_type)

target_2 <- target_2[complete.cases(target_2), ]

dim(target_2)


library(rmcfs)

input <- fix.data(target_2)

sort(input$subtype_Survival..months.)
input$subtype_Survival..months.[1:24] <- "S"
input$subtype_Survival..months.[25:71] <- "L"
input$subtype_Survival..months. <- as.factor(input$subtype_Survival..months.)

set.seed(123)

s <- sample(c(1:3), nrow(input), prob = c(0.33, 0.34, 0.33), replace = TRUE)

saveRDS(input[s == 1, ], "input_2_1.RDS") 
saveRDS(input[s == 2, ], "input_2_2.RDS") 
saveRDS(input[s == 3, ], "input_2_3.RDS") 


#

results <- mcfs(subtype_Survival..months.~., input[s != 3,], featureFreq = 100, splitSetSize = 100, threadsNumber = 10)
saveRDS(results, "results_2_12.RDS")

results <- mcfs(subtype_Survival..months.~., input[s != 2,], featureFreq = 100, splitSetSize = 100, threadsNumber = 10)
saveRDS(results, "results_2_13.RDS")

results <- mcfs(subtype_Survival..months.~., input[s != 1,], featureFreq = 100, splitSetSize = 100, threadsNumber = 10)
saveRDS(results, "results_2_23.RDS")


### ANALIZA 1: budowa modelu przewidującego przeżywalność pacjentów (tu generują się rysunki: 2, 3)

results_1_12 <- readRDS("results_1_12.RDS")
results_1_12$cutoff # algorytmy odcięcia

results_1_13 <- readRDS("results_1_13.RDS")
results_1_13$cutoff 

results_1_23 <- readRDS("results_1_23.RDS")
results_1_23$cutoff 


library(rmcfs)

pdf(file = "rysunek_2.pdf")
plot(results_1_12,  main = "", type = "features", cex = 0.35)
dev.off()

pdf(file = "rysunek_3.pdf")
plot(results_1_12, type = "cv", measure = "wacc", cex = 0.3)
dev.off()


input_1_1 <- readRDS("input_1_1.RDS")
input_1_2 <- readRDS("input_1_2.RDS")
input_1_3 <- readRDS("input_1_3.RDS")


# labels vs wacc

# ranking 12

wacc_IBk_1_12 <- c()
wacc_SMO_1_12 <- c()
wacc_Logistic_1_12 <- c()
wacc_JRip_1_12 <- c()
wacc_nb_1_12 <- c()
wacc_rf_1_12 <- c()

for(i in 1:results_1_12$cutoff$size[3]){
  
  temp <- rbind(input_1_1, input_1_2)
  train <- temp[, c(head(results_1_12$RI$attribute, i), "subtype_Survival..months.")]
  test <- input_1_3[, c(head(results_1_12$RI$attribute, i), "subtype_Survival..months.")]
  
  alg_IBk <- RWeka::IBk(subtype_Survival..months.~., train)   # knn
  testpred_IBk <- predict(alg_IBk, newdata = test)
  cm_IBk <- cbind(test$subtype_Survival..months., testpred_IBk)
  
  alg_SMO <- RWeka::SMO(subtype_Survival..months.~., train)
  testpred_SMO <- predict(alg_SMO, newdata = test)
  cm_SMO <- cbind(test$subtype_Survival..months., testpred_SMO)
  
  alg_Logistic <- RWeka::Logistic(subtype_Survival..months.~., train)
  testpred_Logistic <- predict(alg_Logistic, newdata = test)
  cm_Logistic <- cbind(test$subtype_Survival..months., testpred_Logistic)
  
  alg_JRip <- RWeka::JRip(subtype_Survival..months.~., train)
  testpred_JRip <- predict(alg_JRip, newdata = test)
  cm_JRip <- cbind(test$subtype_Survival..months., testpred_JRip)
  
  alg_nb <- e1071::naiveBayes(subtype_Survival..months.~., train)
  testpred_nb <- predict(alg_nb, newdata = test)
  cm_nb <- cbind(test$subtype_Survival..months., testpred_nb)
  
  alg_rf <- randomForest::randomForest(subtype_Survival..months.~., train)
  testpred_rf <- predict(alg_rf, newdata = test)
  cm_rf <- cbind(test$subtype_Survival..months., testpred_rf)
  
  
  wacc_IBk_1_12[i] <- ((length(subset(cm_IBk[,1], cm_IBk[,1] == 2 & cm_IBk[,2] == 2 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 2))) + 
                         (length(subset(cm_IBk[,1], cm_IBk[,1] == 1 & cm_IBk[,2] == 1 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 1)))) / 2
  
  wacc_SMO_1_12[i] <- ((length(subset(cm_SMO[,1], cm_SMO[,1] == 2 & cm_SMO[,2] == 2 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 2))) + 
                         (length(subset(cm_SMO[,1], cm_SMO[,1] == 1 & cm_SMO[,2] == 1 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 1)))) / 2
  
  wacc_Logistic_1_12[i] <- ((length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2 & cm_Logistic[,2] == 2 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2))) + 
                              (length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1 & cm_Logistic[,2] == 1 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1)))) / 2
  
  wacc_JRip_1_12[i] <- ((length(subset(cm_JRip[,1], cm_JRip[,1] == 2 & cm_JRip[,2] == 2 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 2))) + 
                          (length(subset(cm_JRip[,1], cm_JRip[,1] == 1 & cm_JRip[,2] == 1 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 1)))) / 2
  
  wacc_nb_1_12[i] <- ((length(subset(cm_nb[,1], cm_nb[,1] == 2 & cm_nb[,2] == 2 )) / length(subset(cm_nb[,1], cm_nb[,1] == 2))) + 
                        (length(subset(cm_nb[,1], cm_nb[,1] == 1 & cm_nb[,2] == 1 )) / length(subset(cm_nb[,1], cm_nb[,1] == 1)))) / 2
  
  wacc_rf_1_12[i] <- ((length(subset(cm_rf[,1], cm_rf[,1] == 2 & cm_rf[,2] == 2 )) / length(subset(cm_rf[,1], cm_rf[,1] == 2))) + 
                        (length(subset(cm_rf[,1], cm_rf[,1] == 1 & cm_rf[,2] == 1 )) / length(subset(cm_rf[,1], cm_rf[,1] == 1)))) / 2
  
}


# ranking 13

wacc_IBk_1_13 <- c()
wacc_SMO_1_13 <- c()
wacc_Logistic_1_13 <- c()
wacc_JRip_1_13 <- c()
wacc_nb_1_13 <- c()
wacc_rf_1_13 <- c()

for(i in 1:results_1_13$cutoff$size[3]){
  
  temp <- rbind(input_1_1, input_1_3)
  train <- temp[, c(head(results_1_13$RI$attribute, i), "subtype_Survival..months.")]
  test <- input_1_2[, c(head(results_1_13$RI$attribute, i), "subtype_Survival..months.")]
  
  alg_IBk <- RWeka::IBk(subtype_Survival..months.~., train)   # knn
  testpred_IBk <- predict(alg_IBk, newdata = test)
  cm_IBk <- cbind(test$subtype_Survival..months., testpred_IBk)
  
  alg_SMO <- RWeka::SMO(subtype_Survival..months.~., train)
  testpred_SMO <- predict(alg_SMO, newdata = test)
  cm_SMO <- cbind(test$subtype_Survival..months., testpred_SMO)
  
  alg_Logistic <- RWeka::Logistic(subtype_Survival..months.~., train)
  testpred_Logistic <- predict(alg_Logistic, newdata = test)
  cm_Logistic <- cbind(test$subtype_Survival..months., testpred_Logistic)
  
  alg_JRip <- RWeka::JRip(subtype_Survival..months.~., train)
  testpred_JRip <- predict(alg_JRip, newdata = test)
  cm_JRip <- cbind(test$subtype_Survival..months., testpred_JRip)
  
  alg_nb <- e1071::naiveBayes(subtype_Survival..months.~., train)
  testpred_nb <- predict(alg_nb, newdata = test)
  cm_nb <- cbind(test$subtype_Survival..months., testpred_nb)
  
  alg_rf <- randomForest::randomForest(subtype_Survival..months.~., train)
  testpred_rf <- predict(alg_rf, newdata = test)
  cm_rf <- cbind(test$subtype_Survival..months., testpred_rf)
  
  
  wacc_IBk_1_13[i] <- ((length(subset(cm_IBk[,1], cm_IBk[,1] == 2 & cm_IBk[,2] == 2 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 2))) + 
                         (length(subset(cm_IBk[,1], cm_IBk[,1] == 1 & cm_IBk[,2] == 1 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 1)))) / 2
  
  wacc_SMO_1_13[i] <- ((length(subset(cm_SMO[,1], cm_SMO[,1] == 2 & cm_SMO[,2] == 2 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 2))) + 
                         (length(subset(cm_SMO[,1], cm_SMO[,1] == 1 & cm_SMO[,2] == 1 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 1)))) / 2
  
  wacc_Logistic_1_13[i] <- ((length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2 & cm_Logistic[,2] == 2 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2))) + 
                              (length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1 & cm_Logistic[,2] == 1 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1)))) / 2
  
  wacc_JRip_1_13[i] <- ((length(subset(cm_JRip[,1], cm_JRip[,1] == 2 & cm_JRip[,2] == 2 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 2))) + 
                          (length(subset(cm_JRip[,1], cm_JRip[,1] == 1 & cm_JRip[,2] == 1 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 1)))) / 2
  
  wacc_nb_1_13[i] <- ((length(subset(cm_nb[,1], cm_nb[,1] == 2 & cm_nb[,2] == 2 )) / length(subset(cm_nb[,1], cm_nb[,1] == 2))) + 
                        (length(subset(cm_nb[,1], cm_nb[,1] == 1 & cm_nb[,2] == 1 )) / length(subset(cm_nb[,1], cm_nb[,1] == 1)))) / 2
  
  wacc_rf_1_13[i] <- ((length(subset(cm_rf[,1], cm_rf[,1] == 2 & cm_rf[,2] == 2 )) / length(subset(cm_rf[,1], cm_rf[,1] == 2))) + 
                        (length(subset(cm_rf[,1], cm_rf[,1] == 1 & cm_rf[,2] == 1 )) / length(subset(cm_rf[,1], cm_rf[,1] == 1)))) / 2
  
}


# ranking 23

wacc_IBk_1_23 <- c()
wacc_SMO_1_23 <- c()
wacc_Logistic_1_23 <- c()
wacc_JRip_1_23 <- c()
wacc_nb_1_23 <- c()
wacc_rf_1_23 <- c()

for(i in 1:results_1_23$cutoff$size[3]){
  
  temp <- rbind(input_1_2, input_1_3)
  train <- temp[, c(head(results_1_23$RI$attribute, i), "subtype_Survival..months.")]
  test <- input_1_1[, c(head(results_1_23$RI$attribute, i), "subtype_Survival..months.")]
  
  alg_IBk <- RWeka::IBk(subtype_Survival..months.~., train)   # knn
  testpred_IBk <- predict(alg_IBk, newdata = test)
  cm_IBk <- cbind(test$subtype_Survival..months., testpred_IBk)
  
  alg_SMO <- RWeka::SMO(subtype_Survival..months.~., train)
  testpred_SMO <- predict(alg_SMO, newdata = test)
  cm_SMO <- cbind(test$subtype_Survival..months., testpred_SMO)
  
  alg_Logistic <- RWeka::Logistic(subtype_Survival..months.~., train)
  testpred_Logistic <- predict(alg_Logistic, newdata = test)
  cm_Logistic <- cbind(test$subtype_Survival..months., testpred_Logistic)
  
  alg_JRip <- RWeka::JRip(subtype_Survival..months.~., train)
  testpred_JRip <- predict(alg_JRip, newdata = test)
  cm_JRip <- cbind(test$subtype_Survival..months., testpred_JRip)
  
  alg_nb <- e1071::naiveBayes(subtype_Survival..months.~., train)
  testpred_nb <- predict(alg_nb, newdata = test)
  cm_nb <- cbind(test$subtype_Survival..months., testpred_nb)
  
  alg_rf <- randomForest::randomForest(subtype_Survival..months.~., train)
  testpred_rf <- predict(alg_rf, newdata = test)
  cm_rf <- cbind(test$subtype_Survival..months., testpred_rf)
  
  
  wacc_IBk_1_23[i] <- ((length(subset(cm_IBk[,1], cm_IBk[,1] == 2 & cm_IBk[,2] == 2 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 2))) + 
                         (length(subset(cm_IBk[,1], cm_IBk[,1] == 1 & cm_IBk[,2] == 1 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 1)))) / 2
  
  wacc_SMO_1_23[i] <- ((length(subset(cm_SMO[,1], cm_SMO[,1] == 2 & cm_SMO[,2] == 2 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 2))) + 
                         (length(subset(cm_SMO[,1], cm_SMO[,1] == 1 & cm_SMO[,2] == 1 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 1)))) / 2
  
  wacc_Logistic_1_23[i] <- ((length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2 & cm_Logistic[,2] == 2 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2))) + 
                              (length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1 & cm_Logistic[,2] == 1 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1)))) / 2
  
  wacc_JRip_1_23[i] <- ((length(subset(cm_JRip[,1], cm_JRip[,1] == 2 & cm_JRip[,2] == 2 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 2))) + 
                          (length(subset(cm_JRip[,1], cm_JRip[,1] == 1 & cm_JRip[,2] == 1 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 1)))) / 2
  
  wacc_nb_1_23[i] <- ((length(subset(cm_nb[,1], cm_nb[,1] == 2 & cm_nb[,2] == 2 )) / length(subset(cm_nb[,1], cm_nb[,1] == 2))) + 
                        (length(subset(cm_nb[,1], cm_nb[,1] == 1 & cm_nb[,2] == 1 )) / length(subset(cm_nb[,1], cm_nb[,1] == 1)))) / 2
  
  wacc_rf_1_23[i] <- ((length(subset(cm_rf[,1], cm_rf[,1] == 2 & cm_rf[,2] == 2 )) / length(subset(cm_rf[,1], cm_rf[,1] == 2))) + 
                        (length(subset(cm_rf[,1], cm_rf[,1] == 1 & cm_rf[,2] == 1 )) / length(subset(cm_rf[,1], cm_rf[,1] == 1)))) / 2
  
}

wacc_IBk_1 <- c(wacc_IBk_1_12, wacc_IBk_1_13, wacc_IBk_1_23)
wacc_SMO_1 <- c(wacc_SMO_1_12, wacc_SMO_1_13, wacc_SMO_1_23)
wacc_Logistic_1 <- c(wacc_Logistic_1_12, wacc_Logistic_1_13, wacc_Logistic_1_23)
wacc_JRip_1 <- c(wacc_JRip_1_12, wacc_JRip_1_13, wacc_JRip_1_23)
wacc_nb_1 <- c(wacc_nb_1_12, wacc_nb_1_13, wacc_nb_1_23)
wacc_rf_1 <- c(wacc_rf_1_12, wacc_rf_1_13, wacc_rf_1_23)


# średnie dla poszczególnych algorytmów

mean(wacc_IBk_1)
mean(wacc_SMO_1) # najwyższa średnia (i mediana)
mean(wacc_Logistic_1)
mean(wacc_JRip_1)
mean(wacc_nb_1)
mean(wacc_rf_1)


# SMO 1

min_1_12 <- min(which(wacc_SMO_1_12 > 0.8)) 
min_1_13 <- min(which(wacc_SMO_1_13 > 0.8)) 
min_1_23 <- min(which(wacc_SMO_1_23 > 0.8)) 

logical_sum_1 <- unique(c(results_1_12$RI[1:2,]$attribute, results_1_13$RI[1:3,]$attribute, results_1_23$RI[1:3,]$attribute))

# saveRDS(logical_sum_1, "logical_sum_1.RDS")


# zapytanie do rentrez dla 7 genów (znaleziono jedynie 2 publikacje dla genu GUSBP2, dla pozostałych genów nie znaleziono żadnej publikacji)

pub_ids <- c()

for(i in 1:length(logical_sum_1)){
  search <- entrez_search(db="pubmed", term = ls1[i])
  pub_ids <- c(pub_ids, search$ids)
}

paper_links <- entrez_link(dbfrom="pubmed", id=27217703, cmd="llinks")
linkout_urls(paper_links)

paper_links_2 <- entrez_link(dbfrom="pubmed", id=23733509, cmd="llinks")
linkout_urls(paper_links_2)


### ANALIZA 2: budowa modelu przewidującego przeżywalność pacjentów

results_2_12 <- readRDS("results_2_12.RDS")
results_2_12$cutoff # algorytmy odcięcia

results_2_13 <- readRDS("results_2_13.RDS")
results_2_13$cutoff

results_2_23 <- readRDS("results_2_23.RDS")
results_2_23$cutoff


input_2_1 <- readRDS("input_2_1.RDS")
input_2_2 <- readRDS("input_2_2.RDS")
input_2_3 <- readRDS("input_2_3.RDS")

input_2_1 <- as.data.frame(unclass(input_2_1))
input_2_2 <- as.data.frame(unclass(input_2_2))
input_2_3 <- as.data.frame(unclass(input_2_3))


# labels vs wacc

# ranking 12

wacc_IBk_2_12 <- c()
wacc_SMO_2_12 <- c()
wacc_Logistic_2_12 <- c()
wacc_JRip_2_12 <- c()
wacc_nb_2_12 <- c()
wacc_rf_2_12 <- c()

for(i in 1:results_2_12$cutoff$size[3]){
  
  temp <- rbind(input_2_1, input_2_2)
  train <- temp[, c(head(results_2_12$RI$attribute, i), "subtype_Survival..months.")]
  test <- input_2_3[, c(head(results_2_12$RI$attribute, i), "subtype_Survival..months.")]
  
  alg_IBk <- RWeka::IBk(subtype_Survival..months.~., train)  
  testpred_IBk <- predict(alg_IBk, newdata = test)
  cm_IBk <- cbind(test$subtype_Survival..months., testpred_IBk)
  
  alg_SMO <- RWeka::SMO(subtype_Survival..months.~., train)
  testpred_SMO <- predict(alg_SMO, newdata = test)
  cm_SMO <- cbind(test$subtype_Survival..months., testpred_SMO)
  
  alg_Logistic <- RWeka::Logistic(subtype_Survival..months.~., train)
  testpred_Logistic <- predict(alg_Logistic, newdata = test)
  cm_Logistic <- cbind(test$subtype_Survival..months., testpred_Logistic)
  
  alg_JRip <- RWeka::JRip(subtype_Survival..months.~., train)
  testpred_JRip <- predict(alg_JRip, newdata = test)
  cm_JRip <- cbind(test$subtype_Survival..months., testpred_JRip)
  
  alg_nb <- e1071::naiveBayes(subtype_Survival..months.~., train)
  testpred_nb <- predict(alg_nb, newdata = test)
  cm_nb <- cbind(test$subtype_Survival..months., testpred_nb)
  
  alg_rf <- randomForest::randomForest(subtype_Survival..months.~., train)
  testpred_rf <- predict(alg_rf, newdata = test)
  cm_rf <- cbind(test$subtype_Survival..months., testpred_rf)
  
  
  wacc_IBk_2_12[i] <- ((length(subset(cm_IBk[,1], cm_IBk[,1] == 2 & cm_IBk[,2] == 2 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 2))) + 
                         (length(subset(cm_IBk[,1], cm_IBk[,1] == 1 & cm_IBk[,2] == 1 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 1)))) / 2
  
  wacc_SMO_2_12[i] <- ((length(subset(cm_SMO[,1], cm_SMO[,1] == 2 & cm_SMO[,2] == 2 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 2))) + 
                         (length(subset(cm_SMO[,1], cm_SMO[,1] == 1 & cm_SMO[,2] == 1 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 1)))) / 2
  
  wacc_Logistic_2_12[i] <- ((length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2 & cm_Logistic[,2] == 2 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2))) + 
                              (length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1 & cm_Logistic[,2] == 1 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1)))) / 2
  
  wacc_JRip_2_12[i] <- ((length(subset(cm_JRip[,1], cm_JRip[,1] == 2 & cm_JRip[,2] == 2 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 2))) + 
                          (length(subset(cm_JRip[,1], cm_JRip[,1] == 1 & cm_JRip[,2] == 1 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 1)))) / 2
  
  wacc_nb_2_12[i] <- ((length(subset(cm_nb[,1], cm_nb[,1] == 2 & cm_nb[,2] == 2 )) / length(subset(cm_nb[,1], cm_nb[,1] == 2))) + 
                        (length(subset(cm_nb[,1], cm_nb[,1] == 1 & cm_nb[,2] == 1 )) / length(subset(cm_nb[,1], cm_nb[,1] == 1)))) / 2
  
  wacc_rf_2_12[i] <- ((length(subset(cm_rf[,1], cm_rf[,1] == 2 & cm_rf[,2] == 2 )) / length(subset(cm_rf[,1], cm_rf[,1] == 2))) + 
                        (length(subset(cm_rf[,1], cm_rf[,1] == 1 & cm_rf[,2] == 1 )) / length(subset(cm_rf[,1], cm_rf[,1] == 1)))) / 2
  
}


# ranking 13

wacc_IBk_2_13 <- c()
wacc_SMO_2_13 <- c()
wacc_Logistic_2_13 <- c()
wacc_JRip_2_13 <- c()
wacc_nb_2_13 <- c()
wacc_rf_2_13 <- c()

for(i in 1:results_2_13$cutoff$size[3]){
  
  temp <- rbind(input_2_1, input_2_3)
  train <- temp[, c(head(results_2_13$RI$attribute, i), "subtype_Survival..months.")]
  test <- input_2_2[, c(head(results_2_13$RI$attribute, i), "subtype_Survival..months.")]
  
  alg_IBk <- RWeka::IBk(subtype_Survival..months.~., train)   # knn
  testpred_IBk <- predict(alg_IBk, newdata = test)
  cm_IBk <- cbind(test$subtype_Survival..months., testpred_IBk)
  
  alg_SMO <- RWeka::SMO(subtype_Survival..months.~., train)
  testpred_SMO <- predict(alg_SMO, newdata = test)
  cm_SMO <- cbind(test$subtype_Survival..months., testpred_SMO)
  
  alg_Logistic <- RWeka::Logistic(subtype_Survival..months.~., train)
  testpred_Logistic <- predict(alg_Logistic, newdata = test)
  cm_Logistic <- cbind(test$subtype_Survival..months., testpred_Logistic)
  
  alg_JRip <- RWeka::JRip(subtype_Survival..months.~., train)
  testpred_JRip <- predict(alg_JRip, newdata = test)
  cm_JRip <- cbind(test$subtype_Survival..months., testpred_JRip)
  
  alg_nb <- e1071::naiveBayes(subtype_Survival..months.~., train)
  testpred_nb <- predict(alg_nb, newdata = test)
  cm_nb <- cbind(test$subtype_Survival..months., testpred_nb)
  
  alg_rf <- randomForest::randomForest(subtype_Survival..months.~., train)
  testpred_rf <- predict(alg_rf, newdata = test)
  cm_rf <- cbind(test$subtype_Survival..months., testpred_rf)
  
  
  wacc_IBk_2_13[i] <- ((length(subset(cm_IBk[,1], cm_IBk[,1] == 2 & cm_IBk[,2] == 2 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 2))) + 
                         (length(subset(cm_IBk[,1], cm_IBk[,1] == 1 & cm_IBk[,2] == 1 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 1)))) / 2
  
  wacc_SMO_2_13[i] <- ((length(subset(cm_SMO[,1], cm_SMO[,1] == 2 & cm_SMO[,2] == 2 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 2))) + 
                         (length(subset(cm_SMO[,1], cm_SMO[,1] == 1 & cm_SMO[,2] == 1 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 1)))) / 2
  
  wacc_Logistic_2_13[i] <- ((length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2 & cm_Logistic[,2] == 2 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2))) + 
                              (length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1 & cm_Logistic[,2] == 1 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1)))) / 2
  
  wacc_JRip_2_13[i] <- ((length(subset(cm_JRip[,1], cm_JRip[,1] == 2 & cm_JRip[,2] == 2 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 2))) + 
                          (length(subset(cm_JRip[,1], cm_JRip[,1] == 1 & cm_JRip[,2] == 1 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 1)))) / 2
  
  wacc_nb_2_13[i] <- ((length(subset(cm_nb[,1], cm_nb[,1] == 2 & cm_nb[,2] == 2 )) / length(subset(cm_nb[,1], cm_nb[,1] == 2))) + 
                        (length(subset(cm_nb[,1], cm_nb[,1] == 1 & cm_nb[,2] == 1 )) / length(subset(cm_nb[,1], cm_nb[,1] == 1)))) / 2
  
  wacc_rf_2_13[i] <- ((length(subset(cm_rf[,1], cm_rf[,1] == 2 & cm_rf[,2] == 2 )) / length(subset(cm_rf[,1], cm_rf[,1] == 2))) + 
                        (length(subset(cm_rf[,1], cm_rf[,1] == 1 & cm_rf[,2] == 1 )) / length(subset(cm_rf[,1], cm_rf[,1] == 1)))) / 2
  
}


# ranking 23

wacc_IBk_2_23 <- c()
wacc_SMO_2_23 <- c()
wacc_Logistic_2_23 <- c()
wacc_JRip_2_23 <- c()
wacc_nb_2_23 <- c()
wacc_rf_2_23 <- c()

for(i in 1:results_2_23$cutoff$size[3]){
  
  temp <- rbind(input_2_2, input_2_3)
  train <- temp[, c(head(results_2_23$RI$attribute, i), "subtype_Survival..months.")]
  test <- input_2_1[, c(head(results_2_23$RI$attribute, i), "subtype_Survival..months.")]
  
  alg_IBk <- RWeka::IBk(subtype_Survival..months.~., train)   # knn
  testpred_IBk <- predict(alg_IBk, newdata = test)
  cm_IBk <- cbind(test$subtype_Survival..months., testpred_IBk)
  
  alg_SMO <- RWeka::SMO(subtype_Survival..months.~., train)
  testpred_SMO <- predict(alg_SMO, newdata = test)
  cm_SMO <- cbind(test$subtype_Survival..months., testpred_SMO)
  
  alg_Logistic <- RWeka::Logistic(subtype_Survival..months.~., train)
  testpred_Logistic <- predict(alg_Logistic, newdata = test)
  cm_Logistic <- cbind(test$subtype_Survival..months., testpred_Logistic)
  
  alg_JRip <- RWeka::JRip(subtype_Survival..months.~., train)
  testpred_JRip <- predict(alg_JRip, newdata = test)
  cm_JRip <- cbind(test$subtype_Survival..months., testpred_JRip)
  
  alg_nb <- e1071::naiveBayes(subtype_Survival..months.~., train)
  testpred_nb <- predict(alg_nb, newdata = test)
  cm_nb <- cbind(test$subtype_Survival..months., testpred_nb)
  
  alg_rf <- randomForest::randomForest(subtype_Survival..months.~., train)
  testpred_rf <- predict(alg_rf, newdata = test)
  cm_rf <- cbind(test$subtype_Survival..months., testpred_rf)
  
  
  wacc_IBk_2_23[i] <- ((length(subset(cm_IBk[,1], cm_IBk[,1] == 2 & cm_IBk[,2] == 2 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 2))) + 
                         (length(subset(cm_IBk[,1], cm_IBk[,1] == 1 & cm_IBk[,2] == 1 )) / length(subset(cm_IBk[,1], cm_IBk[,1] == 1)))) / 2
  
  wacc_SMO_2_23[i] <- ((length(subset(cm_SMO[,1], cm_SMO[,1] == 2 & cm_SMO[,2] == 2 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 2))) + 
                         (length(subset(cm_SMO[,1], cm_SMO[,1] == 1 & cm_SMO[,2] == 1 )) / length(subset(cm_SMO[,1], cm_SMO[,1] == 1)))) / 2
  
  wacc_Logistic_2_23[i] <- ((length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2 & cm_Logistic[,2] == 2 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 2))) + 
                              (length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1 & cm_Logistic[,2] == 1 )) / length(subset(cm_Logistic[,1], cm_Logistic[,1] == 1)))) / 2
  
  wacc_JRip_2_23[i] <- ((length(subset(cm_JRip[,1], cm_JRip[,1] == 2 & cm_JRip[,2] == 2 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 2))) + 
                          (length(subset(cm_JRip[,1], cm_JRip[,1] == 1 & cm_JRip[,2] == 1 )) / length(subset(cm_JRip[,1], cm_JRip[,1] == 1)))) / 2
  
  wacc_nb_2_23[i] <- ((length(subset(cm_nb[,1], cm_nb[,1] == 2 & cm_nb[,2] == 2 )) / length(subset(cm_nb[,1], cm_nb[,1] == 2))) + 
                        (length(subset(cm_nb[,1], cm_nb[,1] == 1 & cm_nb[,2] == 1 )) / length(subset(cm_nb[,1], cm_nb[,1] == 1)))) / 2
  
  wacc_rf_2_23[i] <- ((length(subset(cm_rf[,1], cm_rf[,1] == 2 & cm_rf[,2] == 2 )) / length(subset(cm_rf[,1], cm_rf[,1] == 2))) + 
                        (length(subset(cm_rf[,1], cm_rf[,1] == 1 & cm_rf[,2] == 1 )) / length(subset(cm_rf[,1], cm_rf[,1] == 1)))) / 2
  
}

wacc_IBk_2 <- c(wacc_IBk_2_12, wacc_IBk_2_13, wacc_IBk_2_23)
wacc_SMO_2 <- c(wacc_SMO_2_12, wacc_SMO_2_13, wacc_SMO_2_23)
wacc_Logistic_2 <- c(wacc_Logistic_2_12, wacc_Logistic_2_13, wacc_Logistic_2_23)
wacc_JRip_2 <- c(wacc_JRip_2_12, wacc_JRip_2_13, wacc_JRip_2_23)
wacc_nb_2 <- c(wacc_nb_2_12, wacc_nb_2_13, wacc_nb_2_23)
wacc_rf_2 <- c(wacc_rf_2_12, wacc_rf_2_13, wacc_rf_2_23)


# średnie dla poszczególnych algorytmów

mean(wacc_IBk_2)
mean(wacc_SMO_2) # najwyższa średnia 
mean(wacc_Logistic_2)
mean(wacc_JRip_2)
mean(wacc_nb_2)
mean(wacc_rf_2)

min_2_12 <- min(which(wacc_SMO_2_12 > 0.8))
min_2_13 <- min(which(wacc_SMO_2_13 > 0.8))
min_2_23 <- min(which(wacc_SMO_2_23 > 0.8))

logical_sum_2 <- unique(c(results_2_12$RI[1,]$attribute, results_2_13$RI[1,]$attribute, results_2_23$RI[1,]$attribute)) # tylko 1 atrybut: disease_type


### porównanie rankingów pomiędzy analizą 1 i analizą 2

rank_1 <- unique(c(results_1_12$RI$attribute[1:results_1_12$cutoff$size[3]], results_1_13$RI$attribute[1:results_1_13$cutoff$size[3]], 
                   results_1_23$RI$attribute[1:results_1_23$cutoff$size[3]]))

rank_2 <- unique(c(results_2_12$RI$attribute[1:results_2_12$cutoff$size[3]], results_2_13$RI$attribute[1:results_2_13$cutoff$size[3]], 
                   results_2_23$RI$attribute[1:results_2_23$cutoff$size[3]]))  

rank_2 <- rank_2[c(-1, -62, -76, -82)] # odrzucenie cech klinicznych z rankingu dla analizy 2

table(rank_2 %in% rank_1) 


### ANALIZA 1: dalsze obliczenia

m <- (results_1_12$cutoff$size[3] + results_1_13$cutoff$size[3] + results_1_23$cutoff$size[3]) / 3

cross_genes_1 <- results_1_23$RI[1:m,]$attribute[(results_1_23$RI[1:m,]$attribute %in% results_1_13$RI[1:m,]$attribute) == TRUE 
                                                 & (results_1_23$RI[1:m,]$attribute %in% results_1_12$RI[1:m,]$attribute) == TRUE]

# saveRDS(cross_genes_1, "cross_genes_1.RDS")


### zamiana formatu genów z ensemble na "zwykły"

# logical_sum_1 <- readRDS("logical_sum_1.RDS")
# cross_genes_1 <- readRDS("cross_genes_1.RDS")

genes <- readRDS("gene_notations.rds")

ls1 <- sub("\\.1", "", logical_sum_1 )
c1 <- sub("\\.1", "", cross_genes_1)

ls1 <- unique(subset(genes, genes$ensemble %in% ls1))
ls1 <- ls1$name

genes1 <- unique(subset(genes, genes$ensemble %in% c1))

for(i in 1:length(genes1$ensemble)) {
  if(genes1$ensemble[i] == "ENSG00000255963"){
    genes1$ensemble[i] = "ENSG00000255963.1"
  }
  if(genes1$ensemble[i] == "ENSG00000198161"){
    genes1$ensemble[i] = "ENSG00000198161.1"
  }
  if(genes1$ensemble[i] == "ENSG00000255854"){
    genes1$ensemble[i] = "ENSG00000255854.1"
  }
  if(genes1$ensemble[i] == "ENSG00000203832"){
    genes1$ensemble[i] = "ENSG00000203832.1"
  }
}


### heatmap (tu generuje się rysunek 4)

heatmap_table <- target[ , c1]
names(heatmap_table) <- genes1$name
tail(sort(colMeans(heatmap_table)))

sort(target$subtype_Survival..months.)
target$subtype_Survival..months.[1:105] <- "S"
target$subtype_Survival..months.[106:203] <- "L"

heatmap_table <- cbind(target$subtype_Survival..months., heatmap_table)

heatmap_table_S <- subset(heatmap_table, target$subtype_Survival..months. == "S")
heatmap_table_L <- subset(heatmap_table, target$subtype_Survival..months. == "L")

median_S <- sapply(heatmap_table_S[, -1], median)
median_L <- sapply(heatmap_table_L[, -1], median)

mediana <- cbind(median_S, median_L)

median_diff <- sort(abs(median_L - median_S))

library(plotly)

p <- plot_ly(z = mediana, type = "heatmap", colors = colorRamp(c("yellow", "red")), y = c1, x = c("S", "L"))

# export(p, out_file = "rysunek_4.png") # zapisywanie do pliku nie działa, plik z rysunkiem 4 został stworzony manualnie, na podstawie zmiennej 'p'


### geneSummary (tu generują się rysunki 5, 6, 7)

library("GeneSummary")
library("ggplot2")
library("mygene")

# ingenes <- readRDS("cross_genes_1.RDS")

ingenes <- cross_genes_1

mylabel <- 'my_genes'

geneData <- loadGeneData()

gene_summary <- getGeneSummary(ingenes, geneData)

# saveRDS(gene_summary, "gs.RDS")
# gene_summary <- readRDS("gs.RDS")

cnt_items_result <- countListItems(gene_summary$gene_description, mylistElement='gene_ontology_BP', min_plot_freq = 2)

pdf("rysunek_4.pdf")
plot(cnt_items_result$bio_function_plot)
dev.off()


# wordcloud

gdoc <- paste(unlist(lapply(gene_summary$gene_description, wrapDescription)), collapse = " ")

pdf("rysunek_6.pdf")
word_cloud <- getWordcloud(gdoc, type="text", lang="english", excludeWords=c(bioStopWords()), textStemming=F, 
                           colorPalette="Dark2", min.freq=1, max.words=200, wc = 1)
dev.off()


# chromosomy

gene_summary_plots <- getGeneSummaryPlots(gene_summary, geneData, coding_genes = T)

pdf("rysunek_7.pdf")
gene_summary_plots$chromosome_cnt_plot
dev.off()


### zapytanie do rentrez (dla 53 genów)

genes_names <- genes1$name

counts <- c()

library(rentrez)

for(i in 1:length(genes_names)){
  query <- paste(genes_names[i]," AND (Glioma OR Cancer OR Glioblastoma)") 
  search <- entrez_search(db="pubmed", term = query)
  counts <- c(counts, search$count)
}

names_counts <- data.frame(genes = genes_names, frequency = counts)

names_counts <- names_counts[order(names_counts$frequency),]

names_counts_sorted <- data.frame(genes = names_counts$genes, frequency = names_counts$frequency)

# saveRDS(names_counts_sorted, "names_counts_sorted.RDS")

pdf(file = "rysunek_8.pdf")
par(mar=c(3,7,0,3))
barplot(names_counts_sorted$frequency, horiz = TRUE, names.arg = names_counts_sorted$genes, 
        las = 1, cex.names=0.75, xlim = c(0,500))
dev.off()


### utworzenie zapytania do narządzia X2K

genes_1_string <- c()

for(i in 1:length(genes_names)){
  genes_1_string <- cat(paste(genes_1_string, genes_names[i], sep = '\n'))
}


### INTERAKCJE (rmcfs), wartość 25 ustawiona arbitralnie (tu generuje się rysunek 9)

gid <- build.idgraph(results_1_12)

pdf(file = "rysunek_9.pdf")
plot(gid, label_dist = 1, cex = 0.6)
dev.off()


ID_12 <- results_1_12$ID[1:25,]
ID_13 <- results_1_13$ID[1:25,]
ID_23 <- results_1_23$ID[1:25,]

sum_ID <- rbind(results_1_12$ID[1:25,], results_1_13$ID[1:25,], results_1_23$ID[1:25,])


# sprawdzenie, czy pomiędzy trzema zbiorami występują wspólne interakcje

sum_ID$unique <- paste0(sum_ID$edge_a, sum_ID$edge_b) # paste0 skleja bez separatora
length(unique(sum_ID$unique)) # jest 75, więc nie ma pomiędzy podzbiorami wspólnych interakcji
sum_ID <- subset(sum_ID, select = - unique)


# odrzucenie kropki i tego, co po niej występuje

for(i in 10:1){
  sum_ID$edge_a <- sub(paste0("\\.", i), "", sum_ID$edge_a)
}

for(i in 10:1){
  sum_ID$edge_b <- sub(paste0("\\.", i), "", sum_ID$edge_b)
}

# zamiana formatu genów z ensemble na "zwykły"

for(i in 1:length(sum_ID$edge_a)){
  sum_ID$edge_a[i] <- unique(subset(genes, genes$ensemble %in% sum_ID$edge_a[i]))$name
}

for(i in 1:length(sum_ID$edge_b)){
  sum_ID$edge_b[i] <- unique(subset(genes, genes$ensemble %in% sum_ID$edge_b[i]))$name
}

sum_ID <- sum_ID[with(sum_ID,order(weight, decreasing = TRUE)), ]

geny_id <- unique(c(sum_ID$edge_a, sum_ID$edge_b)) 

# wektor 'entrez' został stworzony manualnie, na podstawie wektora 'geny_id', przy użyciu NCBI

entrez <- c(618, 667, 2103, 2502, 2653, 2678, 3167, 3213, 3350, 3563, 104472715, 3885, 5273, 6139, 100526842, 6155, 7319, 7389, 7455, 7643,
            7756, 8365, 8366, 8364, 8294, 8361, 8359, 554313, 8346, 8356, 9168, 9482, 10151, 10276, 10352, 11321, 11325, 22852, 26263, 27332,
            28996, 29097, 29117, 54576, 55073, 58505, 63946, 80133, 80834, 83897, 83955, 727851, 92086, 116093, 128178, 130106, 140690, 146822, 147409, 159296,
            163747, 221710, 245932, 245973, 246184, 259230, 338339, 339501, 342538, 386757, 387036, 389170, 390036, 390144, 390261, 392391, 403315, 440093, 441251, 619418,
            441818, 548645, 552889, 552891, 606551, 642946, 643905, 645455, 645683, 647033, 653505, 653598, 677796, 642659, 728773, 100151683)

# wektor 'interaction' został stworzony manualnie, przy pomocy narzędzia HumanNet

interaction <- c(3, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 2, 5, 3, 3, 5, 2, 5, 3, 5, 3, 5, 3,
                 5, 3, 3, 5, 3, 3, 2, 3, 5, 5, 5, 5, 1, 3, 3, 0, 2, 5, 5, 5, 5, 5, 5, 3, 5,
                 5, 5, 3, 3, 2, 3, 5, 5, 3, 3, 2, 5, 5, 2, 3, 5, 3, 5, 5, 5, 5, 2, 5, 5, 5) 

# 0 - bezpośrednia, 1 - pośrednia przez geny z inputu, 2 - interakcja nie występuje, 3 - disconnected, 4 - not found, 5 - invalid


### ZAŁĄCZNIKI

# tu generuje się kod tex do tabel 5.1 - 5.6

t_1_12 <- data.frame(results_1_12$RI$attribute[1:results_1_12$cutoff$size[3]], results_1_12$RI$RI_norm[1:results_1_12$cutoff$size[3]])
names(t_1_12) <- c("atrybut", "wartość RI")

t_1_13 <- data.frame(results_1_13$RI$attribute[1:results_1_13$cutoff$size[3]], results_1_13$RI$RI_norm[1:results_1_13$cutoff$size[3]])
names(t_1_13) <- c("atrybut", "wartość RI")

t_1_23 <- data.frame(results_1_23$RI$attribute[1:results_1_23$cutoff$size[3]], results_1_23$RI$RI_norm[1:results_1_23$cutoff$size[3]])
names(t_1_23) <- c("atrybut", "wartość RI")

t_2_12 <- data.frame(results_2_12$RI$attribute[1:results_2_12$cutoff$size[3]], results_2_12$RI$RI_norm[1:results_2_12$cutoff$size[3]])
names(t_2_12) <- c("atrybut", "wartość RI")

t_2_13 <- data.frame(results_2_13$RI$attribute[1:results_2_13$cutoff$size[3]], results_2_13$RI$RI_norm[1:results_2_13$cutoff$size[3]])
names(t_2_13) <- c("atrybut", "wartość RI")

t_2_23 <- data.frame(results_2_23$RI$attribute[1:results_2_23$cutoff$size[3]], results_2_23$RI$RI_norm[1:results_2_23$cutoff$size[3]])
names(t_2_23) <- c("atrybut", "wartość RI")

library(xtable)

xtable(t_1_12)
xtable(t_1_13)
xtable(t_1_23)

xtable(t_2_12)
xtable(t_2_13)
xtable(t_2_23)


# tu generuje się kod tex do tabeli 5.7 (potrzebny dataframe 'genes1')

# wektor 'description' został zrobiony manualnie, przy użyciu NCBI

description <- c("Dystonin", "Dynein Cytoplasmic", "Ferritin Heavy Chain 1 Pseudogene 10", "Glycine Cleavage System Protein H", "Ik Cytokine", 
                 "Lactate Dehydrogenase A", "Peptidylprolyl Isomerase A", "Rap1b, Member Of Ras Oncogene Family", "Ribosomal Protein L12", "Ribosomal Protein L17",
                 "Rpl17-C18orf32 Readthrough", "Ubiquitin Conjugating Enzyme E2 A", "Zinc Finger Protein 33a", "Zinc Finger Protein 90", 
                 "Hect And Rld Domain Containing E3 Ubiquitin Protein Ligase 2","Hect And Rld Domain Containing E3 Ubiquitin Protein Ligase Family Member 1", 
                 "Heterogeneous Nuclear Ribonucleoprotein A3 Pseudogene 1", "Poly(A) Binding Protein Interacting Protein 1", "Ankyrin Repeat Domain 26", "Nbpf Member 20",
                 "Methylmalonic Aciduria And Homocystinuria, Cbld Type", "Zinc Finger Protein 638", "Leucine Rich Repeat Containing 37 Member A4, Pseudogene", 
                 "Trinucleotide Repeat Containing 6c", "Oligosaccharyltransferase Complex Non-Catalytic Subunit", "Mitochondrial Ribosomal Protein L36", "Naca Family Member 4, Pseudogene", 
                 "Rab6c, Member Ras Oncogene Family", "Threonyl-Trna Synthetase Like 2", "Edar Associated Death Domain",
                 "Ccctc-Binding Factor Like", "C-X9-C Motif Containing 1", "Atpase Phospholipid Transporting 8b5, Pseudogene", "Family With Sequence Similarity 71 Member D",
                 "Ribosomal Protein L23 Pseudogene 8", "Atpase H+ Transporting V1 Subunit C2", "Cell Division Cycle 26", "High Mobility Group Nucleosomal Binding Domain 2 Pseudogene 46",
                 "Pre-Mrna Processing Factor 31", "Nascent Polypeptide Associated Complex Subunit Alpha 2",
                 "Solute Carrier Family 6 Member 10, Pseudogene", "Glucuronidase, Beta Pseudogene 2", "H3 Histone Family Member 3c", "Family With Sequence Similarity 90 Member A21, Pseudogene", 
                 "Ww Domain Binding Protein 11 Pseudogene 1", "Ubiquitin Conjugating Enzyme E2 M Pseudogene 1", "Centrosomal Protein 170 Pseudogene 1", "Ribosomal Protein L13a Pseudogene 3",
                 "Proliferation-Associated 2g4 Pseudogene 4", "Peptidylprolyl Isomerase A Like 4a",
                 "Peptidylprolyl Isomerase A Like 4a", "Peptidylprolyl Isomerase A Like 4c", "Heterogeneous Nuclear Ribonucleoprotein A1 Pseudogene 48")


for(i in 1:length(genes1$ensemble)) {
  temp_12 <- subset(results_1_12$RI, attribute == c1[i])
  temp_13 <- subset(results_1_13$RI, attribute == c1[i])
  temp_23 <- subset(results_1_23$RI, attribute == c1[i])
  genes1$RI[i] <- (temp_12$RI_norm + temp_13$RI_norm + temp_23$RI_norm) / 3
}

geny_53 <- data.frame(genes1$name, description, genes1$RI)
names(geny_53) <- c("nazwa", "opis (wg NCBI)", "wartość RI")

library(xtable)
xtable(geny_53)


# tu generuje się kod tex do tabeli 5.8 (potrzebny dataframe 'sum_ID')

sum_ID$interaction <- interaction

interakcje <- data.frame(sum_ID$edge_a, sum_ID$edge_b, sum_ID$weight, sum_ID$interaction)

names(interakcje) <- c("gen1", "gen2", "wartość", "interakcja")

for (i in 1:length(interaction)) {
  if (interaction[i] == 0) {
    interaction[i] = "bezpośrednia"
  }
  else if (interaction[i] == 1) {
    interaction[i] = "pośrednia"
  }  
  else {
    interaction[i] = "brak"
  }
}

interakcje$interakcja <- as.factor(interaction)

library(xtable)
xtable(interakcje)




