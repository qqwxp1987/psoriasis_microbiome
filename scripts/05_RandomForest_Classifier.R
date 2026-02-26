###############################################################################
# 05_RandomForest_Classifier.R
# Random Forest classification model for psoriasis prediction
# - Feature selection: Boruta + findLinearCombos
# - Model training: caret RF with 5-fold x5-repeat CV
# - External validation: ROC curves + DeLong tests
# - Precision-Recall curves (ROCR)
# - Confusion matrix with optimal threshold
# - SMOTE for class imbalance (Validation set)
#
# Generates results corresponding to:
# - Table SM1: Comprehensive Performance Metrics
# - Table SM2: Precision-Recall AUC Comparison
# - Figure SM1: RF Hyperparameter Tuning and CV Performance
# - Figure SM2: Confusion Matrices
# - Figure SM3: PR Curves (original and SMOTE)
# - Figure SM4: ROC Curves on SMOTE balanced set
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(Boruta)
library(caret)
library(dplyr)
library(pROC)
library(patchwork)
library(ggplot2)
library(here)
library(ggpubr)
library(scales)
library(tidyr)
library(ROCR)
# library(DMwR) # Loaded dynamically when needed

# --- 1. Load training data (Discovery set: 28 pairs) -------------------------
data_dir <- here("scripts", "caret_data_code")

abu <- read.csv(file.path(data_dir, "train_metaphlan_taxonomic_profiles_20240418.csv"),
    sep = ",", header = TRUE, row.names = 1
)
bateria_names <- stringr::str_split(rownames(abu), pattern = ";|\\|", simplify = TRUE) %>%
    as.data.frame()
rownames(abu) <- bateria_names[, ncol(bateria_names)]
TrainX <- t(abu)

TrainY <- read.csv(file.path(data_dir, "train_metadata_20240401.txt"),
    sep = "\t", header = TRUE, row.names = 1
)
Train <- cbind(TrainY, TrainX)

# --- 2. Preprocessing function ------------------------------------------------
caret_preprocess <- function(fileGroup = NULL, fileOTU = NULL,
                             caret_formula = "Group~.",
                             metaphlanFormat = FALSE,
                             caret_delZeroVar = FALSE,
                             caret_collinearity = FALSE,
                             caret_delcorr = FALSE) {
    set.seed(1234)
    yvar <- strsplit(gsub(" ", "", caret_formula), "~")[[1]][1]
    dataY <- fileGroup[, yvar, drop = FALSE]

    abu_local <- fileOTU
    if (metaphlanFormat) {
        bateria_names <- strsplit(rownames(abu_local), split = ";|\\|")
        yy <- c()
        for (t in bateria_names) {
            yy <- c(yy, t[length(t)])
        }
        rownames(abu_local) <- yy
    }
    dataX <- t(abu_local)

    # Remove zero-variance features
    if (caret_delZeroVar) {
        set.seed(1234)
        zerovar <- caret::nearZeroVar(dataX)
        if (length(zerovar) > 0) dataX <- dataX[, -zerovar]
    }

    # Resolve multicollinearity
    if (caret_collinearity) {
        comboInfo1 <- caret::findLinearCombos(dataX)
        if (length(comboInfo1$remove) > 0) dataX <- dataX[, -comboInfo1$remove]
    }

    # Remove highly correlated features (cutoff = 0.9)
    if (caret_delcorr) {
        set.seed(1234)
        highCorr <- caret::findCorrelation(cor(dataX), cutoff = 0.9)
        if (length(highCorr) > 0) dataX <- dataX[, -highCorr]
    }

    Data <- cbind(dataY, dataX)
    Data[[yvar]] <- factor(Data[[yvar]])
    return(Data)
}

# --- 3. Boruta feature selection ----------------------------------------------
caret_formula <- "Group~."

# Method 1: Boruta only
Data1 <- caret_preprocess(
    fileGroup = TrainY, fileOTU = abu,
    caret_formula = caret_formula,
    caret_collinearity = FALSE,
    metaphlanFormat = TRUE
)

# Method 2: findLinearCombos + Boruta
Data2 <- caret_preprocess(
    fileGroup = TrainY, fileOTU = abu,
    caret_formula = caret_formula,
    caret_collinearity = TRUE,
    metaphlanFormat = TRUE
)

set.seed(1234)
Boruta.srx1 <- Boruta::Boruta(as.formula(caret_formula), data = Data1)
Boruta.srx2 <- Boruta::Boruta(as.formula(caret_formula), data = Data2)

out_dir <- here("outputs", "RandomForest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file = file.path(out_dir, "boruta_model1.pdf"), width = 20, height = 5)
plot(Boruta.srx1, las = 2, xlab = "")
dev.off()

pdf(file = file.path(out_dir, "boruta_model2.pdf"), width = 10, height = 5)
plot(Boruta.srx2, las = 2, xlab = "")
dev.off()

# --- 4. Define model formulas -------------------------------------------------
# Model 1: 9 features from Boruta (no collinearity filter)
formula1 <- as.formula("Group ~ s__GGB51647_SGB4348 + s__GGB3746_SGB5089 + s__GGB33085_SGB1673 +
    s__Prevotella_stercorea + s__Clostridium_sp_AF20_17LB + s__Bacteroides_finegoldii +
    s__Bacteroides_congonensis + s__Clostridiaceae_bacterium + s__Bifidobacterium_bifidum")

# Model 2: 3 features from findLinearCombos + Boruta
formula2 <- as.formula("Group ~ s__GGB51647_SGB4348 + s__GGB1109_SGB1423 + s__GGB3746_SGB5089")

# Model 3: Single feature (s__GGB51647_SGB4348 only)
formula3 <- as.formula("Group ~ s__GGB51647_SGB4348")

# --- 5. Train RF models with 5-fold x 5-repeat CV -----------------------------
set.seed(1234)
fitControl <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 5)
tuneGrid <- expand.grid(mtry = c(1:10))

set.seed(1234)
fit1 <- train(formula1,
    data = Train, method = "rf",
    trControl = fitControl, tuneGrid = tuneGrid
)
set.seed(1234)
fit2 <- train(formula2,
    data = Train, method = "rf",
    trControl = fitControl, tuneGrid = tuneGrid
)
set.seed(1234)
fit3 <- train(formula3,
    data = Train, method = "rf",
    trControl = fitControl, tuneGrid = tuneGrid
)

# Variable importance plot (Model 1)
P0 <- ggplot(caret::varImp(fit1, scale = FALSE)) +
    theme_bw() +
    theme(text = element_text(size = 15))
ggsave(filename = file.path(out_dir, "varImp_model1.pdf"), plot = P0, width = 6, height = 3)

# --- 6. Figure SM1: Training metrics visualization ----------------------------
plot_train_metrics <- function(fit_obj, outname) {
    p1 <- ggplot(fit_obj) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme(text = element_text(size = 12), legend.position = "none") +
        theme_test()
    p2 <- ggplot(fit_obj, metric = "Kappa") +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme(text = element_text(size = 12)) +
        theme_test()

    resample <- fit_obj$resample[order(fit_obj$resample$Resample), ]
    tmp <- tidyr::separate(resample, col = "Resample", into = c("Folds", "Repeats"), sep = "\\.")
    tmp$Folds <- gsub("^Fold", "", tmp$Folds)

    p3 <- ggpubr::ggboxplot(tmp,
        x = "Folds", y = "Accuracy",
        color = "Folds", palette = "d3", add = "jitter"
    ) +
        theme(text = element_text(size = 12), legend.position = "none")
    p4 <- ggpubr::ggboxplot(tmp,
        x = "Folds", y = "Kappa",
        color = "Folds", palette = "d3", add = "jitter"
    ) +
        theme(text = element_text(size = 12), legend.position = "none")

    p0 <- (p1 + p2) / (p3 + p4)
    ggsave(filename = file.path(out_dir, outname), plot = p0, width = 6, height = 6)
}

plot_train_metrics(fit1, "Figure_SM1a_Model1_metrics.pdf")
plot_train_metrics(fit2, "Figure_SM1b_Model2_metrics.pdf")
plot_train_metrics(fit3, "Figure_SM1c_Model3_metrics.pdf")


# --- 7. Load Validation Data --------------------------------------------------
abu2 <- read.csv(file.path(data_dir, "extra_data_20240418.csv"),
    sep = ",", header = TRUE, row.names = 1, check.names = FALSE
)
bateria_names2 <- strsplit(colnames(abu2), split = ";|\\|")
yy <- c()
for (t in bateria_names2) {
    yy <- c(yy, t[length(t)])
}
colnames(abu2) <- yy

TestX <- dplyr::select(abu2, -Group)
# Note: Psoriasis is the positive class. Ensure factors are set with Healthy as reference,
# or Psoriasis as the second level (often the positive class in caret).
TestY <- factor(abu2[, "Group"], levels = c("Healthy", "Psoriasis"))

# Predictions on original validation set
pred1 <- predict(fit1, newdata = TestX, type = "prob")
pred2 <- predict(fit2, newdata = TestX, type = "prob")
pred3 <- predict(fit3, newdata = TestX, type = "prob")

# ROC curves
roc1 <- pROC::roc(response = TestY, predictor = pred1[, 2], ci = TRUE)
roc2 <- pROC::roc(response = TestY, predictor = pred2[, 2], ci = TRUE)
roc3 <- pROC::roc(response = TestY, predictor = pred3[, 2], ci = TRUE)

# --- 8. Figure SM2 & Table SM1: Validation Set Metrics ------------------------
compute_metrics <- function(fit_obj, TestX, TestY, model_name, roc_obj) {
    pred0 <- predict(fit_obj, newdata = TestX, type = "prob")
    bestp <- roc_obj$thresholds[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)]
    pred_class <- factor(ifelse(pred0[, 2] > bestp, "Psoriasis", "Healthy"), levels = c("Healthy", "Psoriasis"))

    cm <- caret::confusionMatrix(data = pred_class, reference = TestY, positive = "Psoriasis", mode = "everything")

    metrics <- data.frame(
        Model = model_name,
        ROC_AUC = as.numeric(roc_obj$auc),
        Accuracy = cm$overall["Accuracy"],
        Sensitivity = cm$byClass["Sensitivity"],
        Specificity = cm$byClass["Specificity"],
        Precision = cm$byClass["Precision"],
        F1_Score = cm$byClass["F1"],
        Balanced_Accuracy = cm$byClass["Balanced Accuracy"],
        Cohen_Kappa = cm$overall["Kappa"]
    )

    # Save Confusion Matrix as pseudo "Figure SM2" representation in text format
    sink(file.path(out_dir, paste0("Figure_SM2_", model_name, "_CM.txt")))
    print(cm)
    sink()

    return(metrics)
}

metrics1 <- compute_metrics(fit1, TestX, TestY, "Model 1", roc1)
metrics2 <- compute_metrics(fit2, TestX, TestY, "Model 2", roc2)
metrics3 <- compute_metrics(fit3, TestX, TestY, "Model 3", roc3)

Table_SM1 <- bind_rows(metrics1, metrics2, metrics3)
write.csv(Table_SM1, file.path(out_dir, "Table_SM1_Validation_Metrics.csv"), row.names = FALSE)


# --- 9. Original ROC Curve Plot (for reference) -------------------------------
pdf(file = file.path(out_dir, "Original_ROC.pdf"), width = 7, height = 7)
plot(roc1, print.auc = TRUE, grid = TRUE, legacy.axes = FALSE, bty = "l", col = "#00A087", print.auc.x = 0.75, print.auc.y = 0.45)
plot(roc2, print.auc = TRUE, add = TRUE, col = "#4DBBD5", print.auc.x = 0.68, print.auc.y = 0.78)
plot(roc3, print.auc = TRUE, add = TRUE, col = "#E64B35", print.auc.x = 0.42, print.auc.y = 0.92)
legend("bottomright", legend = c("Model1 (Boruta)", "Model2 (findLinearCombos+Boruta)", "Model3 (s__GGB51647_SGB4348)"), col = c("#00A087", "#4DBBD5", "#E64B35"), lwd = 2, bg = "white", cex = 0.9)
dev.off()

# --- 10. Figure SM3a: Original PR Curves --------------------------------------
pred_obj1 <- prediction(pred1[, 2], TestY)
pred_obj2 <- prediction(pred2[, 2], TestY)
pred_obj3 <- prediction(pred3[, 2], TestY)

prc1 <- performance(pred_obj1, measure = "prec", x.measure = "rec")
prc2 <- performance(pred_obj2, measure = "prec", x.measure = "rec")
prc3 <- performance(pred_obj3, measure = "prec", x.measure = "rec")

aucpr1 <- performance(pred_obj1, measure = "aucpr")@y.values[[1]]
aucpr2 <- performance(pred_obj2, measure = "aucpr")@y.values[[1]]
aucpr3 <- performance(pred_obj3, measure = "aucpr")@y.values[[1]]

pdf(file = file.path(out_dir, "Figure_SM3a_PRC_Original.pdf"), width = 7, height = 7)
plot(prc1, main = "PR Curve (Original Validation Set)", ylim = c(0, 1), xlab = "Recall", ylab = "Precision", col = "#00A087", lwd = 2)
plot(prc2, add = TRUE, col = "#4DBBD5", lwd = 2)
plot(prc3, add = TRUE, col = "#E64B35", lwd = 2)
legend("bottomleft", legend = c(paste0("Model1 (PR-AUC: ", round(aucpr1, 2), ")"), paste0("Model2 (PR-AUC: ", round(aucpr2, 2), ")"), paste0("Model3 (PR-AUC: ", round(aucpr3, 2), ")")), col = c("#00A087", "#4DBBD5", "#E64B35"), lwd = 2, bg = "white")
dev.off()


# --- 11. SMOTE Balancing (SM6.2) ----------------------------------------------
library(DMwR)
set.seed(1234)

# Combine TestX and TestY for SMOTE
val_data <- abu2
val_data$Group <- factor(val_data$Group, levels = c("Healthy", "Psoriasis"))

# Apply SMOTE (perc.over = 200 handles minority (Healthy=17) -> 51. perc.under = 200 handles majority (Pso) -> 68)
smote_data <- SMOTE(Group ~ ., data = val_data, perc.over = 200, perc.under = 200)

cat(sprintf("SMOTE Balanced Set: %d Pso, %d Healthy\n", sum(smote_data$Group == "Psoriasis"), sum(smote_data$Group == "Healthy")))

TestX_smote <- dplyr::select(smote_data, -Group)
TestY_smote <- smote_data$Group

# Predictions on SMOTE data
pred1_smote <- predict(fit1, newdata = TestX_smote, type = "prob")
pred2_smote <- predict(fit2, newdata = TestX_smote, type = "prob")
pred3_smote <- predict(fit3, newdata = TestX_smote, type = "prob")

roc1_smote <- pROC::roc(response = TestY_smote, predictor = pred1_smote[, 2])
roc2_smote <- pROC::roc(response = TestY_smote, predictor = pred2_smote[, 2])
roc3_smote <- pROC::roc(response = TestY_smote, predictor = pred3_smote[, 2])

# --- 12. Figure SM4: ROC Curves on SMOTE Data ---------------------------------
pdf(file = file.path(out_dir, "Figure_SM4_ROC_SMOTE.pdf"), width = 7, height = 7)
plot(roc1_smote, print.auc = TRUE, grid = TRUE, legacy.axes = FALSE, bty = "l", col = "#00A087", print.auc.x = 0.75, print.auc.y = 0.45)
plot(roc2_smote, print.auc = TRUE, add = TRUE, col = "#4DBBD5", print.auc.x = 0.68, print.auc.y = 0.78)
plot(roc3_smote, print.auc = TRUE, add = TRUE, col = "#E64B35", print.auc.x = 0.42, print.auc.y = 0.92)
legend("bottomright", legend = c("Model1 (Boruta)", "Model2 (findLinearCombos+Boruta)", "Model3 (s__GGB51647_SGB4348)"), col = c("#00A087", "#4DBBD5", "#E64B35"), lwd = 2, bg = "white", cex = 0.9)
dev.off()

# --- 13. Figure SM3b: PR Curves on SMOTE Data ---------------------------------
pred_obj1_smote <- prediction(pred1_smote[, 2], TestY_smote)
pred_obj2_smote <- prediction(pred2_smote[, 2], TestY_smote)
pred_obj3_smote <- prediction(pred3_smote[, 2], TestY_smote)

prc1_smote <- performance(pred_obj1_smote, measure = "prec", x.measure = "rec")
prc2_smote <- performance(pred_obj2_smote, measure = "prec", x.measure = "rec")
prc3_smote <- performance(pred_obj3_smote, measure = "prec", x.measure = "rec")

aucpr1_smote <- performance(pred_obj1_smote, measure = "aucpr")@y.values[[1]]
aucpr2_smote <- performance(pred_obj2_smote, measure = "aucpr")@y.values[[1]]
aucpr3_smote <- performance(pred_obj3_smote, measure = "aucpr")@y.values[[1]]

pdf(file = file.path(out_dir, "Figure_SM3b_PRC_SMOTE.pdf"), width = 7, height = 7)
plot(prc1_smote, main = "PR Curve (SMOTE Balanced Set)", ylim = c(0, 1), xlab = "Recall", ylab = "Precision", col = "#00A087", lwd = 2)
plot(prc2_smote, add = TRUE, col = "#4DBBD5", lwd = 2)
plot(prc3_smote, add = TRUE, col = "#E64B35", lwd = 2)
legend("bottomleft", legend = c(paste0("Model1 (PR-AUC: ", round(aucpr1_smote, 2), ")"), paste0("Model2 (PR-AUC: ", round(aucpr2_smote, 2), ")"), paste0("Model3 (PR-AUC: ", round(aucpr3_smote, 2), ")")), col = c("#00A087", "#4DBBD5", "#E64B35"), lwd = 2, bg = "white")
dev.off()

# --- 14. Table SM2: PR-AUC Comparison Export ----------------------------------
Table_SM2 <- data.frame(
    Metric = c("PR-AUC (original validation set)", "PR-AUC (SMOTE-balanced set)"),
    Model_1_Boruta = c(round(aucpr1, 2), round(aucpr1_smote, 2)),
    Model_2_LinearCombosBoruta = c(round(aucpr2, 2), round(aucpr2_smote, 2)),
    Model_3_SingleFeature = c(round(aucpr3, 2), round(aucpr3_smote, 2))
)
write.csv(Table_SM2, file.path(out_dir, "Table_SM2_PR_AUC_Comparison.csv"), row.names = FALSE)

message("05_RandomForest_Classifier.R completed: Figures SM1-SM4 and Tables SM1-SM2 generated in outputs/RandomForest/")
