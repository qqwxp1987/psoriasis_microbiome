###############################################################################
# 05_RandomForest_Classifier.R
# Random Forest classification model for psoriasis prediction
# - Feature selection: Boruta + findLinearCombos
# - Model training: caret RF with 5-fold x5-repeat CV
# - External validation: ROC curves + DeLong tests
# - Precision-Recall curves (ROCR)
# - Confusion matrix with optimal threshold
# - SMOTE for class imbalance (optional)
#
# Required input files (in scripts/caret_data_code/):
#   train_metaphlan_taxonomic_profiles_20240418.csv  (Discovery set CLR abundance)
#   train_metadata_20240401.txt                      (Discovery set metadata)
#   extra_data_20240418.csv                          (Validation set CLR abundance)
#
# Required R packages:
#   Boruta, caret, dplyr, pROC, patchwork, ggpubr, scales, tidyr, ROCR, DMwR
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(Boruta)
library(caret)
library(dplyr)
library(pROC)
library(patchwork)
library(ggplot2)
library(here)

# --- 1. Helper function: plot single ROC curve --------------------------------
plot1ROC <- function(roc) {
    P <- plot(roc,
        print.auc = TRUE,
        auc.polygon = TRUE,
        grid = TRUE,
        max.auc.polygon = TRUE,
        max.auc.polygon.col = "white",
        auc.polygon.col = "#C9DBEF",
        print.thres = TRUE,
        legacy.axes = FALSE,
        bty = "l"
    )
    return(list(plot = P, ci.auc = ci.auc(roc)))
}

# --- 2. Load training data (Discovery set: 28 pairs) -------------------------
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

# --- 3. Preprocessing function ------------------------------------------------
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

# --- 4. Boruta feature selection-----------------------------------------------
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

# --- 5. Define model formulas -------------------------------------------------
# Model 1: 9 features from Boruta (no collinearity filter)
formula1 <- as.formula("Group ~ s__GGB51647_SGB4348 + s__GGB3746_SGB5089 + s__GGB33085_SGB1673 +
    s__Prevotella_stercorea + s__Clostridium_sp_AF20_17LB + s__Bacteroides_finegoldii +
    s__Bacteroides_congonensis + s__Clostridiaceae_bacterium + s__Bifidobacterium_bifidum")

# Model 2: 3 features from findLinearCombos + Boruta
formula2 <- as.formula("Group ~ s__GGB51647_SGB4348 + s__GGB1109_SGB1423 + s__GGB3746_SGB5089")

# Model 3: Single feature (s__GGB51647_SGB4348 only)
formula3 <- as.formula("Group ~ s__GGB51647_SGB4348")

# --- 6. Train RF models with 5-fold x 5-repeat CV -----------------------------
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

# --- 7. Training metrics visualization ----------------------------------------
library(ggpubr)
library(scales)
library(tidyr)

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

plot_train_metrics(fit1, "vars9_trainset_metric.pdf")
plot_train_metrics(fit2, "vars3_trainset_metric.pdf")
plot_train_metrics(fit3, "vars1_trainset_metric.pdf")

# --- 8. External validation (Validation set) ----------------------------------
abu2 <- read.csv(file.path(data_dir, "extra_data_20240418.csv"),
    sep = ",", header = TRUE, row.names = 1, check.names = FALSE
)
bateria_names2 <- strsplit(colnames(abu2), split = ";|\\|")
yy <- c()
for (t in bateria_names2) {
    yy <- c(yy, t[length(t)])
}
colnames(abu2) <- yy

# Use original (unbalanced) test data
TestX <- dplyr::select(abu2, -Group)
TestY <- factor(abu2[, "Group"], levels = c("Healthy", "Psoriasis"))

# Predictions
pred1 <- predict(fit1, newdata = TestX, type = "prob")
pred2 <- predict(fit2, newdata = TestX, type = "prob")
pred3 <- predict(fit3, newdata = TestX, type = "prob")

# ROC curves
roc1 <- pROC::roc(response = TestY, predictor = pred1[, 2], ci = TRUE)
roc2 <- pROC::roc(response = TestY, predictor = pred2[, 2], ci = TRUE)
roc3 <- pROC::roc(response = TestY, predictor = pred3[, 2], ci = TRUE)

# DeLong tests
message("DeLong test: Model1 vs Model2")
roc.test(roc1, roc2)
message("DeLong test: Model1 vs Model3")
roc.test(roc1, roc3)
message("DeLong test: Model2 vs Model3")
roc.test(roc2, roc3)

# --- 9. ROC curve plot --------------------------------------------------------
pdf(file = file.path(out_dir, "AUC-ROC.pdf"), width = 7, height = 7)
plot(roc1,
    print.auc = TRUE, grid = TRUE,
    max.auc.polygon = TRUE, max.auc.polygon.col = "white",
    print.auc.x = 0.75, print.auc.y = 0.45,
    legacy.axes = FALSE, bty = "l", col = "#00A087"
)
plot(roc2,
    print.auc = TRUE, print.auc.x = 0.68, print.auc.y = 0.78,
    add = TRUE, col = "#4DBBD5"
)
plot(roc3,
    print.auc = TRUE, print.auc.x = 0.42, print.auc.y = 0.92,
    add = TRUE, col = "#E64B35"
)
legend("bottomright",
    legend = c(
        "Model1 (Boruta)", "Model2 (findLinearCombos+Boruta)",
        "Model3 (s__GGB51647_SGB4348)"
    ),
    col = c("#00A087", "#4DBBD5", "#E64B35"),
    lwd = 2, bg = "white", cex = 0.9
)
dev.off()

# --- 10. Precision-Recall curve -----------------------------------------------
library(ROCR)
pred_obj1 <- prediction(pred1[, 2], TestY)
pred_obj2 <- prediction(pred2[, 2], TestY)
pred_obj3 <- prediction(pred3[, 2], TestY)

prc1 <- performance(pred_obj1, measure = "prec", x.measure = "rec")
prc2 <- performance(pred_obj2, measure = "prec", x.measure = "rec")
prc3 <- performance(pred_obj3, measure = "prec", x.measure = "rec")

v1 <- performance(pred_obj1, measure = "aucpr")@y.values[[1]]
v2 <- performance(pred_obj2, measure = "aucpr")@y.values[[1]]
v3 <- performance(pred_obj3, measure = "aucpr")@y.values[[1]]

pdf(file = file.path(out_dir, "AUC-PRC.pdf"), width = 7, height = 7)
plot(prc1,
    main = "Precision-Recall Curve (AUC-PR)", ylim = c(0, 1),
    xlab = "Recall", ylab = "Precision", col = "#00A087", lwd = 2
)
plot(prc2, add = TRUE, col = "#4DBBD5", lwd = 2)
plot(prc3, add = TRUE, col = "#E64B35", lwd = 2)
legend("bottomleft",
    legend = c(
        paste0("Model1 (Boruta, AUC-PRC: ", round(v1, 4), ")"),
        paste0("Model2 (findLinearCombos+Boruta, AUC-PRC: ", round(v2, 4), ")"),
        paste0("Model3 (s__GGB51647_SGB4348, AUC-PRC: ", round(v3, 4), ")")
    ),
    col = c("#00A087", "#4DBBD5", "#E64B35"),
    lwd = 2, bg = "white", cex = 0.9
)
dev.off()

# --- 11. Confusion matrix (using optimal Youden threshold) --------------------
compute_cm <- function(fit_obj, TestX, TestY, model_name) {
    pred0 <- predict(fit_obj, newdata = TestX, type = "prob")
    roc_tmp <- pROC::roc(response = TestY, predictor = pred0[, 2])
    bestp <- roc_tmp$thresholds[which.max(roc_tmp$sensitivities + roc_tmp$specificities - 1)]
    pred <- as.factor(ifelse(pred0[, 2] > bestp,
        levels(as.factor(TestY))[2],
        levels(as.factor(TestY))[1]
    ))
    cm <- caret::confusionMatrix(
        data = pred,
        reference = as.factor(TestY),
        positive = levels(as.factor(TestY))[2],
        mode = "everything"
    )
    message(paste("--- Confusion Matrix for", model_name, "---"))
    print(cm)
    return(cm)
}

cm1 <- compute_cm(fit1, TestX, TestY, "Model1")
cm2 <- compute_cm(fit2, TestX, TestY, "Model2")
cm3 <- compute_cm(fit3, TestX, TestY, "Model3")
