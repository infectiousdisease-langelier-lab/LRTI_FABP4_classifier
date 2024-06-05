# Here I will analyze mBAL data to test how good FABP4 is as a single-gene classifier of LRTI

library(tidyverse)
library(ggplot2)
library(patchwork)

library(limma)
library(DESeq2)

library(pROC)

# Import data -----
## Adults -----
# Counts
counts.adult <- read.csv(
  "./adult_count/merged_counts.csv",
  row.names=1, check.names=FALSE)

# Metadata
metadata.adult <- read.csv(
  "./adult_count/merged_metadata.csv")
metadata.adult$patient_id <- as.character(metadata.adult$patient_id)
metadata.adult$batch <- as.factor(metadata.adult$batch)

# Number of samples per LRTI status
print(table(metadata.adult$LRTI_adj))

gene.symbol <- counts.adult[,"gene_symbol",drop=FALSE]

## Pediatrics -----
counts.ped <- read.csv(
  "GSE212532_gene_counts_GEO.csv",
  row.names=1)
stopifnot(rownames(counts.ped)==rownames(counts))

# Metadata
metadata.ped <- read.csv(
  "GSE212532_sample_metadata_GEO.csv"
)

# Number of samples per LRTI status
print(table(metadata.ped$LRTI_adj))

# limma DE -----
## Children -----
# Only keep definite and no evidence for DE
metadata.ped.DE <- metadata.ped %>%
  subset(LRTI_adjudication %in% c("Definite","No Evidence"))
metadata.ped.DE$LRTI <- factor(ifelse(
  metadata.ped.DE$LRTI_adjudication=="Definite", "1", "0"),
  levels=c("0","1"))
counts.ped.DE <- counts.ped[,metadata.ped.DE$sample_name]

# Number of samples per LRTI status
print(table(metadata.ped.DE$LRTI))

# Gene filtering
keep <- rowSums(counts.ped.DE>=10) >= (0.2*ncol(counts.ped.DE))

# Convert to CPM for plotting later
cpm.ped <- edgeR::cpm(
  counts.ped.DE[keep,], log=TRUE
)

design <- model.matrix(
  ~ LRTI,
  data = metadata.ped.DE)
print(colnames(design))

# limma-voom
vwts <- voom(counts.ped.DE[keep,], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
res.ped <- topTable(
  vfit,
  coef = "LRTI1", sort.by = "none", 
  number = Inf)
print("LRTI vs no LRTI in children")
print(sprintf("Number of DE genes with FDR < 0.05: %d", sum(res.ped$adj.P.Val<0.05)))

# Add gene symbols
res.ped$gene_symbol <- gene.symbol[rownames(res.ped),"gene_symbol"]

# Check FABP4 gene
print(res.ped %>%
        subset(gene_symbol=="FABP4"))

# Volcano plot
ggplot(data=res.ped,
       aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=factor(adj.P.Val<0.05)), size=0.3) +
  geom_point(data=. %>% subset(gene_symbol=="FABP4"),
             color="#d95f02", size=1.5) +
  scale_color_manual(values=c("TRUE"="black","FALSE"="gray75"),
                     guide="none") +
  scale_x_continuous(limits=c(-5.5,5.5), expand=c(0,0)) +
  labs(x="log2(fold change)",
       y="-log10(adjusted P-value)") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    # panel.border = element_rect(linewidth=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )
ggsave(
  "volcano_children.png",
  width=3.3, height=3.4, units="in", dpi=1200)

write.csv(
  res.ped,
  "DE_children.csv"
)

## Adults -----
# Only keep definite and no evidence for DE
metadata.adult.DE <- metadata.adult %>%
  subset(LRTI_adj %in% c("Definite","No Evidence"))
metadata.adult.DE$LRTI <- factor(ifelse(
  metadata.adult.DE$LRTI_adj=="Definite", "1", "0"),
  levels=c("0","1"))
counts.adult.DE <- counts.adult[,metadata.adult.DE$patient_id]

# Number of samples per LRTI status
print(table(metadata.adult.DE$LRTI))

keep <- rowSums(counts.adult.DE>=10) >= (0.2*ncol(counts.adult.DE))

# Convert to CPM for plotting later
cpm.adult <- edgeR::cpm(
  counts.adult.DE[keep,], log=TRUE
)

design <- model.matrix(
  ~ LRTI + batch,
  data = metadata.adult.DE)
print(colnames(design))

# limma-voom
vwts <- voom(counts.adult.DE[keep,], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
res.adult <- topTable(
  vfit,
  coef = "LRTI1", sort.by = "none", 
  number = Inf)
print("LRTI vs no LRTI in adults")
print(sprintf("Number of DE genes with FDR < 0.05: %d", sum(res.adult$adj.P.Val<0.05)))

# Add gene symbols
res.adult$gene_symbol <- gene.symbol[rownames(res.adult),"gene_symbol"]

# Check FABP4 gene
print(res.adult %>%
        subset(gene_symbol=="FABP4"))

write.csv(
  res.adult,
  "DE_adults.csv"
)

# Children 5-fold CV -----
## Generate 5 folds -----
# Generate folds for 5-fold cross-validation
# Try to approximate the ratio of LRTI to no LRTI sample in each fold
min.LRTI <- floor(sum(metadata.ped.DE$LRTI=="0")/5) # minimum number of no LRTI sample per fold
set.seed(0)
while (TRUE) {
  # Generate 5 fold
  cv.folds.ped <- metadata.ped.DE %>%
    select(sample_name, LRTI) %>%
    mutate(fold=sample(rep(1:5, length.out=nrow(.))))
  
  # Count number of LRTI samples per fold
  cv.folds.ped.table <- cv.folds.ped %>%
    group_by(fold) %>%
    dplyr::count(LRTI)
  
  print("Generated CV folds with following counts:")
  print(cv.folds.ped.table)
  
  if (min(cv.folds.ped.table[cv.folds.ped.table$LRTI=="0","n"]) < min.LRTI) {
    print("At least one fold has too few LRTI samples. Regenerating CV folds...")
  } else {
    break
  }
}

## 5-fold CV normalized FABP4 expression -----
# Get ensembl ID of FABP4
gene.id <- rownames(gene.symbol %>% subset(gene_symbol=="FABP4"))

# For storing ROC results
fold.roc.ped <- list()
# For storing cutoff value, sensitivity and specificity according to Youden
youden.ped <- list()
# For storing cutoff value, sensitivity and specificity at 90% sensitivity
sens90.ped <- list()
for (k in c(1:5)) {
  test.fold <- cv.folds.ped[cv.folds.ped$fold==k,]
  train.folds <- cv.folds.ped[cv.folds.ped$fold!=k,]
  
  # Train data
  counts.train <- counts.ped.DE[,train.folds$sample_name]
  keep.train <- rowSums(counts.train>=10) >= (0.2*ncol(counts.train))
  dds.train <- DESeqDataSetFromMatrix(
    countData = counts.train[keep.train,],
    colData = train.folds,
    design = ~1)
  dds.train <- estimateSizeFactors(dds.train)
  dds.train <- estimateDispersions(dds.train)
  vsd.train <- varianceStabilizingTransformation(dds.train) %>% 
    assay %>% 
    round(., digits=2)
  # Add FABP4 expression
  train.folds$FABP4 <- vsd.train[gene.id,]
  
  # Test data
  counts.test <- counts.ped.DE[,test.fold$sample_name]
  dds.test <- DESeqDataSetFromMatrix(
    countData = counts.test[keep.train,],
    colData = test.fold,
    design = ~1)
  dds.test <- estimateSizeFactors(dds.test)
  dispersionFunction(dds.test) <- dispersionFunction(dds.train) # assign the dispersion function from the training data alone
  vsd.test <- varianceStabilizingTransformation(dds.test, blind=FALSE) %>% 
    assay %>% 
    round(., digits=2)
  # Add FABP4 expression
  test.fold$FABP4 <- vsd.test[gene.id,]
  
  fold.roc.ped[[k]] <- pROC::roc(
    test.fold$LRTI ~ test.fold$FABP4,
    levels=c("0","1"),
    direction=">",
    plot=TRUE, print.auc=TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", k, fold.roc.ped[[k]]$auc))
  
  # Youden's index
  youden.ped[[k]] <- pROC::coords(
    fold.roc.ped[[k]], x="best", best.method="youden",
    ret=c("threshold","specificity","sensitivity","accuracy","precision"))
  
  # 90% sensitivity
  sens90.ped[[k]] <- pROC::coords(
    fold.roc.ped[[k]], x=0.9, input="sensitivity",
    ret=c("threshold","specificity","sensitivity","accuracy","precision"))
}

# Get the list of AUCs
fold.auc.ped <- unlist(lapply(
  fold.roc.ped,
  FUN=function(x) x$auc))

# Mean: 0.9006159
print(mean(fold.auc.ped))
# SD: 0.06795188
print(sd(fold.auc.ped))

# Convert the youden list to dataframe
youden.ped.df <- Reduce(rbind, youden.ped)
rownames(youden.ped.df) <- 1:5

# Convert the sens90 list to dataframe
sens90.ped.df <- Reduce(rbind, sens90.ped)
rownames(sens90.ped.df) <- 1:5

## Mean ROC curve -----
# Interpolate the ROC to calculate the mean AUC, taking inspiration from https://stats.stackexchange.com/a/187003
# Reverse the FPR and TPR order of pROC::roc's output, so that it works with approx(ties="ordered")
roc.ped.approx <- data.frame(
  fpr.out=seq(0, 1, length.out=100),
  tpr1=0,
  tpr2=0,
  tpr3=0,
  tpr4=0,
  tpr5=0)
for (k in 1:5) {
  fpr <- rev(1-fold.roc.ped[[k]]$specificities)
  tpr <- rev(fold.roc.ped[[k]]$sensitivities)
  tpr.out <- approx(fpr, tpr,
                    xout=roc.ped.approx$fpr.out,
                    method="linear", ties="ordered")
  roc.ped.approx[,k+1] <- tpr.out$y
  plot(fpr, tpr, main=k)
  points(tpr.out$x, tpr.out$y, col="red")
}

roc.ped.approx$tpr.mean <- rowMeans(roc.ped.approx[,2:6])
roc.ped.approx[1,"tpr.mean"] <- 0 # Force the mean ROC to start at (0,0)
plot(roc.ped.approx$fpr.out, roc.ped.approx$tpr.mean)

# Adults 5-fold CV -----
## Generate 5 folds -----
# Generate folds for 5-fold cross-validation
# Try to approximate the ratio of LRTI to no LRTI sample in each fold
min.LRTI <- floor(sum(metadata.adult.DE$LRTI=="1")/5) # minimum number of LRTI sample per fold
set.seed(0)
while (TRUE) {
  # Generate 5 fold
  cv.folds <- metadata.adult.DE %>%
    select(patient_id, batch, LRTI) %>%
    mutate(fold=sample(rep(1:5, length.out=nrow(.))))
  
  # Count number of LRTI samples per fold
  cv.folds.table <- cv.folds %>%
    group_by(fold) %>%
    dplyr::count(LRTI)
  
  print("Generated CV folds with following counts:")
  print(cv.folds.table)
  
  if (min(cv.folds.table[cv.folds.table$LRTI=="1","n"]) < min.LRTI) {
    print("At least one fold has too few LRTI samples. Regenerating CV folds...")
  } else {
    break
  }
}
# LRTI status vs batch: all three batches are well presented in each of the 5 folds
print(cv.folds %>%
        group_by(fold) %>%
        dplyr::count(batch))

## 5-fold CV with normalized FABP4 expression -----
# Follows https://stackoverflow.com/a/56525428 and https://daviddalpiaz.github.io/r4sl/logistic-regression.html#roc-curves

# Get ensembl ID of FABP4
gene.id <- rownames(gene.symbol %>% subset(gene_symbol=="FABP4"))

# For storing ROC results
fold.roc.adult <- list()
# For storing cutoff value, sensitivity and specificity according to Youden
youden <- list()
# For storing cutoff value, sensitivity and specificity at 90% sensitivity
sens90 <- list()
for (k in c(1:5)) {
  test.fold <- cv.folds[cv.folds$fold==k,]
  train.folds <- cv.folds[cv.folds$fold!=k,]
  
  # Train data
  counts.train <- counts.adult.DE[,train.folds$patient_id]
  keep.train <- rowSums(counts.train>=10) >= (0.2*ncol(counts.train))
  dds.train <- DESeqDataSetFromMatrix(
    countData = counts.train[keep.train,],
    colData = train.folds,
    design = ~1)
  dds.train <- estimateSizeFactors(dds.train)
  dds.train <- estimateDispersions(dds.train)
  vsd.train <- varianceStabilizingTransformation(dds.train) %>% 
    assay %>% 
    round(., digits=2)
  # Add FABP4 expression
  train.folds$FABP4 <- vsd.train[gene.id,]
  
  # Test data
  counts.test <- counts.adult.DE[,test.fold$patient_id]
  dds.test <- DESeqDataSetFromMatrix(
    countData = counts.test[keep.train,],
    colData = test.fold,
    design = ~1)
  dds.test <- estimateSizeFactors(dds.test)
  dispersionFunction(dds.test) <- dispersionFunction(dds.train) # assign the dispersion function from the training data alone
  vsd.test <- varianceStabilizingTransformation(dds.test, blind=FALSE) %>% 
    assay %>% 
    round(., digits=2)
  # Add FABP4 expression
  test.fold$FABP4 <- vsd.test[gene.id,]
  
  fold.roc.adult[[k]] <- pROC::roc(
    test.fold$LRTI ~ test.fold$FABP4,
    levels=c("0","1"),
    direction=">",
    plot=TRUE, print.auc=TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", k, fold.roc.adult[[k]]$auc))
  
  # Youden's index
  youden[[k]] <- pROC::coords(
    fold.roc.adult[[k]], x="best", best.method="youden",
    ret=c("threshold","specificity","sensitivity","accuracy","precision"))
  
  # 90% sensitivity
  sens90[[k]] <- pROC::coords(
    fold.roc.adult[[k]], x=0.9, input="sensitivity",
    ret=c("threshold","specificity","sensitivity","accuracy","precision"))
}

# Get the list of AUCs
fold.auc.adult <- unlist(lapply(
  fold.roc.adult,
  FUN=function(x) x$auc))

# Mean: 0.8527778
print(mean(fold.auc.adult))
# SD: 0.1205382
print(sd(fold.auc.adult))

# Convert the youden list to dataframe
youden.adult.df <- Reduce(rbind, youden)
rownames(youden.adult.df) <- 1:5

# Convert the sens90 list to dataframe
sens90.adult.df <- Reduce(rbind, sens90)
rownames(sens90.adult.df) <- 1:5

## Mean ROC curve -----
# Interpolate the ROC to calculate the mean AUC, taking inspiration from https://stats.stackexchange.com/a/187003
# Reverse the FPR and TPR order of pROC::roc's output, so that it works with approx(ties="ordered")
roc.adult.approx <- data.frame(
  fpr.out=seq(0, 1, length.out=100),
  tpr1=0,
  tpr2=0,
  tpr3=0,
  tpr4=0,
  tpr5=0)
for (k in 1:5) {
  fpr <- rev(1-fold.roc.adult[[k]]$specificities)
  tpr <- rev(fold.roc.adult[[k]]$sensitivities)
  tpr.out <- approx(fpr, tpr,
                    xout=roc.adult.approx$fpr.out,
                    method="linear", ties="ordered")
  roc.adult.approx[,k+1] <- tpr.out$y
  plot(fpr, tpr, main=k)
  points(tpr.out$x, tpr.out$y, col="red")
}

roc.adult.approx$tpr.mean <- rowMeans(roc.adult.approx[,2:6])
roc.adult.approx[1,"tpr.mean"] <- 0 # Force the mean ROC to start at (0,0)
plot(roc.adult.approx$fpr.out, roc.adult.approx$tpr.mean)



# Pretty plot 5-fold CV -----
p <- ggplot()

# Adults
for (k in 1:5) {
  p <- p + geom_line(
    data=data.frame(
      x=rev(1-fold.roc.adult[[k]]$specificities),
      y=rev(fold.roc.adult[[k]]$sensitivities)),
    aes(x=x,y=y), color="#1b9e77", alpha=0.4, linewidth=0.4)
}
p <- p + geom_line(
  data=roc.adult.approx,
  aes(x=fpr.out, y=tpr.mean),
  col="#1b9e77", linewidth=1)

# Pediatrics
for (k in 1:5) {
  p <- p + geom_line(
    data=data.frame(
      x=rev(1-fold.roc.ped[[k]]$specificities),
      y=rev(fold.roc.ped[[k]]$sensitivities)),
    aes(x=x,y=y), color="#d95f02", alpha=0.4, linewidth=0.4)
}
p <- p + geom_line(
  data=roc.ped.approx,
  aes(x=fpr.out, y=tpr.mean),
  col="#d95f02", linewidth=1)

# Add the diagonal line
p <- p +
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               linetype="dashed", col="black", linewidth=0.4)

# Formatting
pp <- p +
  labs(x="False positive rate\n(1 - specificity)",
       y="True positive rate\n(sensitivity)") +
  xlim(-0.01,1.01) +
  ylim(-0.01,1.01) +
  expand_limits(x = 0, y = 0) +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )
ggsave(
  "mean_roc.svg",
  plot=pp,
  width=3.65, height=3.6, units="in")

# Classify Suspected & Indeterminate with logistic regression -----
# Get ensembl ID of FABP4
gene.id <- rownames(gene.symbol %>% subset(gene_symbol=="FABP4"))

## Children -----
# Extract Suspected and Indeterminate
metadata.ped.test <- metadata.ped %>%
  subset(LRTI_adjudication %in% c("Suspected","Indeterminate"))
print(table(metadata.ped.test$LRTI_adjudication))

# Train data: all Definite & No Evidence
counts.train <- counts.ped[,metadata.ped.DE$sample_name]
keep.train <- rowSums(counts.train>=10) >= (0.2*ncol(counts.train))
dds.train <- DESeqDataSetFromMatrix(
  countData = counts.train[keep.train,],
  colData = metadata.ped.DE,
  design = ~1)
dds.train <- estimateSizeFactors(dds.train)
dds.train <- estimateDispersions(dds.train)
vsd.train <- varianceStabilizingTransformation(dds.train) %>% 
  assay %>% 
  round(., digits=2)
# Save FABP4's normalized expression
metadata.ped.DE$FABP4 <- vsd.train[gene.id,]

# Fit logistic regression
mod.fit <- glm(
  LRTI ~ FABP4,
  data=metadata.ped.DE,
  family="binomial")

# Test data:
counts.test <- counts.ped[,metadata.ped.test$sample_name]
dds.test <- DESeqDataSetFromMatrix(
  countData = counts.test[keep.train,],
  colData = metadata.ped.test,
  design = ~1)
dds.test <- estimateSizeFactors(dds.test)
dispersionFunction(dds.test) <- dispersionFunction(dds.train) # assign the dispersion function from the training data alone
vsd.test <- varianceStabilizingTransformation(dds.test, blind=FALSE) %>% 
  assay %>% 
  round(., digits=2)
# Save FABP4's normalized expression
metadata.ped.test$FABP4 <- vsd.test[gene.id,]

# Predicted probability of LRTI
metadata.ped.test$prob <- predict(
  mod.fit,
  newdata=metadata.ped.test,
  type="response")
metadata.ped.test$prob.LRTI <- ifelse(
  metadata.ped.test$prob>0.5, "1", "0") # 1 is predicted LRTI

# Number of predicted LRTI for Suspected and Indeterminate
print(table(metadata.ped.test$LRTI_adjudication,
            metadata.ped.test$prob.LRTI))

## Adults -----
# Extract Suspected and Indeterminate
metadata.adult.test <- metadata.adult %>%
  subset(LRTI_adj %in% c("Suspected","Indeterminate"))
print(table(metadata.adult.test$LRTI_adj))

# Train data: all Definite & No Evidence
counts.train <- counts.adult[,metadata.adult.DE$patient_id]
keep.train <- rowSums(counts.train>=10) >= (0.2*ncol(counts.train))
dds.train <- DESeqDataSetFromMatrix(
  countData = counts.train[keep.train,],
  colData = metadata.adult.DE,
  design = ~1)
dds.train <- estimateSizeFactors(dds.train)
dds.train <- estimateDispersions(dds.train)
vsd.train <- varianceStabilizingTransformation(dds.train) %>% 
  assay %>% 
  round(., digits=2)
# Save FABP4's normalized expression
metadata.adult.DE$FABP4 <- vsd.train[gene.id,]

# Fit logistic regression
mod.fit <- glm(
  LRTI ~ FABP4,
  data=metadata.adult.DE,
  family="binomial")

# Test data:
counts.test <- counts.adult[,metadata.adult.test$patient_id]
dds.test <- DESeqDataSetFromMatrix(
  countData = counts.test[keep.train,],
  colData = metadata.adult.test,
  design = ~1)
dds.test <- estimateSizeFactors(dds.test)
dispersionFunction(dds.test) <- dispersionFunction(dds.train) # assign the dispersion function from the training data alone
vsd.test <- varianceStabilizingTransformation(dds.test, blind=FALSE) %>% 
  assay %>% 
  round(., digits=2)
# Save FABP4's normalized expression
metadata.adult.test$FABP4 <- vsd.test[gene.id,]

# Predicted probability of LRTI
metadata.adult.test$prob <- predict(
  mod.fit,
  newdata=metadata.adult.test,
  type="response")
metadata.adult.test$prob.LRTI <- ifelse(
  metadata.adult.test$prob>0.5, "1", "0") # 1 is predicted LRTI

# Number of predicted LRTI for Suspected and Indeterminate
print(table(metadata.adult.test$LRTI_adj,
            metadata.adult.test$prob.LRTI))

# Verify that predicted group "1" corresponds to LRTI
test <- predict(
  mod.fit,
  newdata=metadata.adult.DE,
  type="response")
print(table(metadata.adult.DE$LRTI_adj, test>0.5))
