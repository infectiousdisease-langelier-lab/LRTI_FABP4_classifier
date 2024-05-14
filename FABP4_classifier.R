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
counts <- read.csv(
  "adult_count/merged_counts.csv",
  row.names=1, check.names=FALSE)

# Metadata
metadata <- read.csv(
  "adult_count/merged_metadata.csv")
metadata$patient_id <- as.character(metadata$patient_id)
metadata$LRTI <- factor(metadata$LRTI, levels=c("0","1"))
metadata$batch <- factor(metadata$batch)

# Number of samples per LRTI status
print(table(metadata$LRTI))

gene.symbol <- counts[,"gene_symbol",drop=FALSE]

## Pediatrics -----
counts.ped <- read.csv(
  "GSE212532_gene_counts_GEO.csv",
  row.names=1)
stopifnot(rownames(counts.ped)==rownames(counts))

# Metadata
metadata.ped <- read.csv(
  "GSE212532_sample_metadata_GEO.csv"
)

# Only keep definite and no evidence
metadata.ped <- metadata.ped %>%
  subset(LRTI_adjudication %in% c("Definite","No Evidence"))
metadata.ped$LRTI <- factor(ifelse(
  metadata.ped$LRTI_adjudication=="Definite", "1", "0"),
  levels=c("0","1"))
counts.ped <- counts.ped[,metadata.ped$sample_name]

# Number of samples per LRTI status
print(table(metadata.ped$LRTI))

# limma DE -----
## Children -----
keep <- rowSums(counts.ped>=10) >= (0.2*ncol(counts.ped))

design <- model.matrix(
  ~ LRTI,
  data = metadata.ped)
print(colnames(design))

# limma-voom
vwts <- voom(counts.ped[keep,], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
res.1 <- topTable(
  vfit,
  coef = "LRTI1", sort.by = "none", 
  number = Inf)
print("LRTI vs no LRTI in children")
print(sprintf("Number of DE genes with FDR < 0.05: %d", sum(res.1$adj.P.Val<0.05)))

# Add gene symbols
res.1$gene_symbol <- gene.symbol[rownames(res.1),"gene_symbol"]

# Check FABP4 gene
print(res.1 %>%
        subset(gene_symbol=="FABP4"))

# Volcano plot
ggplot(data=res.1,
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
  res.1,
  "DE_children.csv"
)

## Adults -----
counts.all <- counts[,metadata$patient_id]
keep <- rowSums(counts.all>=10) >= (0.2*ncol(counts.all))

design <- model.matrix(
  ~ LRTI + batch,
  data = metadata)
print(colnames(design))

# limma-voom
vwts <- voom(counts.all[keep,], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
res.2 <- topTable(
  vfit,
  coef = "LRTI1", sort.by = "none", 
  number = Inf)
print("LRTI vs no LRTI in adults")
print(sprintf("Number of DE genes with FDR < 0.05: %d", sum(res.2$adj.P.Val<0.05)))

# Add gene symbols
res.2$gene_symbol <- gene.symbol[rownames(res.2),"gene_symbol"]

# Check FABP4 gene
print(res.2 %>%
        subset(gene_symbol=="FABP4"))

write.csv(
  res.2,
  "DE_adults.csv"
)

# Adults -----
## Generate 5 folds -----
# Generate folds for 5-fold cross-validation
# Try to approximate the ratio of LRTI to no LRTI sample in each fold
min.LRTI <- floor(sum(metadata$LRTI=="1")/5) # minimum number of LRTI sample per fold
set.seed(0)
while (TRUE) {
  # Generate 5 fold
  cv.folds <- metadata %>%
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

## 5-fold CV direct ROC with VST -----
# Follows https://stackoverflow.com/a/56525428 and https://daviddalpiaz.github.io/r4sl/logistic-regression.html#roc-curves

# Get ensembl ID of FABP4
gene.id <- rownames(gene.symbol %>% subset(gene_symbol=="FABP4"))

# For storing ROC results
fold.roc <- list()
# For storing cutoff value, sensitivity and specificity according to Youden
youden <- list()
# For storing cutoff value, sensitivity and specificity at 90% sensitivity
sens90 <- list()
for (k in c(1:5)) {
  test.fold <- cv.folds[cv.folds$fold==k,]
  train.folds <- cv.folds[cv.folds$fold!=k,]
  
  # Train data
  counts.train <- counts.all[,train.folds$patient_id]
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
  counts.test <- counts.all[,test.fold$patient_id]
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
  
  fold.roc[[k]] <- pROC::roc(
    test.fold$LRTI ~ test.fold$FABP4,
    levels=c("0","1"),
    direction=">",
    plot=TRUE, print.auc=TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", k, fold.roc[[k]]$auc))
  
  # Youden's index
  youden[[k]] <- pROC::coords(
    fold.roc[[k]], x="best", best.method="youden",
    ret=c("threshold","specificity","sensitivity","accuracy","precision"))
  
  # 90% sensitivity
  sens90[[k]] <- pROC::coords(
    fold.roc[[k]], x=0.9, input="sensitivity",
    ret=c("threshold","specificity","sensitivity","accuracy","precision"))
}

# Get the list of AUCs
fold.auc <- unlist(lapply(
  fold.roc,
  FUN=function(x) x$auc))

# Mean
print(mean(fold.auc))
# SD
print(sd(fold.auc))

# Convert the youden list to dataframe
youden.df <- Reduce(rbind, youden)
rownames(youden.df) <- 1:5

# Convert the sens90 list to dataframe
sens90.df <- Reduce(rbind, sens90)
rownames(sens90.df) <- 1:5

## Mean ROC curve -----
# Interpolate the ROC to calculate the mean AUC, taking inspiration from https://stats.stackexchange.com/a/187003
# Reverse the FPR and TPR order of pROC::roc's output, so that it works with approx(ties="ordered")
roc.approx <- data.frame(
  fpr.out=seq(0, 1, length.out=100),
  tpr1=0,
  tpr2=0,
  tpr3=0,
  tpr4=0,
  tpr5=0)
for (k in 1:5) {
  fpr <- rev(1-fold.roc[[k]]$specificities)
  tpr <- rev(fold.roc[[k]]$sensitivities)
  tpr.out <- approx(fpr, tpr, xout=roc.approx$fpr.out,
                    method="linear", ties="ordered")
  roc.approx[,k+1] <- tpr.out$y
  plot(fpr, tpr, main=k)
  points(tpr.out$x, tpr.out$y, col="red")
}

roc.approx$tpr.mean <- rowMeans(roc.approx[,2:6])
roc.approx[1,"tpr.mean"] <- 0 # Force the mean ROC to start at (0,0)
plot(roc.approx$fpr.out, roc.approx$tpr.mean)

# Children -----
## Generate 5 folds -----
# Generate folds for 5-fold cross-validation
# Try to approximate the ratio of LRTI to no LRTI sample in each fold
min.LRTI <- floor(sum(metadata.ped$LRTI=="0")/5) # minimum number of no LRTI sample per fold
set.seed(0)
while (TRUE) {
  # Generate 5 fold
  cv.folds.ped <- metadata.ped %>%
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

## 5-fold CV direct ROC with VST -----
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
  counts.train <- counts.ped[,train.folds$sample_name]
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
  counts.test <- counts.ped[,test.fold$sample_name]
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

# Mean
print(mean(fold.auc.ped))
# SD
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
  tpr.out <- approx(fpr, tpr, xout=roc.ped.approx$fpr.out,
                    method="linear", ties="ordered")
  roc.ped.approx[,k+1] <- tpr.out$y
  plot(fpr, tpr, main=k)
  points(tpr.out$x, tpr.out$y, col="red")
}

roc.ped.approx$tpr.mean <- rowMeans(roc.ped.approx[,2:6])
roc.ped.approx[1,"tpr.mean"] <- 0 # Force the mean ROC to start at (0,0)
plot(roc.ped.approx$fpr.out, roc.ped.approx$tpr.mean)

# Pretty plot -----
p <- ggplot()

# Adults
for (k in 1:5) {
  p <- p + geom_line(
    data=data.frame(
      x=rev(1-fold.roc[[k]]$specificities),
      y=rev(fold.roc[[k]]$sensitivities)),
    aes(x=x,y=y), color="#1b9e77", alpha=0.4, linewidth=0.4)
}
p <- p + geom_line(
  data=roc.approx,
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
