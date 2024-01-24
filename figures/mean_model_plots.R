set.seed(10)
source("../analysis-code/classifier/coefplot.R")
library(caret)

## Read in dataset

all_with_extrinsic <- read.table("../analysis-code/classifier/mean_all_with_ext_features.tsv", sep = "\t", head = T)
all_without_extrinsic <- read.table("../analysis-code/classifier/mean_all_features.tsv", sep = "\t", head = T)
print(nrow(all_with_extrinsic))
print(nrow(all_without_extrinsic))

## Create holdout set

holdout <- read.table("../analysis-code/classifier/mean_LP3_intrinsic_holdout.tsv", head = T, sep = "\t")
holdout <- holdout[holdout$top_bottom_mean != 3, ]
print(table(holdout$top_bottom_mean))

## Restrict to high low noise locations

print(table(all_with_extrinsic$top_bottom_mean))
print(table(all_without_extrinsic$top_bottom_mean))
all_with_extrinsic <- all_with_extrinsic[all_with_extrinsic$top_bottom_mean != 3, ]
all_without_extrinsic <- all_without_extrinsic[all_without_extrinsic$top_bottom_mean != 3, ]
print(table(all_with_extrinsic$top_bottom_mean))
print(table(all_without_extrinsic$top_bottom_mean))

## Plot histograms of features


## All Intrinsic

print("All intrinsic")
m02 <- glm(top_bottom_mean ~ 1 + H2AZ+ H3K27ac+ H3K27me3+ H3K36me3+ H3K4me1+ H3K4me2+ H3K4me3+ H3K79me2+ H3K9ac+ H3K9me3+ H4K20me1+ ZNF76_HUMAN.H11MO.0.C + RARA_HUMAN.H11MO.0.A + OSR2_HUMAN.H11MO.0.C + NR1D1_HUMAN.H11MO.0.B + TAF1_HUMAN.H11MO.0.A + NFIA_HUMAN.H11MO.0.C + RFX5_HUMAN.H11MO.0.A + ZBTB6_HUMAN.H11MO.0.C + E2F7_HUMAN.H11MO.0.B + MAFF_HUMAN.H11MO.0.B + SMAD4_HUMAN.H11MO.0.B + RORA_HUMAN.H11MO.0.C + BACH1_HUMAN.H11MO.0.A + SMAD3_HUMAN.H11MO.0.B + MTF1_HUMAN.H11MO.0.C + E2F3_HUMAN.H11MO.0.A + ZFX_HUMAN.H11MO.0.A + FEV_HUMAN.H11MO.0.B + ZEB1_HUMAN.H11MO.0.A + MITF_HUMAN.H11MO.0.A + RXRB_HUMAN.H11MO.0.C + SNAI1_HUMAN.H11MO.0.C + IKZF1_HUMAN.H11MO.0.C + SOX5_HUMAN.H11MO.0.C + GC_content+ atac_count+ in_gene+ n_enhancers, family = "binomial", data = all_without_extrinsic)
print(summary(m02))
sum(round(m02$fitted.values) == all_without_extrinsic$top_bottom_mean)/nrow(all_without_extrinsic)
holdout_predictions <- predict(m02, holdout, type = "response")
sum(round(holdout_predictions) == holdout$top_bottom_mean)/nrow(holdout)
plot_coeffs_S(m02, "SuppFigure3F.pdf", title = "mean: intrinsic")

## Cross validation on Intrinsic only model

#Cross validation on Intrinsic only model

#specify the cross-validation method
ctrl <- trainControl(method = "LOOCV")

all_without_extrinsic$mean_level[all_without_extrinsic$top_bottom_mean==1] <- "High Mean"
all_without_extrinsic$mean_level[all_without_extrinsic$top_bottom_mean==0] <- "Low Mean"
all_without_extrinsic$mean_level <- factor(all_without_extrinsic$mean_level, levels = c("Low Mean", "High Mean"))

#fit a regression model and use LOOCV to evaluate performance
all_without_extrinsic_model <- train(mean_level ~ 1 + H2AZ+ H3K27ac+ H3K27me3+ H3K36me3+ H3K4me1+ H3K4me2+ H3K4me3+ H3K79me2+ H3K9ac+ H3K9me3+ H4K20me1+ ZNF76_HUMAN.H11MO.0.C + RARA_HUMAN.H11MO.0.A + OSR2_HUMAN.H11MO.0.C + NR1D1_HUMAN.H11MO.0.B + TAF1_HUMAN.H11MO.0.A + NFIA_HUMAN.H11MO.0.C + RFX5_HUMAN.H11MO.0.A + ZBTB6_HUMAN.H11MO.0.C + E2F7_HUMAN.H11MO.0.B + MAFF_HUMAN.H11MO.0.B + SMAD4_HUMAN.H11MO.0.B + RORA_HUMAN.H11MO.0.C + BACH1_HUMAN.H11MO.0.A + SMAD3_HUMAN.H11MO.0.B + MTF1_HUMAN.H11MO.0.C + E2F3_HUMAN.H11MO.0.A + ZFX_HUMAN.H11MO.0.A + FEV_HUMAN.H11MO.0.B + ZEB1_HUMAN.H11MO.0.A + MITF_HUMAN.H11MO.0.A + RXRB_HUMAN.H11MO.0.C + SNAI1_HUMAN.H11MO.0.C + IKZF1_HUMAN.H11MO.0.C + SOX5_HUMAN.H11MO.0.C + GC_content+ atac_count+ in_gene+ n_enhancers, data = all_without_extrinsic, method = "glm", trControl = ctrl)

#view summary of LOOCV               
print(all_without_extrinsic_model)
print(summary(all_without_extrinsic_model))





## Extrinsic alone

m1 <- glm(top_bottom_mean ~ 1 + s_prop+ g2_prop+ cd24_prop, family = "binomial", data = all_with_extrinsic)
print(summary(m1))
sum(round(m1$fitted.values) == all_with_extrinsic$top_bottom_mean)/nrow(all_with_extrinsic)
#plot_coeffs_S(m1, "mean_coeffs_extrinsic_features.pdf", title = "mean: extrinsic")

## Cross validation on Extrinsic only model

#Cross validation on Extrinsic only model

#specify the cross-validation method
ctrl <- trainControl(method = "LOOCV")

all_with_extrinsic$mean_level[all_with_extrinsic$top_bottom_mean==1] <- "High Mean"
all_with_extrinsic$mean_level[all_with_extrinsic$top_bottom_mean==0] <- "Low Mean"
all_with_extrinsic$mean_level <- factor(all_with_extrinsic$mean_level, levels = c("Low Mean", "High Mean"))

#fit a regression model and use LOOCV to evaluate performance
all_only_extrinsic_model <- train(mean_level ~ 1 + s_prop+ g2_prop+ cd24_prop, data = all_with_extrinsic, method = "glm", trControl = ctrl)

#view summary of LOOCV               
print(all_only_extrinsic_model)
print(summary(all_only_extrinsic_model))



## Intrinsic plus extrinsic

m2 <- glm(top_bottom_mean ~ 1 + s_prop + g2_prop + cd24_prop + H3K27me3 + H3K4me2 +  ZNF76_HUMAN.H11MO.0.C +  E2F7_HUMAN.H11MO.0.B +  BACH1_HUMAN.H11MO.0.A +  SMAD3_HUMAN.H11MO.0.B +  E2F3_HUMAN.H11MO.0.A +  SOX5_HUMAN.H11MO.0.C +  atac_count, family = "binomial", data = all_with_extrinsic)
print(summary(m2))
sum(round(m2$fitted.values) == all_with_extrinsic$top_bottom_mean)/nrow(all_with_extrinsic)
plot_coeffs_S(m2, "SuppFigure5F.pdf", title = "mean: intrinsic + extrinsic")

## Cross validation on Intrinsic + Extrinsic  model

#Cross validation on Intrinsic + Extrinsic  model

#specify the cross-validation method
ctrl <- trainControl(method = "LOOCV")

#fit a regression model and use LOOCV to evaluate performance
all_intrinsic_plus_extrinsic_model <- train(mean_level ~ 1 + s_prop + g2_prop + cd24_prop + H3K27me3 + H3K4me2 +  ZNF76_HUMAN.H11MO.0.C +  E2F7_HUMAN.H11MO.0.B +  BACH1_HUMAN.H11MO.0.A +  SMAD3_HUMAN.H11MO.0.B +  E2F3_HUMAN.H11MO.0.A +  SOX5_HUMAN.H11MO.0.C +  atac_count, data = all_with_extrinsic, method = "glm", trControl = ctrl)

#view summary of LOOCV               
print(all_intrinsic_plus_extrinsic_model)
print(summary(all_intrinsic_plus_extrinsic_model))

mean_accuracies <- read.table("dat/mean_accuracies.txt", head = T)
ggplot(data=mean_accuracies, aes(x=Category, y=Accuracy)) + geom_bar(stat="identity") + ylim(0, 1) + xlab("") + ggtitle("Mean")
ggsave("SuppFigure5E.pdf")

