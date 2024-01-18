require(ggplot2)
reps <- read.table("dat/SL4_denoise_comparereplicates_10_5_2_10_merged.tsv", head = T, sep = "\t")

pdf("fig1_1.pdf", useDingbats=FALSE)
    ggplot(reps) + geom_point(aes(x = log10(mean_x), y = log10(mean_y))) +  xlab("Replicate 1 mean (log10)") +  ylab("Replicate 2 mean (log10)")
    #abline(lm(reps$mean_y ~ reps$mean_x))
#ggsave("101921run_rep1_rep2_mean.pdf")

     ggplot(reps) + geom_point(aes(x = log10(var_x), y = log10(var_y))) +  xlab("Replicate 1 variance (log10)") +  ylab("Replicate 2 variance (log10)")
    #abline(lm(log10(reps$var_y) ~ log10(reps$var_x)))
dev.off()
#ggsave("101921run_rep1_rep2_variance.pdf")
