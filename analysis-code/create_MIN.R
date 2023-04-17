args <- commandArgs(trailing = T)
if (length(args) < 2) {
    print("Usage Rscript create_MIN.R input op")
    stop()
}
t1 <- read.table(args[1], head = T, sep = "\t")
t1 <- t1[t1$var!= 0, ]
print(head(t1[, 1:4]))
print(t1[t1$var ==0, ])
m1 <- lm(log2(t1$var) ~ log2(t1$mean))
t1$pool <- args[2]
t1$MIN <- m1$residuals
plot(m1$residuals, log2(t1$mean))
dev.off()


write.table(t1, file = args[3], sep = "\t", quote = F, row.names = F)
