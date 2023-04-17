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
t1$MIN <- m1$residuals
t1$twopower_MIN <- 2^m1$residuals
t1$fano <- t1$var/t1$mean
t1$cv2 <- t1$var/t1$mean^2

t1$pool <- args[2]
t1$mean_z <- (t1$mean - mean(t1$mean))/sd(t1$mean)
t1$var_z <- (t1$var - mean(t1$var))/sd(t1$var)
t1$fano_z <- (t1$fano - mean(t1$fano))/sd(t1$fano)
t1$cv2_z <- (t1$cv2 - mean(t1$cv2))/sd(t1$cv2)

pdf("residuals.pdf")
plot(m1$residuals, log2(t1$mean))
plot(hist(t1$mean_z))
plot(hist(t1$var_z))
plot(t1$mean_z, t1$var_z)
plot(t1$mean_z, t1$fano_z)
plot(t1$mean_z, t1$cv2_z)
dev.off()


write.table(t1, file = args[3], sep = "\t", quote = F, row.names = F)
