library(survival)
setwd("~/projects/hits-clip/data/20130322/")
df = read.table("fixed_matrix.txt", header=TRUE, row.names=1)

attach(df)
erpos = df[df$ER==1,]
erneg = df[df$ER==1,]

plotsurv <- function(fit, miR, ER, lty, names, p){
    pdf(paste0(miR, "__", ER, ".pdf"), width=6, height=8)
    plot(fit, lty=lty, mark.time=TRUE)
    legend(.1, .1, names, lty=lty)
    title(main=paste0(miR, " (ER ", ER, ")"), sub=paste('p =', p), xlab="Time to Relapse")
    dev.off()
}

# redoing stuff in order to get significant p-value
# setwd("~/projects/hits-clip/data/20130328/")
# df = read.table("mir9-erpos.txt", header=TRUE, row.names=1)
# 
# mir = "hsa-miR-9-5p__logrank"
# fit = survfit(Surv(Time, Relapse)~miR.9.5, data=df, rho=0)
# sdf = survdiff(Surv(Time, Relapse)~miR.9.5, data=df, rho=0)
# sdf
# p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)
# 
# fit = survfit(Surv(Time, Relapse)~miR.9.5, data=df)
# sdf = survdiff(Surv(Time, Relapse)~miR.9.5, data=df)
# sdf
# p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# ILMN_3167194 = hsa-miR-9-3p
mir = "hsa-miR-9-3p"
fit = survfit(Surv(Time, Relapse)~ILMN_3167194, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167194, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3167194, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167194, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# ILMN_3167447 = hsa-miR-9-5p
mir = "hsa-miR-9-5p"
fit = survfit(Surv(Time, Relapse)~ILMN_3167447, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3167447, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# again but with log-rank
mir = "hsa-miR-9-5p__petopeto"
fit = survfit(Surv(Time, Relapse)~ILMN_3167447, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447, data=erpos, rho=1)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3167447, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# ILMN_3167944 = hsa-miR-193b-3p
mir = "hsa-miR-193b-3p"
fit = survfit(Surv(Time, Relapse)~ILMN_3167944, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167944, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3167944, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167944, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# ILMN_3168366 = hsa-miR-193a-3p
mir = "hsa-miR-193a-3p"
fit = survfit(Surv(Time, Relapse)~ILMN_3168366, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3168366, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3168366, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3168366, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# ILMN_3167681 = hsa-miR-221-3p
mir = "hsa-miR-221-3p"
fit = survfit(Surv(Time, Relapse)~ILMN_3167681, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167681, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3167681, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167681, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# ILMN_3168429 = hsa-miR-34a-5p
mir = "hsa-miR-34a-5p"
fit = survfit(Surv(Time, Relapse)~ILMN_3168429, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3168429, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

fit = survfit(Surv(Time, Relapse)~ILMN_3168429, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3168429, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:2, c("low", "high"), p)

# 2. mir9-5p together with 193a-3p
mir = "hsa-miR-9-5p_hsa-miR-193a-3p"
names=c("low-low", "low-high", "high-low", "high-high")
fit = survfit(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168366, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168366, data=erpos, rho=1)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:4, names, p)

# trying to plot only h-h and l-l
fit
plot(fit, lty=1:4)
legend('bottomleft', names, lty=1:4)
title(main=paste0(miR, " (ER ", ER, ")"), sub=paste('p =', p), xlab="Time to Relapse")



fit = survfit(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168366, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168366, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:4, names, p)

# mir9-5/and mir17-92 cluster (18a and 19a)
# ILMN_3167447 = hsa-miR-9-5p
# ILMN_3168282 = hsa-miR-18a-5p
### ILMN_3167178 = hsa-miR-18a-3p
# ILMN_3167787 hsa-miR-19a-3p
mir = "hsa-miR-9-5p_hsa-miR-18a-5p_hsa-miR-19a-3p"
names = c("l-l-l", "l-l-h", "l-h-l", "l-h-h", "h-l-l", "h-l-h", "h-h-l", "h-h-h")
fit = survfit(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168282 + ILMN_3167787, data=erpos)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168282 + ILMN_3167787, data=erpos)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:8, names, p)

fit = survfit(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168282 + ILMN_3167787, data=erneg)
sdf = survdiff(Surv(Time, Relapse)~ILMN_3167447 + ILMN_3168282 + ILMN_3167787, data=erneg)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "neg", lty=1:8, names, p)


# 20130328
setwd("~/projects/hits-clip/data/20130328/")
df = read.table("ditzel_highlow.txt", header=TRUE, row.names=1, na.strings="", sep="\t")

attach(df)

# plotsurv <- function(fit, miR, lty, names, p){
#     pdf(paste0(miR, ".pdf"), width=6, height=8)
#     plot(fit, lty=lty, mark.time=TRUE)
#     legend(.1, .1, names, lty=lty)
#     title(main=miR, sub=paste('p =', p), xlab="Time to Relapse")
#     dev.off()
# }
# 
# mir <- "Ditzel__hsa-miR-9-5p"
# names <- c("low", "high")
# fit <- survfit(Surv(timerec, relapse)~hsa.miR.9.5p, data=df)
# sdf <- survdiff(Surv(timerec, relapse)~hsa.miR.9.5p, data=df)
# sdf
# p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# plotsurv(fit, mir, lty=1:2, names, p)

df = read.table("combined_erpos.txt", header=TRUE, row.names=1, sep="\t")
attach(df)

plotsurv <- function(fit, miR, ER, lty, names, p){
    pdf(paste0(miR, "__", ER, ".pdf"), width=6, height=8)
    plot(fit, lty=lty, mark.time=TRUE)
    legend(.1, .1, names, lty=lty)
    title(main=paste0(miR, " (ER ", ER, ")"), sub=paste('p =', p), xlab="Time to Relapse")
    dev.off()
}

mir = "combined_all_hsa-miR-9-5p"
fit = survfit(Surv(Time, Relapse)~hsa.miR.9.5p, data=df)
sdf = survdiff(Surv(Time, Relapse)~hsa.miR.9.5p, data=df)
sdf
p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
plotsurv(fit, mir, "pos", lty=1:2, c("low", "high"), p)

# playing with binomial test for genes that are regulated by miR
fail <- 53 + (105 - 51)
success <- (105 - 53) + 51

binom.test(c(107,103), p=0.5)
