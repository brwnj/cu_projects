library(survival)
setwd("~/projects/hits-clip/data/20130404/")
df = read.table("GSE22216_highlow.txt", header=TRUE, row.names=1)

attach(df)
# erpos = df[df$ER==1,]
# erneg = df[df$ER==1,]

plotsurv <- function(fit, miR, lty, names, p){
    pdf(paste0(miR, ".pdf"), width=6, height=8)
    plot(fit, lty=lty, mark.time=TRUE)
    legend('bottomleft', names, lty=lty)
    title(main=miR, sub=paste('p =', p), xlab="Time to Relapse")
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


# 20130404 working on new miRs and clusters
library(survival)
setwd("~/projects/hits-clip/data/20130404/")
df = read.table("GSE22216_highlow.txt", header=TRUE, row.names=1)
attach(df)

# plot(fit[1], col="blue", conf.int=FALSE)
# lines(fit[2], col="red")
# colors=c("blue", "red")
# names=c("low","high")
# 
# mir = "hsa-miR-9-5p"
# fit = survfit(Surv(Time, Relapse)~hsa.miR.9.5p)
# sdf = survdiff(Surv(Time, Relapse)~hsa.miR.9.5p, rho=1)
# sdf
# p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# pdf(paste0(mir, ".pdf"), width=6, height=8)
# plot(fit, col=colors)
# legend('bottomleft', names, col=colors, lty=c(1,1), lwd=c(2,2))
# title(main=mir, sub=paste('p =', p), xlab="Time")
# dev.off()

#single mir plot
single_mir <- function(df, mir){
    colors=c("blue", "red")
    names=c("low","high")
    fit = survfit(as.formula(paste0("Surv(Time, Relapse)~",mir)))
    sdf = survdiff(as.formula(paste0("Surv(Time, Relapse)~",mir)), rho=1)
    # sdf
    p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    pdf(paste0(mir, ".pdf"), width=6, height=8)
    plot(fit, col=colors)
    legend('bottomleft', names, col=colors, lty=c(1,1), lwd=c(2,2))
    title(main=mir, sub=paste('p =', p), xlab="Time")
    dev.off()
}

single_mir(df, "hsa.miR.9.5p")
single_mir(df, "hsa.miR.9.3p")
single_mir(df, "hsa.miR.193b.3p")
single_mir(df, "hsa.miR.193a.3p")
single_mir(df, "hsa.miR.221.3p")
single_mir(df, "hsa.miR.34a.5p")
single_mir(df, "hsa.miR.18a.5p")
single_mir(df, "hsa.miR.19b.3p")

many_mirs <- function(df, mirs){
    colors <- c("blue", "red")
    names <- c("all low", "all high")
    fit <- survfit(as.formula(paste("Surv(Time, Relapse)~",paste(mirs, collapse="+"))))
    sdf <- survdiff(as.formula(paste("Surv(Time, Relapse)~",paste(mirs, collapse="+"))), rho=0)
    p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    pdf(paste0(paste(mirs, collapse="_"), ".pdf"), width=6, height=8)
    plot(fit[1], col="blue", conf.int=FALSE)
    lines(fit[length(fit$n)], col="red")
    legend('bottomleft', names, col=colors, lty=c(1,1), lwd=c(2,2))
    title(main=paste(mirs, collapse=" "), sub=paste('p =', p), xlab="Time")
    dev.off()
}

many_mirs(df, c("hsa.miR.18a.5p", "hsa.miR.19b.3p"))
many_mirs(df, c("hsa.miR.9.5p", "hsa.miR.18a.5p", "hsa.miR.19b.3p"))
many_mirs(df, c("hsa.miR.9.5p", "hsa.miR.193a.3p"))
many_mirs(df, c("hsa.miR.9.5p", "hsa.miR.18a.5p"))
many_mirs(df, c("hsa.miR.9.5p", "hsa.miR.19b.3p"))
many_mirs(df, c("hsa.miR.9.5p", "hsa.miR.221.3p"))
many_mirs(df, c("hsa.miR.9.5p", "hsa.miR.221.3p", "hsa.miR.18a.5p", "hsa.miR.19b.3p"))

xlabs = "Time"
ylabs = "Survival Probability"
xlims = c(0,max(fit$time))
ylims = c(0,1)

times <- seq(0, max(fit$time), by = 1)
subs1 <- 1:length(levels(summary(fit)$strata))
subs2 <- 1:length(summary(fit,censored=T)$strata)
subs3 <- 1:length(summary(fit,times = times,extend = TRUE)$strata)

gdf <- data.frame(
    time = fit$time[subs2],
    n.risk = fit$n.risk[subs2],
    n.event = fit$n.event[subs2],
    surv = fit$surv[subs2],
    strata = factor(summary(fit, censored=TRUE)$strata[subs2]),
    upper = fit$upper[subs2],
    lower = fit$lower[subs2]
)

p <- ggplot(gdf, aes(time, surv)) +
    geom_step(aes(color=strata)) +
    scale_color_manual(values=c("blue", "red"))
#     scale_color_hue(breaks=levels(gdf$strata), labels=c("Low", "High")) + 
    xlab("Time") +
    ylab("Survival Probability") + 
    ggtitle(mir) + 
    theme_bw()
p