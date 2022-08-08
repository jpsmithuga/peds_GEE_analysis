###########################################################
#' Figures for unadjusted and adjusted models: 
#'      "Clinical and radiologic factors associated with 
#'      detection of Mycobacterium tuberculosis in 
#'      children under 5 years old using invasive and 
#'      noninvasive sample collection techniques - Kenya"
#'
#' Author: Jonathan Smith, PhD, MPH
#'         Global Tuberculosis Branch
#'         Division of Global HIV and Tuberculosis
#'         Centers for Disease Control and Prevention
###########################################################
# rm(list = ls())

# READ IN MODEL RESULTS
#' Change filepath to where data and code are 
#' saved from github. Default to working directory.
flpth <- getwd() 
source(paste0(flpth, "/Models - unadjusted and adjusted analysis.R"))


###################################################
## Unadjusted Figure
###################################################
unadj_figdata_toto <- data.frame(unadj_logit_table, unadj_gee_table)
unadj_figdata_toto$factor <- c("Female","No or mild\nimmunosuppression","Moderate or high\nimmunosuppression","1 to <2 Years", "2-5 Years", 
                               "Malnutrition","Positive TST or IGRA","History of Exposure",
                               "Prolonged Cough", "Prolonged Fever", "Prolonged Lethargy",
                               "Two Symptoms","Three Symptoms",
                               "CXR Consistent with TB", "IGRA Only","TST Only")
unadj_figdata_toto$any.txt <- paste0(sprintf("%.1f", round(unadj_figdata_toto[,"any.est"],1)), 
                                     " (", sprintf("%.1f", round(unadj_figdata_toto[,"any.lwr"],1)),"-",
                                     sprintf("%.1f", round(unadj_figdata_toto[,"any.upr"],1)),")")
unadj_figdata_toto$gsasp.txt <- paste0(sprintf("%.1f", round(unadj_figdata_toto[,"gsasp.est"],1)), 
                                       " (", sprintf("%.1f", round(unadj_figdata_toto[,"gsasp.lwr"],1)),"-",
                                       sprintf("%.1f", round(unadj_figdata_toto[,"gsasp.upr"],1)),")")
unadj_figdata_toto$inspt.txt <- paste0(sprintf("%.1f", round(unadj_figdata_toto[,"inspt.est"],1)), 
                                       " (", sprintf("%.1f", round(unadj_figdata_toto[,"inspt.lwr"],1)),"-",
                                       sprintf("%.1f", round(unadj_figdata_toto[,"inspt.upr"],1)),")")
unadj_figdata_toto$npasp.txt <- paste0(sprintf("%.1f", round(unadj_figdata_toto[,"npasp.est"],1)), 
                                       " (", sprintf("%.1f", round(unadj_figdata_toto[,"npasp.lwr"],1)),"-",
                                       sprintf("%.1f", round(unadj_figdata_toto[,"npasp.upr"],1)),")")
unadj_figdata_toto$stool.txt <- paste0(sprintf("%.1f", round(unadj_figdata_toto[,"stool.est"],1)), 
                                       " (", sprintf("%.1f", round(unadj_figdata_toto[,"stool.lwr"],1)),"-",
                                       sprintf("%.1f", round(unadj_figdata_toto[,"stool.upr"],1)),")")
unadj_figdata_toto$strng.txt <- paste0(sprintf("%.1f", round(unadj_figdata_toto[,"strng.est"],1)), 
                                       " (", sprintf("%.1f", round(unadj_figdata_toto[,"strng.lwr"],1)),"-",
                                       sprintf("%.1f", round(unadj_figdata_toto[,"strng.upr"],1)),")")

# reorder
unadj_figdata_toto <- unadj_figdata_toto[c(1:5,7,15:16,8,14,6,9:13),]

unadj_figdata_toto$inspt.txt <- ifelse(is.na(unadj_figdata_toto$inspt.est), "", unadj_figdata_toto$inspt.txt)

# remove number of symptoms
nn <- nrow(unadj_figdata_toto)
unadj_figdata_toto <- unadj_figdata_toto[-c(nn, nn - 1), ]

seqx <- 2^c(-3:6)
unadj_figdata_toto[grep("lwr", names(unadj_figdata_toto))] <- lapply(unadj_figdata_toto[grep("lwr", names(unadj_figdata_toto))],
                                                                     function(x) ifelse(x < min(seqx), min(seqx), x))
unadj_figdata_toto[grep("upr", names(unadj_figdata_toto))] <- lapply(unadj_figdata_toto[grep("upr", names(unadj_figdata_toto))],
                                                                     function(x) ifelse(x > max(seqx), max(seqx), x))
## setup figure
xx <- 0.0005
maxx <- 1e10
ymax <- nrow(unadj_figdata_toto)
yy <- ymax:1
yh <- seq(-0.2, 0.2, length.out = 6)
cex.adj <- 1.0
poss <- 3
cex.text <- 0.85
cix <- seq(ceiling(log(max(seqx)))+1, floor(log(maxx)), length.out = 6)
madj <- -5
cols <- c("#BF457E", "#5E9ABF", "#638C5D", "#F2B33D","#902423")

# Global to save figures or not. Will print PDF if TRUE; default FALSE
savefigs <- FALSE

wdth <- 12
hght <- 8
if(savefigs){
  pdf(paste0("OR_table_clinical_unadjusted_all_", Sys.Date(), ".pdf"),
     width = wdth, height = hght)
 }

#set up blank plot and axes
plot(log(range(xx, maxx)), range(0.9, ymax+0.1), type='n', axes = FALSE, xlab='', ylab='')
axis(side = 1, at = log(seqx), labels = seqx, cex.axis = cex.text, mgp = c(2, 1.5, 1))
mtext(side = 1, at = log(median(seqx)), 'Odds Ratio (log)', line = 2.5, cex = cex.text, padj = 1)
axis(side = 3, at = log(seqx), labels = seqx, cex.axis = cex.text, mgp = c(2, 1.5, 1))
mtext(side = 3, at = log(median(seqx)), 'Odds Ratio (log)', line = 3.5, cex = cex.text, padj = 1)
par(xpd = TRUE)
segments(x0 = log(1), y0 = ymax+0.5, y1 = 0.5, lty = 3, lwd = 0.5)

# Add estimates
for (i in 1:nrow(unadj_figdata_toto)){
  segments(log(unadj_figdata_toto[i, "strng.lwr"]), yy[i]+yh[1], log(unadj_figdata_toto[i,"strng.upr"]), yy[i]+yh[1])
  segments(log(unadj_figdata_toto[i, "stool.lwr"]), yy[i]+yh[2], log(unadj_figdata_toto[i,"stool.upr"]), yy[i]+yh[2])
  segments(log(unadj_figdata_toto[i, "npasp.lwr"]), yy[i]+yh[3], log(unadj_figdata_toto[i,"npasp.upr"]), yy[i]+yh[3])
  segments(log(unadj_figdata_toto[i, "inspt.lwr"]), yy[i]+yh[4], log(unadj_figdata_toto[i,"inspt.upr"]), yy[i]+yh[4])
  segments(log(unadj_figdata_toto[i, "gsasp.lwr"]), yy[i]+yh[5], log(unadj_figdata_toto[i,"gsasp.upr"]), yy[i]+yh[5])
  segments(log(unadj_figdata_toto[i, "any.lwr"]), yy[i]+yh[6]+0.1, log(unadj_figdata_toto[i,"any.upr"]), yy[i]+yh[6]+0.1, lty = 3)
  
  points(log(unadj_figdata_toto[i,"strng.est"]), yy[i]+yh[1], pch = 22, bg = cols[5], cex = cex.adj)
  points(log(unadj_figdata_toto[i,"stool.est"]), yy[i]+yh[2], pch = 22, bg = cols[4], cex = cex.adj)
  points(log(unadj_figdata_toto[i,"npasp.est"]), yy[i]+yh[3], pch = 22, bg = cols[3], cex = cex.adj)
  points(log(unadj_figdata_toto[i,"inspt.est"]), yy[i]+yh[4], pch = 22, bg = cols[2], cex = cex.adj)
  points(log(unadj_figdata_toto[i,"gsasp.est"]), yy[i]+yh[5], pch = 22, bg = cols[1], cex = cex.adj)
  points(log(unadj_figdata_toto[i,"any.est"]), yy[i]+yh[6]+0.1, pch = 25, bg = "black", cex = cex.adj/1.5)
  
  text(xx-2.5, yy[i]+mean(yh[1:6]), unadj_figdata_toto[i,"factor"], cex = cex.text, pos = 2)
  
  hh <- 0.25
  text(cix[1], yy[i]+mean(yh[1:5])-hh, unadj_figdata_toto[i,"any.txt"], cex = cex.text, pos = poss)
  text(cix[2], yy[i]+mean(yh[1:5])-hh, unadj_figdata_toto[i,"gsasp.txt"], cex = cex.text, pos = poss)
  text(cix[3], yy[i]+mean(yh[1:5])-hh, unadj_figdata_toto[i,"inspt.txt"], cex = cex.text, pos = poss)
  text(cix[4], yy[i]+mean(yh[1:5])-hh, unadj_figdata_toto[i,"npasp.txt"], cex = cex.text, pos = poss)
  text(cix[5], yy[i]+mean(yh[1:5])-hh, unadj_figdata_toto[i,"stool.txt"], cex = cex.text, pos = poss)
  text(cix[6], yy[i]+mean(yh[1:5])-hh, unadj_figdata_toto[i,"strng.txt"], cex = cex.text, pos = poss)
}  

# Add additional text/comments
coladj <- 0.75
text(x = mean(cix), y = ymax + 2, "Odds Ratio (95% Confidence Interval)", font = 2, cex = cex.text)

points(cix, rep(ymax + coladj-0.25, length(cix)), pch = c(25, rep(22,5)), bg = c("black",cols), cex = c(0.80, rep(1.5,5))) 
text(cix[1], ymax + coladj, "Any Positive\nSpecimen", cex = cex.text, pos = poss)
text(cix[2], ymax + coladj, "Gastric\nAspirate", cex = cex.text, pos = poss)
text(cix[3], ymax + coladj, "Induced\nSputum", cex = cex.text, pos = poss)
text(cix[4], ymax + coladj, "Nasopharyngeal\nAspirate", cex = cex.text, pos = poss)
text(cix[5], ymax + coladj, "\nStool", cex = cex.text, pos = poss)
text(cix[6], ymax + coladj, "\nString", cex = cex.text, pos = poss)

# HIV
text(x = -9.5,  y = yy[3] + 0.5, "HIV Positive\n(Reference:\nHIV negative)", cex = cex.text)
segments(x0 = -7.8, y0 = yy[4] + 0.75, y1 = yy[2] + 0.25, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[4] + 0.75, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[2] + 0.25, lty = 1)

# Age
text(x = -9.5,  y = yy[5] + 0.5, "Age\n(Reference:\n<1 Year)", cex = cex.text)
segments(x0 = -7.8, y0 = yy[6] + 0.75, y1 = yy[4] + 0.25, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[6] + 0.75, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[4] + 0.25, lty = 1)
#  text(x = log(0.0006), y = 7.5, "Age", srt = 90, cex = cex.text)

# IE/TST/IGRA
segments(x0 = -6.2, y0 = yy[6] - 0.50, y1 = yy[8])
arrows(x0 = -6.2,  x1 = -5.2, y0 = yy[7], length = 0.07)
arrows(x0 = -6.2,  x1 = -5.2, y0 = yy[8], length = 0.07)

# TB symptoms
text(x = -9.5,  y = 2.5, "TB-like\nSymptoms", cex = cex.text)
segments(x0 = -7.8, y0 = yy[11] + 0.25, y1 = 0.75, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[11] + 0.25, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = 0.75, lty = 1)

# Legend
legend(x = 8, y = 0.25, c("Per-Patient Analysis (Standard Logistic Regression)","Per-Specimen Analysis (GEE Logistic Models)", "95% Confidence Interval"),
       pch = c(25, 22, NA), pt.bg = c("black", "darkgrey", "black"), lty = c(rep(NA, 2),1), bty = "n",
       cex = 0.75, pt.cex = rev(c(1, 0.75)))
if(savefigs) {
  dev.off()
}


###################################################
## Adjusted Figure
###################################################
adj_figdata_toto <- data.frame(adj_logit_table, adj_gee_table)
adj_figdata_toto$factor <- c("Female","No or mild\nimmunosuppression","Moderate or high\nimmunosuppression","1 to <2 Years", "2-5 Years", 
                             "Malnutrition","Positive TST or IGRA","History of Exposure",
                             "Prolonged Cough", "Prolonged Fever", "Prolonged Lethargy",
                             "Two Symptoms","Three Symptoms",
                             "CXR Consistent with TB", "IGRA Only","TST Only")
adj_figdata_toto$any.txt <- paste0(sprintf("%.1f", round(adj_figdata_toto[,"any.est"],1)), 
                                   " (", sprintf("%.1f", round(adj_figdata_toto[,"any.lwr"],1)),"-",
                                   sprintf("%.1f", round(adj_figdata_toto[,"any.upr"],1)),")")
adj_figdata_toto$gsasp.txt <- paste0(sprintf("%.1f", round(adj_figdata_toto[,"gsasp.est"],1)), 
                                     " (", sprintf("%.1f", round(adj_figdata_toto[,"gsasp.lwr"],1)),"-",
                                     sprintf("%.1f", round(adj_figdata_toto[,"gsasp.upr"],1)),")")
adj_figdata_toto$inspt.txt <- paste0(sprintf("%.1f", round(adj_figdata_toto[,"inspt.est"],1)), 
                                     " (", sprintf("%.1f", round(adj_figdata_toto[,"inspt.lwr"],1)),"-",
                                     sprintf("%.1f", round(adj_figdata_toto[,"inspt.upr"],1)),")")
adj_figdata_toto$npasp.txt <- paste0(sprintf("%.1f", round(adj_figdata_toto[,"npasp.est"],1)), 
                                     " (", sprintf("%.1f", round(adj_figdata_toto[,"npasp.lwr"],1)),"-",
                                     sprintf("%.1f", round(adj_figdata_toto[,"npasp.upr"],1)),")")
adj_figdata_toto$stool.txt <- paste0(sprintf("%.1f", round(adj_figdata_toto[,"stool.est"],1)), 
                                     " (", sprintf("%.1f", round(adj_figdata_toto[,"stool.lwr"],1)),"-",
                                     sprintf("%.1f", round(adj_figdata_toto[,"stool.upr"],1)),")")
adj_figdata_toto$strng.txt <- paste0(sprintf("%.1f", round(adj_figdata_toto[,"strng.est"],1)), 
                                     " (", sprintf("%.1f", round(adj_figdata_toto[,"strng.lwr"],1)),"-",
                                     sprintf("%.1f", round(adj_figdata_toto[,"strng.upr"],1)),")")

# reorder
adj_figdata_toto <- adj_figdata_toto[c(1:5,7,15:16,8,14,6,9:13),]

adj_figdata_toto$inspt.txt <- ifelse(is.na(adj_figdata_toto$inspt.est), "", adj_figdata_toto$inspt.txt)

# remove number of symptoms
adj_figdata_toto <- adj_figdata_toto[-c(nn, nn-1),]


adj_figdata_toto[grep("lwr", names(adj_figdata_toto))] <- lapply(adj_figdata_toto[grep("lwr", names(adj_figdata_toto))],
                                                                 function(x) ifelse(x < min(seqx), min(seqx), x))

adj_figdata_toto[grep("upr", names(adj_figdata_toto))] <- lapply(adj_figdata_toto[grep("upr", names(adj_figdata_toto))],
                                                                 function(x) ifelse(x > max(seqx), max(seqx), x))


if(savefigs){
  pdf(paste0("OR_table_clinical_unadjusted_all_", Sys.Date(), ".pdf"),
      width = wdth, height = hght)
}
plot(log(range(xx, maxx)), range(0.9, ymax+0.1), type='n', axes = FALSE, xlab='', ylab='')
axis(side = 1, at = log(seqx), labels = seqx, cex.axis = cex.text, mgp = c(2, 1.5, 1))
mtext(side = 1, at = log(median(seqx)), 'Odds Ratio (log)', line = 2.5, cex = cex.text, padj = 1)
axis(side = 3, at = log(seqx), labels = seqx, cex.axis = cex.text, mgp = c(2, 1.5, 1))
mtext(side = 3, at = log(median(seqx)), 'Odds Ratio (log)', line = 3.5, cex = cex.text, padj = 1)

par(xpd = TRUE)
segments(x0 = log(1), y0 = ymax+0.5, y1 = 0.5, lty = 3, lwd = 0.5)

for (i in 1:nrow(adj_figdata_toto)){
  segments(log(adj_figdata_toto[i, "strng.lwr"]), yy[i]+yh[1], log(adj_figdata_toto[i,"strng.upr"]), yy[i]+yh[1])
  segments(log(adj_figdata_toto[i, "stool.lwr"]), yy[i]+yh[2], log(adj_figdata_toto[i,"stool.upr"]), yy[i]+yh[2])
  segments(log(adj_figdata_toto[i, "npasp.lwr"]), yy[i]+yh[3], log(adj_figdata_toto[i,"npasp.upr"]), yy[i]+yh[3])
  segments(log(adj_figdata_toto[i, "inspt.lwr"]), yy[i]+yh[4], log(adj_figdata_toto[i,"inspt.upr"]), yy[i]+yh[4])
  segments(log(adj_figdata_toto[i, "gsasp.lwr"]), yy[i]+yh[5], log(adj_figdata_toto[i,"gsasp.upr"]), yy[i]+yh[5])
  segments(log(adj_figdata_toto[i, "any.lwr"]), yy[i]+yh[6]+0.1, log(adj_figdata_toto[i,"any.upr"]), yy[i]+yh[6]+0.1, lty = 3)
  
  points(log(adj_figdata_toto[i,"strng.est"]), yy[i]+yh[1], pch = 22, bg = cols[5], cex = cex.adj)
  points(log(adj_figdata_toto[i,"stool.est"]), yy[i]+yh[2], pch = 22, bg = cols[4], cex = cex.adj)
  points(log(adj_figdata_toto[i,"npasp.est"]), yy[i]+yh[3], pch = 22, bg = cols[3], cex = cex.adj)
  points(log(adj_figdata_toto[i,"inspt.est"]), yy[i]+yh[4], pch = 22, bg = cols[2], cex = cex.adj)
  points(log(adj_figdata_toto[i,"gsasp.est"]), yy[i]+yh[5], pch = 22, bg = cols[1], cex = cex.adj)
  points(log(adj_figdata_toto[i,"any.est"]), yy[i]+yh[6]+0.1, pch = 25, bg = "black", cex = cex.adj/1.5)
  
  text(xx-2.5, yy[i]+mean(yh[1:6]), adj_figdata_toto[i,"factor"], cex = cex.text, pos = 2)
  
  hh <- 0.25
  text(cix[1], yy[i]+mean(yh[1:5])-hh, adj_figdata_toto[i,"any.txt"], cex = cex.text, pos = poss)
  text(cix[2], yy[i]+mean(yh[1:5])-hh, adj_figdata_toto[i,"gsasp.txt"], cex = cex.text, pos = poss)
  text(cix[3], yy[i]+mean(yh[1:5])-hh, adj_figdata_toto[i,"inspt.txt"], cex = cex.text, pos = poss)
  text(cix[4], yy[i]+mean(yh[1:5])-hh, adj_figdata_toto[i,"npasp.txt"], cex = cex.text, pos = poss)
  text(cix[5], yy[i]+mean(yh[1:5])-hh, adj_figdata_toto[i,"stool.txt"], cex = cex.text, pos = poss)
  text(cix[6], yy[i]+mean(yh[1:5])-hh, adj_figdata_toto[i,"strng.txt"], cex = cex.text, pos = poss)
}  
coladj <- 0.75
text(x = mean(cix), y = ymax + 2.0, "Odds Ratio (95% Confidence Interval)", font = 2, cex = cex.text)

points(cix, rep(ymax + coladj-0.25, length(cix)), pch = c(25, rep(22,5)), bg = c("black",cols), cex = c(0.80, rep(1.5,5))) 
text(cix[1], ymax + coladj, "Any Positive\nSpecimen", cex = cex.text, pos = poss)
text(cix[2], ymax + coladj, "Gastric\nAspirate", cex = cex.text, pos = poss)
text(cix[3], ymax + coladj, "Induced\nSputum", cex = cex.text, pos = poss)
text(cix[4], ymax + coladj, "Nasopharyngeal\nAspirate", cex = cex.text, pos = poss)
text(cix[5], ymax + coladj, "\nStool", cex = cex.text, pos = poss)
text(cix[6], ymax + coladj, "\nString", cex = cex.text, pos = poss)

# HIV
text(x = -9.5,  y = yy[3] + 0.5, "HIV Positive\n(Reference:\nHIV negative)", cex = cex.text)
segments(x0 = -7.8, y0 = yy[4] + 0.75, y1 = yy[2] + 0.25, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[4] + 0.75, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[2] + 0.25, lty = 1)

# Age
text(x = -9.5,  y = yy[5] + 0.5, "Age\n(Reference:\n<1 Year)", cex = cex.text)
segments(x0 = -7.8, y0 = yy[6] + 0.75, y1 = yy[4] + 0.25, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[6] + 0.75, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[4] + 0.25, lty = 1)
#  text(x = log(0.0006), y = 7.5, "Age", srt = 90, cex = cex.text)

# IE/TST/IGRA
segments(x0 = -6.2, y0 = yy[6] - 0.50, y1 = yy[8])
arrows(x0 = -6.2,  x1 = -5.2, y0 = yy[7], length = 0.07)
arrows(x0 = -6.2,  x1 = -5.2, y0 = yy[8], length = 0.07)

# TB symtoms
text(x = -9.5,  y = 2.5, "TB-like\nSymptoms", cex = cex.text)
segments(x0 = -7.8, y0 = yy[11] + 0.25, y1 = 0.75, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = yy[11] + 0.25, lty = 1)
segments(x0 = -7.8, x1 = -7.6, y0 = 0.75, lty = 1)

legend(x = 8, y = 0.25, c("Per-Patient Analysis (Standard Logistic Regression)","Per-Specimen Analysis (GEE Logistic Models)", "95% Confidence Interval"),
       pch = c(25, 22, NA), pt.bg = c("black", "darkgrey", "black"), lty = c(rep(NA, 2),1), bty = "n",
       cex = 0.75, pt.cex = rev(c(1, 0.75)))

if(savefigs){dev.off()}