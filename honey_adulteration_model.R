#mengatur direktori kerja
setwd("D:\\GitHub\\Madu")

#import excel
library(readxl)
df1 <- as.data.frame(read_excel("ftirhoney.xlsx", sheet="Lembar1"))
View(df1)

# visualize raw spectra
x1 <- df1[,-c(1,2,3)]
y1 <- df1[,c(1,2,3)]
df2 <- data.frame(y1, I(x1))
View(df2)
wavenumber <- read_excel("wavenumber.xlsx")
matplot(wavenumber, t(df2$x[c(20,1,25),]),
        type='l', ylab="absorbance (a.u)",
        xlab="wavenumber [1/cm]",
        xlim=c(max(wavenumber), min(wavenumber)),
        ylim=c(-.1,.8), col=c(1,1,1),
        lty=c(1,2,3)
)
legend("topleft", c("honey", "cane sugar", "corn syrup"),
       col=c(1,1,1), lty=c(1,2,3), bty="n")

# smoothing
library(prospectr)
sgvec <- savitzkyGolay(X = df2$x1[25,], p=3, w=19, m=0)
matplot(wavenumber[c(19:2262),], t(sgvec), type='l',
        ylab="absorbance (a.u)",
        xlab="wavenumber [1/cm]",
        xlim=c(max(wavenumber), min(wavenumber)))
sg <- savitzkyGolay(X = df2$x1, p=3, w=19, m=0)
df3 <- data.frame(y1, I(sg))
View(df3)
par(mar=c(4.5,4,1,1))
wavenumber2 <- wavenumber[c(19:2262),]
matplot(wavenumber2, t(df3$sg[c(20,1,25),]),
        type='l',ylab="absorbance (a.u)",
        xlab="wavenumber [1/cm]",
        xlim=c(max(wavenumber2), min(wavenumber2)),
        col=c(1,1,1),
        lty=c(1,2,3)
        )
legend("topleft", c("honey", "cane sugar", "corn syrup"),
       col=c(1,1,1), lty=c(1,2,3), bty="n")

# wavenumber selection for cane sugar

matplot(wavenumber2, t(df3$sg[c(20,1),]),
        type='l',ylab="absorbance (a.u)",
        xlab="wavenumber [1/cm]",
        xlim=c(max(wavenumber2), min(wavenumber2)),
        col=c(1,1,1),
        lty=c(1,2,3)
)
legend("topleft", c("honey", "cane sugar"),
       col=c(1,1), lty=c(1,2), bty="n", cex=.8)

abline(v=c(2951, 2901, 1431, 1383, 1360, 1269, 1161, 1070, 1040,
           885, 835, 791), col="gray", lty=2)

# wavenumber selection for corn syrup
matplot(wavenumber[c(19:2262),], t(df3$sg[c(20,25),]), 
        type='l', ylab='absorbance [a.u]', 
        xlab="wavelength [1/cm]",
        xlim=c(max(wavenumber2), min(wavenumber2)),
        col=c(1,1),
        lty=c(1,2))
legend("topleft", c("honey", "corn syrup"),
       col=c(1,1),lty=c(1,2),
       bty="n")
abline(v=c(3275,2900, 2950, 1260, 1355, 1435, 835, 790), col='gray', lty=2)

# principal component analysis -- without unkown samples data
x.pca <- as.data.frame(df3$sg[1:58, c(226,249,275,355,371,418,474,521, 533,
                                      558, 1320,1346, 560,519,469,1514)])
pca_data <- data.frame(df2$group[1:58], x.pca)
colnames(pca_data)[colnames(pca_data) == 'df2.group.1.58.'] <- 'group'
View(pca_data)

library(mixOmics)


tune.pca.multi <- tune.pca(as.matrix(x.pca[c(1:43),]), ncomp = 10, scale = TRUE)
plot(tune.pca.multi)
final.pca.multi <- pca(as.matrix(x.pca[c(1:43),]), ncomp = 2, 
                       center = TRUE, scale = TRUE)# final.pca.multi  
# Lists possible outputs
final.pca.multi$var.tot

plotIndiv(final.pca.multi,
          comp = c(1, 2),   # Specify components to plot
          ind.names = FALSE, # Show row names of samples
          group = pca_data$group[c(1:43)],
          title = 'Honey, PCA comp 1 - 2',
          legend = TRUE, legend.title = 'sample')

#Partial Least Square Discriminant Analysis
datatest <- pca_data[c(1, 3, 6, 9, 12, 15, 21, 23, 27, 29, 31, 38, 42),]
datatrain <- pca_data[-c(1, 3, 6, 9, 12, 15, 21, 23, 27, 29, 31, 38, 42),]
View(datatest)
plsda.honey <- plsda(as.matrix(datatrain$x),datatrain$group, 
                     ncomp = 10)
set.seed(090924) 
perf.plsda.honey <- perf(plsda.honey, validation = 'Mfold', folds = 3, 
                         progressBar = FALSE,  
                         nrepeat = 100)        
plot(perf.plsda.honey, sd = TRUE)

final.plsda.honey <- plsda(as.matrix(datatrain$x),datatrain$group,
                           ncomp = 2)
plotIndiv(final.plsda.honey, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = F, 
          title = 'PLS-DA on Honey comp 1-2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')
background.max <- background.predict(final.plsda.honey, 
                                     comp.predicted = 2,
                                     dist = 'mahalanobis.dist')
plotIndiv(final.plsda.honey, comp = 1:2, group = datatrain$group,
          ind.names = F, title = 'Mahalanobis distance',
          legend = TRUE,  background = background.max)
auc.srbct <- auroc(final.plsda.honey, newdata = datatest$x, 
                   outcome.test = as.factor(datatest$group))

predict.plsda.honey <- predict(final.plsda.honey, as.matrix(datatrain$x),
                               dist = "mahalanobis.dist")

dataprediction <- data.frame(predict.plsda.honey$class, Truth = datatrain$group)
confusion_matrix <- table(true=dataprediction$Truth,
                          predicted = dataprediction$mahalanobis.dist.comp2)
print(confusion_matrix)
write.csv(confusion_matrix, file="confusiontrain.csv")

# external validation plsda
predict.plsda.honey <- predict(final.plsda.honey, as.matrix(datatest$x),
                               dist = "mahalanobis.dist")

dataprediction <- data.frame(predict.plsda.honey$class, Truth = datatest$group)
confusion_matrix <- table(true=dataprediction$Truth,
                          predicted = dataprediction$mahalanobis.dist.comp2)
print(confusion_matrix)
write.csv(confusion_matrix, file="confusiontest.csv")

# predcition of unknown samples
x.sampel <- as.data.frame(df3$sg[44:58, c(226,249,275,355,371,418,474,521, 533,
                                      558, 1320,1346, 560,519,469,1514)])
sample_data <- data.frame(df2$group[44:58], x.sampel)
colnames(sample_data)[colnames(sample_data) == 'df2.group.44.58.'] <- 'group'
View(sample_data)
predict.plsda.honey <- predict(final.plsda.honey, as.matrix(sample_data$x),
                               dist = "mahalanobis.dist")

dataprediction <- data.frame(predict.plsda.honey$class, 
                             Truth = sample_data$group)
print(dataprediction)


