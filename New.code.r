library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)
library(seqinr)
library(ggplot2)
library(phylopath)
library(plyr)

dnds.alignment <- read.alignment("shrew.fas", format= "fasta")
kaks.1 <- kaks(dnds.alignment)
dn <- kaks.1$ka
ds <- kaks.1$ks
dn.ds <- dn/ds # dn.ds is assembled in Dataset.txt

dat <- read.csv("scaled.dat.csv")

primtree <- read.tree(file = "10k.tree.nwk", text = NULL, tree.names = NULL, skip = 0,
	comment.char = "#", keep.multi = FALSE)

tree.1 <- drop.tip(primtree, "Mirza_coquereli", trim.internal = TRUE)

tree <- drop.tip(tree.1, "Lepilemur_mustelinus", trim.internal = TRUE)

dat$Taxa<-gsub(" ", "_", dat$Taxa)
rownames(dat) <- dat$Taxa
dat$Taxa <- as.factor(dat$Taxa)


primate.bmr <- comparative.data(phy = tree, data = dat, 
names.col = Taxa, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#DAG Models
models <- list(
  one   = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ GC),
  two   = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ dNdS + GC),
  three = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BRNM),
  four  = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BM + BRNM),
  five  = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BM + BRNM + GC),
  six   = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM, BMR ~ BM),
  seven = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM, BMR ~ dNdS + BM),
  eight = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM),
  nine  = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM, BMR ~ dNdS),
  ten   = DAG(dNdS ~ BM + BRNM, BRNM ~ BM, GC ~ BRNM, BMR ~ dNdS + GC),
  eleven = DAG(BM ~ BRNM, dNdS ~ BRNM, GC ~ BM, BMR ~ BM),
  twelve = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ GC),
  thirteen = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ dNdS + GC),
  fourteen  = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ BRNM + BM),
  fifteen  = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ BRNM + BM + GC),
  sixteen   = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM, BMR ~ BRNM),
  seventeen = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM, BMR ~ dNdS + BRNM),
  eighteen = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM),
  nineteen  = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM, BMR ~ dNdS)
)




#DAG Models
models <- list(
  one   = DAG(dNdS ~ BM, BRNM ~ BM, GL ~ BRNM, BMR ~ GL),
  two   = DAG(dNdS ~ BM, BRNM ~ BM, GL ~ BRNM, BMR ~ dNdS + GL),
  three = DAG(dNdS ~ BM, BRNM ~ BM, GL ~ BRNM, BMR ~ BRNM),
  four  = DAG(dNdS ~ BM, BRNM ~ BM, GL ~ BRNM, BMR ~ BM + BRNM),
  five  = DAG(dNdS ~ BM, BRNM ~ BM, GL ~ BRNM, BMR ~ BM + BRNM + GL),
  six   = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GL ~ BRNM, BMR ~ BM),
  seven = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GL ~ BRNM, BMR ~ dNdS + BM),
  eight = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GL ~ BRNM),
  nine  = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GL ~ BRNM, BMR ~ dNdS),
  ten   = DAG(dNdS ~ BM + BRNM, BRNM ~ BM, GL ~ BRNM, BMR ~ dNdS + GL),
  eleven = DAG(BM ~ BRNM, dNdS ~ BRNM, GL ~ BM, BMR ~ BM),
  twelve = DAG(dNdS ~ BRNM, BM ~ BRNM, GL ~ BM, BMR ~ GL),
  thirteen = DAG(dNdS ~ BRNM, BM ~ BRNM, GL ~ BM, BMR ~ dNdS + GL),
  fourteen  = DAG(dNdS ~ BRNM, BM ~ BRNM, GL ~ BM, BMR ~ BRNM + BM),
  fifteen  = DAG(dNdS ~ BRNM, BM ~ BRNM, GL ~ BM, BMR ~ BRNM + BM + GL),
  sixteen   = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GL ~ BM, BMR ~ BRNM),
  seventeen = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GL ~ BM, BMR ~ dNdS + BRNM),
  eighteen = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GL ~ BM),
  nineteen  = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GL ~ BM, BMR ~ dNdS)
)
result <- phylo_path(models, data = dat, tree = tree, 
	order = c('BM', 'BRNM', 'GC', 'dNdS', 'BMR'))

summary(result)

(best_model <- best(result))
plot(best_model)

average_model <- average(result)
plot(average_model)	

plot(est_DAG(models$one, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$two, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$three, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$four, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$five, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$six, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$seven, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$eight, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$nine, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$ten, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$eleven, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twelve, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$thirteen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$fourteen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$fifteen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$sixteen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$seventeen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$eighteen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$nineteen, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.one, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.two, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.three, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.four, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.five, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.six, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.seven, data=dat, cor_fun=ape::corPagel, tree))
plot(est_DAG(models$twenty.eight, data=dat, cor_fun=ape::corPagel, tree))


(best_model <- best(result))
plot(best_model)

average_model <- average(result)
plot(average_model)	

coef_plot(best_model)

coef_plot(average_model, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_grey()
  
coef_plot(best_model, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_grey()
  
a <- (est_DAG(models$ten, data=dat, cor_fun=ape::corPagel, tree))  
  
coef_plot(a, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_grey()
  
b <- (est_DAG(models$four, data=dat, cor_fun=ape::corPagel, tree))  
  
coef_plot(b, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_grey()
  
c <- (est_DAG(models$eleven, data=dat, cor_fun=ape::corPagel, tree))  
  
coef_plot(c, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_grey()
  
 


  
	




#DAG Models
models <- list(
  one   = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ GC),
  two   = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ dNdS + GC),
  three = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BRNM),
  four  = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BM + BRNM),
  five  = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BM + BRNM + GC),
  six   = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM, BMR ~ BM),
  seven = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM, BMR ~ dNdS + BM),
  eight = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM),
  nine  = DAG(dNdS ~ BM, BRNM ~ BM + BMR, GC ~ BRNM, BMR ~ dNdS),
  ten   = DAG(dNdS ~ BM + BRNM, BRNM ~ BM, GC ~ BRNM, BMR ~ dNdS + GC),
  eleven = DAG(BM ~ BRNM, dNdS ~ BRNM, GC ~ BM, BMR ~ BM),
  twelve = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ GC),
  thirteen = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ dNdS + GC),
  fourteen  = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ BRNM + BM),
  fifteen  = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ BRNM + BM + GC),
  sixteen   = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM, BMR ~ BRNM),
  seventeen = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM, BMR ~ dNdS + BRNM),
  eighteen = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM),
  nineteen  = DAG(dNdS ~ BRNM, BM ~ BRNM + BMR, GC ~ BM, BMR ~ dNdS),
  twenty  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ BRNM),
  twenty.one  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ BM + BRNM),
  twenty.two  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ GC),
  twenty.three  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ dNdS + GC),
  twenty.four  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ dNdS + GC + BRNM),
  twenty.five  = DAG(BM ~ dNdS, GC ~ dNdS + BMR, BRNM ~ GC, BMR ~ dNdS),
  twenty.six  = DAG(BM ~ dNdS, GC ~ dNdS + BMR, BRNM ~ GC, BMR ~ BM + dNdS),
  twenty.seven  = DAG(BM ~ dNdS, GC ~ dNdS + BMR, BRNM ~ GC),
  twenty.eight  = DAG(BM ~ dNdS, GC ~ dNdS + BMR, BRNM ~ GC, BMR ~ BM)
)


#DAG Models
models <- list(
  one   = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ GC),
  two   = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ dNdS + GC),
  three = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BRNM),
  four  = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BM + BRNM),
  five  = DAG(dNdS ~ BM, BRNM ~ BM, GC ~ BRNM, BMR ~ BM + BRNM + GC),
  six   = DAG(dNdS ~ BM + BRNM, BRNM ~ BM, GC ~ BRNM, BMR ~ dNdS + GC),
  seven = DAG(BM ~ BRNM, dNdS ~ BRNM, GC ~ BM, BMR ~ BM),
  eight = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ GC),
  nine = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ dNdS + GC),
  ten  = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ BRNM + BM),
  eleven  = DAG(dNdS ~ BRNM, BM ~ BRNM, GC ~ BM, BMR ~ BRNM + BM + GC),
  twelve  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ BRNM),
  thirteen  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ BM + BRNM),
  fourteen  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ GC),
  fifteen  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ dNdS + GC),
  sixteen  = DAG(BM ~ dNdS, GC ~ dNdS, BRNM ~ GC, BMR ~ dNdS + GC + BRNM)
  )




tree = chronos(tree, lambda = 0, model = "correlated")

dat <- read.csv("shrew.compare.csv", header = T)
plot(log10(dat$BMR) ~ log10(dat$BM), xlab = "log10(Body Mass)", ylab = "log10(BMR)")
plot(log10(dat$BMR) ~ log10(dat$BRNM), xlab = "log10(Brain Mass)", ylab = "log10(BMR)")
plot(log10(dat$BMR) ~ log10(dat$dNdS), xlab = "log10(dN/dS)", ylab = "log10(BMR)")
plot(log10(dat$BMR) ~ log10(dat$GC), xlab = "log10(GC content)", ylab = "log10(BMR)")

gc <- dat[, "GC"]
bmr <- dat[, "BMR"]

# Give them names
names(bmr) <- names(gc) <- rownames(dat)

# Calculate PICs
gcPic <- pic(gc, tree)
names(gcPIC) <- rownames(dat)
bmrPic <- pic(bmr, tree)

# Make a model
picModel <- lm(bmrPic ~ gcPic - 1)
plot(bmrPic ~ gcPic)
abline(picModel)

pglsModel <- gls(log10(BMR) ~ log10(GC), correlation = corBrownian(phy=tree),
    data = dat, method = "ML")
summary(pglsModel)
plot(log10(dat$BMR) ~ log10(dat$GC), xlab = "log10(GC content)", ylab = "log10(BMR)")
abline(pglsModel)
text(x = -0.429, y = 3.7, labels = "Slope = 11.591\nIntercept = 7.195\nP-value = 0.034")

pglsModel.1 <- gls(log10(BMR) ~ log10(dNdS), correlation = corBrownian(phy=tree),
    data = dat, method = "ML")
summary(pglsModel.1)
plot(log10(dat$BMR) ~ log10(dat$dNdS), xlab = "log10(dN/dS)", ylab = "log10(BMR)")
abline(pglsModel.1)
text(x = -0.429, y = 3.7, labels = "Slope = 16.963\nIntercept = 9.498\nP-value = 0.013")

pglsModel.2 <- gls(log10(BMR) ~ log10(BRNM), correlation = corBrownian(phy=tree),
    data = dat, method = "ML")
summary(pglsModel.2)
plot(log10(dat$BMR) ~ log10(dat$BRNM), xlab = "log10(Brain Mass)", ylab = "log10(BMR)")
abline(pglsModel.2)
text(x = 0.5, y = 3.7, labels = "Slope = 0.941\nIntercept = 1.482\nP-value < 0.001")

pglsModel.3 <- gls(log10(BMR) ~ log10(BM), correlation = corBrownian(phy=tree),
    data = dat, method = "ML")
summary(pglsModel.3)
plot(log10(dat$BMR) ~ log10(dat$BM), xlab = "log10(Body Mass)", ylab = "log10(BMR)")
abline(pglsModel.3)
text(x = 2.0, y = 3.7, labels = "Slope = 0.739\nIntercept = 0.414\nP-value < 0.001")
ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
}



qplot(log10(BMR),log10(GC), data=dat,xlab = "log10(GC content)", ylab = "log10(BMR)",  size = I(7))+theme_grey()+theme(text = element_text(size=30))+ geom_abline(intercept=0.4143634, slope=0.7388666)


bmr <- log10(dat[, "BMR"])
gc <- log10(dat[, "GC"])
dNdS <- log(dat[, "dNdS"])
bm <- log10(dat[, "BM"])
brnm <- log10(dat[, "BRNM"])

names(bmr) <- names(gc) <- names(dNdS) <- names(bm) <- names(brnm) <- rownames(dat)

gcPic <- pic(gc, tree)
dndsPic <- pic(dNdS, tree)
bmrPic <- pic(bmr, tree)
bmPic <- pic(bm, tree)
brnmPic <- pic(brnm, tree)

plot(bmrPic ~ gcPic)
plot(bmrPic ~ bmPic)
plot(bmrPic ~ brnmPic)
plot(bmrPic ~ dndsPic)

a <- sma(BMR ~ GC, data = dat, log = "xy")
b <- sma(BMR ~ dNdS, data = dat, log = "xy")
c <- sma(BMR ~ BRNM, data = dat, log = "xy")
d <- sma(BMR ~ BM, data= dat, log = "xy")


e <- ma(BMR ~ GC, data = dat, log = "xy")
f <- ma(BMR ~ dNdS, data = dat, log = "xy")
g <- ma(BMR ~ BRNM, data = dat, log = "xy")
h <- ma(BMR ~ BM, data= dat, log = "xy")


par(mfrow=c(2,2))
plot(a)
text(x = 0.375, y = 4000, labels = "Slope = 27.955\nIntercept = 13.789\nR-Squared = 0.438\nP-value < 0.001")
plot(b)
text(x = 0.375, y = 4000, labels = "Slope = 34.545\nIntercept = 16.660\nR-Squared = 0.506\nP-value < 0.001")
plot(c)
text(x = 3.5, y = 4000, labels = "Slope = 1.020\nIntercept = 1.389\nR-Squared = 0.949\nP-value < 0.001")
plot(d)
text(x = 150, y = 4000, labels = "Slope = 0.822\nIntercept = 0.152\nR-Squared = 0.940\nP-value < 0.001")







