### FUNCTIONAL REGRESSION ###

# Load necessary packages ----
library(fda)
library(fastDummies)
library(ggplot2) #needs version 3.4.4 of ggplot, doesn't work on newer version
library(gridExtra)
library(cowplot)

# Load the data ----
load("SW_fd.RData")
expl.dat <- readRDS("Data/SW19-23explanatoryvar.RData")

# Pre-processing ----
expl.dat$Groupsize <- expl.dat$Groupsize/max(expl.dat$Groupsize)
expl.dat$Depth <- -expl.dat$Depth/max(abs(expl.dat$Depth))

# Split ID in groups >10, 4-9, 1-3
#sw_id10 <- names(which(table(expl.dat$ID_SW) > 9))
sw_id3 <- names(which(table(expl.dat$ID_SW) < 4))
sw_id4_9 <- names(which(table(expl.dat$ID_SW) >= 4 & table(expl.dat$ID_SW) <= 9))
#expl.dat$ID_SW[expl.dat$ID_SW %in% sw_id10] <- "many"
expl.dat$ID_SW[expl.dat$ID_SW %in% sw_id3] <- "few"
expl.dat$ID_SW[expl.dat$ID_SW %in% sw_id4_9] <- "several"

# Create dummy variables
expl.dat <- dummy_cols(expl.dat, select_columns = c("ID_SW"), remove_selected_columns = TRUE)
str(expl.dat)
summary(expl.dat)
expl.dat$ID_SW_few <- NULL

# Create xfdlist, independent variables
expl.dat <- cbind(rep(1, nrow(expl.dat)), expl.dat) 
colnames(expl.dat) <- c("intercept", colnames(expl.dat[,-1]))
p <- ncol(expl.dat)
explList <- vector("list", p)

for (i in 1:p) {
  explList[[i]] <- as.numeric(expl.dat[,i])
}

# Create BetaList, functional parameter objects
betafdPar <- fdPar(basis)
betaList <- vector("list", p)
for (i in 1:p){
  betaList[[i]] <- betafdPar
} 


# Regression ----
fregr <- fRegress(fdobj.fd, explList, betaList)

# estimated regression coefficient functions
betaList.est <- fregr$betaestlist

# Compute standard errors of coefficient functions
SigmaE <- var(t(eval.fd(time, fregr$yfdobj) - eval.fd(time, fregr$yhatfdobj)))
cifregr <- fRegress.stderr(fregr, fdobj$y2cMap, SigmaE)

betaList_upic <- betaList_lowic <- betaList.est <- fregr$betaestlist

for(i in 1:ncol(expl.dat)){
  betaList_lowic[[i]]$fd$coefs <- betaList.est[[i]]$fd$coefs - 1.96*cifregr$betastderrlist[[i]]$coefs
  betaList_upic[[i]]$fd$coefs <- betaList.est[[i]]$fd$coefs + 1.96*cifregr$betastderrlist[[i]]$coefs
}

# Plots with confidence bands 
plt <- list()

for (i in 1:p) {
  time_vector <- seq(betaList.est[[i]]$fd$basis$rangeval[1], betaList.est[[i]]$fd$basis$rangeval[2], 1)
  df_est <- data.frame(time = time_vector, value = predict(betaList.est[[i]], newdata = time_vector))
  df_upic <- data.frame(time = time_vector, value = predict(betaList_upic[[i]], newdata = time_vector))
  df_lowic <- data.frame(time = time_vector, value = predict(betaList_lowic[[i]], newdata = time_vector))
  
  plt[[i]] <- ggplot() +
    geom_line(data = df_est, aes(x = time, y = value), color = "black", linewidth = 1) +
    geom_ribbon(data = data.frame(time = df_upic$time, 
                                  ymin = df_lowic$value, 
                                  ymax = df_upic$value), 
                aes(x = time, ymin = ymin, ymax = ymax), 
                fill = "#3399CC", alpha = 0.3) +
    geom_line(data = df_upic, aes(x = time, y = value), color = "#3399CC", linewidth = 1, linetype = "dashed") +
    geom_line(data = df_lowic, aes(x = time, y = value), color = "#3399CC", linewidth = 1, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlim(range(df_est$time)) +
    ylim(-12, 13) +
    labs(x = "Time (ms)", y = "Coefficient", title = colnames(expl.dat)[i]) +
    theme_minimal()
}


grform <- c()
for (i in 1:p) {
  gr_i <- paste0("plt[[", i, "]]")
  grform <- c(grform,gr_i)
}

grform <- paste(grform, collapse= ",")
plot_list <- eval(parse(text = paste0("list(", grform, ")")))
do.call(grid.arrange, plot_list)


# Model evaluation ----

# Evaluating the overall statistical significance of the regression by means of a Permutation F-Test
fres <- Fperm.fd(fdobj.fd, explList, betaList, nperm=10)

Fperm.fd.legend <- function (yfdPar, xfdlist, betalist, wt = NULL, nperm = 5, 
                             argvals = NULL, q = 0.05, plotres = TRUE, ...) 
{
  Fnull = rep(0, nperm)
  Fnullvals = c()
  q <- 1 - q
  begin <- proc.time()
  fRegressList <- fRegress(yfdPar, xfdlist, betalist)
  elapsed.time <- max(proc.time() - begin, na.rm = TRUE)
  if (elapsed.time > 30/nperm) {
    print(paste("Estimated Computing time =", round(nperm * 
                                                      elapsed.time), "seconds."))
  }
  yhat <- fRegressList$yhatfdobj.
  if (is.null(yhat)) 
    yhat = fRegressList$yhat
  if (is.list(yhat) && ("fd" %in% names(yhat))) 
    yhat <- yhat$fd
  tFstat <- Fstat.fd(yfdPar, yhat, argvals)
  Fvals <- tFstat$F
  Fobs <- max(Fvals)
  argvals <- tFstat$argvals
  if (is.vector(yfdPar)) {
    n <- length(yfdPar)
  }
  else {
    n <- ncol(yfdPar$coefs)
  }
  for (i in 1:nperm) {
    tyfdPar <- yfdPar[sample(n)]
    fRegressList <- fRegress(tyfdPar, xfdlist, betalist)
    if (is.fd(yhat)) {
      yhat <- fRegressList$yhatfdobj
      if (is.list(yhat) && ("fd" %in% names(yhat))) 
        yhat <- yhat$fd
    }
    else {
      yhat <- fRegressList$yhat
    }
    tFstat <- Fstat.fd(tyfdPar, yhat, argvals)
    Fnullvals <- cbind(Fnullvals, tFstat$F)
    Fnull[i] <- max(Fnullvals[, i])
  }
  pval <- mean(Fobs < Fnull)
  qval <- quantile(Fnull, q)
  pvals.pts <- apply(Fvals < Fnullvals, 1, mean)
  qvals.pts <- apply(Fnullvals, 1, quantile, q)
  if (plotres) {
    if (is.fd(yfdPar)) {
      ylims <- c(min(c(Fvals, qval, qvals.pts)), max(c(Fobs, 
                                                       qval)))
      if (is.null(names(yhat$fdnames))) {
        xlab <- "argvals"
      }
      else {
        xlab <- names(yhat$fdnames)[1]
      }
      plot(argvals, Fvals, type = "l", ylim = ylims, col = 2, 
           lwd = 2, xlab = xlab, ylab = "F-statistic", 
           main = "Permutation F-Test", ...)
      lines(argvals, qvals.pts, lty = 3, col = 4, lwd = 2)
      abline(h = qval, lty = 2, col = 4, lwd = 2)
    }
    else {
      xlims <- c(min(c(Fnull, Fobs)), max(c(Fnull, Fobs)))
      hstat <- hist(Fnull, xlim = xlims, lwd = 2, xlab = "F-value", 
                    main = "Permutation F-Test")
      abline(v = Fobs, col = 2, lwd = 2)
      abline(v = qval, col = 4, lty = 2, lwd = 2)
    }
  }
  return(list(pval = pval, qval = qval, Fobs = Fobs, Fnull = Fnull, 
              Fvals = Fvals, Fnullvals = Fnullvals, pvals.pts = pvals.pts, 
              qvals.pts = qvals.pts, fRegressList = fRegressList, 
              argvals = argvals))
}

Fperm.fd.legend(fdobj.fd, explList, betaList)


# Additional model diagnostics

model_diagnostics <- function(fit){
  data_diagn = data.frame(y = as.vector(fit$yfdobj$coefs),
                          y_hat = as.vector(fit$yhatfdobj$coefs),
                          res = as.vector(fit$yfdobj$coefs)-as.vector(fit$yhatfdobj$coefs))
  
  g1 = ggplot(data=data_diagn, aes(x=y_hat, y=y))+
    theme_minimal()+
    geom_point(size=.6, alpha=0.5)+
    xlab("Predicted values")+
    ylab("Real values")+
    geom_abline(intercept = 0,
                slope = 1,
                color="red",
                linetype="dashed",
                size=1.5)
  
  g2 <- ggplot(data=data_diagn, aes(x=y_hat, y=res))+
    geom_point(size=.6, alpha=0.5)+
    geom_abline(intercept = 0,
                slope = 0,
                color="red")+
    xlab("Predicted values")+
    ylab("Residuals")+
    theme_minimal()+
    geom_hline(aes(yintercept = -1.96*sqrt(var(res))),
               col="blue",
               linetype="dashed")+
    geom_hline(aes(yintercept = 1.96*sqrt(var(res))),
               col="blue",
               linetype="dashed")
  
  g3 <- ggplot(data_diagn, aes(x=res))+
    geom_histogram(binwidth = 0.1,aes(y=..density..),
                   colour="black", fill = "lightblue")+
    stat_function(fun = dnorm, n = 101,
                  args = list(mean = 0,
                              sd = sqrt(var(data_diagn$res))),
                  size=1)+
    geom_vline(aes(xintercept=mean(res)),
               color="blue", linetype="dashed", size=.3)+
    xlab("Residuals")+ theme_minimal()+ ylab("Density")
  
  g4 <- ggplot(data_diagn, aes(sample=res))+
    stat_qq(size = 0.3) +
    theme_minimal()+
    stat_qq_line(col="red",
                 linetype = "dashed",
                 size = 1)
  return(grafici = grid.arrange(grobs = list(g1, g2, g3, g4),
                                nrow = 2))
}


model_diagnostics(fregr)




