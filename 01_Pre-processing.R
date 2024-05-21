### PRE-PROCESSING ###

# Load necessary packages ----
library(Rraven)
library(dplyr)
library(ggplot2) 

# Load the data ----

# import all selection tables in the wd 
# (sounds previously segmented and quality graded in Raven Pro)
rvn.dat <- imp_raven(path = "Data/", all.data = TRUE, name.from.file = TRUE, ext.case = "lower")

# OR import data frame containing all selection tables
rvn.dat <- readRDS("Data/rvn.dat.RData")

# select only signature whistles (SW) with quality >= 2 
rvn.dat <- filter(rvn.dat, Category == "SW", Quality >= 2)

# FM pattern extraction ----

# extract the "Peak Freq Contour (Hz)" column
# every row contains sound file, selection number, and PFC values
peakfreq <- extract_ts(rvn.dat, ts.column = "Peak Freq Contour (Hz)", equal.length = FALSE)

# create time vector in ms
time <- seq(0, 2.6*(ncol(peakfreq)-3), by=2.6)

# adapt the dataframe & convert Hz to kHz
peakfreq <- t(peakfreq[,-(1:2)]/1000)
SW <- cbind(time, peakfreq)

# take every 2nd time point
SW <- SW[seq(1, nrow(SW), 2), ]
time <- SW[,1]


# Treatment of outliers ----

# Create median filters
# with t = 2.8 
medianfilter <- function(w, thresh = 2.8, n = 5){ 
  lengthSW <- (which(is.na(w))[1])-1
  if(is.na(lengthSW)){lengthSW<-length(w)}
  for (i in (n+1):(lengthSW-1)){
    median1 <- median(w[(i-n):(i-1)], na.rm = TRUE)
    difftemp <- abs(w[i]-median1)
    ifelse(difftemp>thresh, w[i]<-mean(w[i-1], w[i+1]), next)
  }
  return(w)
}

# with t = sd of the input vector (for shorter sounds)
medianfilter_sd <- function(w,n = 5){ 
  thresh <- sd(w, na.rm = TRUE)
  lengthSW <- (which(is.na(w))[1])-1
  if(is.na(lengthSW)){lengthSW<-length(w)}
  for (i in (n+1):(lengthSW-1)){
    median1 <- median(w[(i-n):(i-1)], na.rm = TRUE)
    difftemp <- abs(w[i]-median1)
    ifelse(difftemp>thresh, w[i]<-mean(w[i-1], w[i+1]), next)
  }
  return(w)
}


# plot all sounds unprocessed, with median filter t=2.8 and t=sd

pdf("Plots_SW.pdf")
par(mfrow=c(3,2))
for (i in 2:ncol(SW)){
  plot(x=time[!is.na(SW[,i])], y=na.omit(SW[,i]), xlab = i)
}
dev.off() 

SW_mf <- apply(SW, 2, medianfilter)

pdf("Plots_SW_MEDIAN.pdf")
par(mfrow=c(3,2))
for (i in 2:ncol(SW_mf)){
  plot(x=time[!is.na(SW_mf[,i])], y=na.omit(SW_mf[,i]), xlab = i)
}
dev.off() 


SW_sd <- apply(SW, 2, medianfilter_sd)

pdf("Plots_SW_SD.pdf")
par(mfrow=c(3,2))
for (i in 2:ncol(SW_sd)){
  plot(x=time[!is.na(SW_sd[,i])], y=na.omit(SW_sd[,i]), xlab = i)
}
dev.off() 


# define sounds that don't need median filtering
without <- c(5,12,28,29,31,32,33,35,36,38,43,51,52,53,54,56,60,62,63,64,74,75,107,148,156,157,158,159,
             160,161,164,165,166,171,189,198,203,207,208,209,210,212,213,229,230,239)

# define sounds which should be processed with median filter t=sd
with_sd <- c(5,6,16,23,24,41,45,59,162,163,186,187,190,191,200,206,214,215,217,221,223,240,241)


# apply median filters
ind <- 2:ncol(SW)
not <- c(without, with_sd)
with <- setdiff(ind, not)
SW[,with] <- apply(SW[,with], 2, medianfilter)
SW[,with_sd] <- apply(SW[,with_sd], 2, medianfilter_sd)

rm(i, ind, medianfilter, medianfilter_sd, not, SW_mf, SW_sd, with, with_sd, without)


# Zero-filling ----

# replace NA values with 0's
SW[is.na(SW)] <- 0

saveRDS(SW, "SW_00.RData")
save.image("sw_data_00_matrix.RData")



