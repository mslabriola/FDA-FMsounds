### FUNCTIONAL CLUSTERING ###

# Load necessary packages ----
library(dplyr)
library(funFEM)
library(fda)
library(gridExtra)
library(ggplot2)
library(grid)
library(tidyr)

# Load the data ----
load("sw_data_00_matrix.RData")

# Pre-processing ----

# Select sounds to use in clustering
SW.clust <- SW[,-1]
rvn.clust <- rvn.dat

# Split data into 3 groups based on duration (DeltaTime)
# and remove duration outliers
boxplot(rvn.clust$`Delta Time (s)`)
summary(rvn.clust$`Delta Time (s)`) 

bstats <- boxplot.stats(rvn.clust$`Delta Time (s)`)
q1 <- bstats$stats[2]
q3 <- bstats$stats[4]
upwh <- bstats$stats[5]


long_ind <- c(rvn.clust$`Delta Time (s)`>q3&
                rvn.clust$`Delta Time (s)`<= upwh)
short_ind <- rvn.clust$`Delta Time (s)`<= q1
out_ind <- rvn.clust$`Delta Time (s)`> upwh
SW_long <- cbind(time , SW.clust[,long_ind ])
SW_short <- cbind(time, SW.clust[,short_ind])
SW_medium <- cbind(time , SW.clust[,!c( long_ind | short_ind |out_ind )])
rvn.clust_long <- rvn.clust[long_ind,]
rvn.clust_short <- rvn.clust[short_ind,]
rvn.clust_medium <- rvn.clust[!c(long_ind|short_ind|out_ind),]


# Adapt max length of each group to remove additional 0s
# the vectors should be as long as the longest sound in the group
max_delta_ind_short <- which (rvn.clust_short$`Delta Time (s)`
                                   == max(rvn.clust_short$`Delta Time (s)`))[1]
max_delta_ind_medium <- which (rvn.clust_medium$`Delta Time (s)`
                                    == max(rvn.clust_medium$`Delta Time (s)`))[1]
max_delta_ind_long <- which (rvn.clust_long$`Delta Time (s)`
                                 == max(rvn.clust_long$`Delta Time (s)`))[1]

max_SW_short <- SW_short[, max_delta_ind_short +1]
max_SW_medium <- SW_medium[, max_delta_ind_medium +1]
max_SW_long <- SW_long[, max_delta_ind_long +1]
SW_short <- SW_short[1:( which(max_SW_short == 0)[1] -1) , ]
SW_medium <- SW_medium[1:( which (max_SW_medium == 0)[1] -1) , ]
SW_long <- SW_long [1:( which (max_SW_long == 0)[1] -1) , ]
time_short <- SW_short[ ,1]
time_medium <- SW_medium[ ,1]
time_long <- SW_long[ ,1]

sw_id_long <- names(table(rvn.clust_long$ID_Category))
sw_id_medium <- names(table(rvn.clust_medium$ID_Category))
sw_id_short <- names(table(rvn.clust_short$ID_Category))
sw_id3_long <- sw_id_long[table(rvn.clust_long$ID_Category)>2]
sw_id3_medium <- sw_id_medium[table(rvn.clust_medium$ID_Category) >2]
sw_id3_short <- sw_id_short[table(rvn.clust_short$ID_Category) >2]
SW_long_clust <- SW_long[,-1]
SW_medium_clust <- SW_medium[,-1]
SW_short_clust <- SW_short[,-1]
SW_long <- cbind(time_long,
                    SW_long_clust[,rvn.clust_long$ID_Category %in% sw_id3_long])
SW_medium <- cbind(time_medium,
                 SW_medium_clust[,rvn.clust_medium$ID_Category %in% sw_id3_medium])
SW_short <- cbind(time_short,
                 SW_short_clust[,rvn.clust_short$ID_Category %in% sw_id3_short])

rvn.dat_long <- rvn.clust_long[rvn.clust_long$ID_Category %in% sw_id3_long ,]
rvn.dat_medium <- rvn.clust_medium[rvn.clust_medium$ID_Category %in% sw_id3_medium ,]
rvn.dat_short <- rvn.clust_short[rvn.clust_short$ID_Category %in% sw_id3_short ,]


# number of unique contours in each group
length(unique(rvn.dat_long$ID_Category)) 
length(unique(rvn.dat_medium$ID_Category)) 
length(unique(rvn.dat_short$ID_Category)) 

# total number of unique contours
length(unique(c(unique(rvn.dat_long$ID_Category),unique(rvn.dat_medium$ID_Category),
         unique(rvn.dat_short$ID_Category))))

rm(long_ind,max_delta_ind_long, max_delta_ind_medium, max_delta_ind_short,
   max_SW_long,max_SW_medium, max_SW_short, out_ind, short_ind, sw_id3_long,
   sw_id3_short, sw_id3_medium, sw_id_long, sw_id_medium, sw_id_short, SW_long_clust,
   SW_medium_clust, SW_short_clust)

save.image(file = "sw_data_deltatime.RData")

# Tuning procedure ----

#authors suggestion: run this section on a high-performance supercomputer
#(we used TeraStat, owned by the Department of Statistical Sciences of 
#Sapienza University of Rome)

#the section reports only the procedure for short whistles and the model AB, 
#but the exact same code can be used for medium and long whistles and all the  
#other funFEM models (ABk, AjB, AjBk, AkB, AkBk, AkjB, AkjBk, DB, DBk, DkB, DkBk)
#by replacing the appropriate terms

load("sw_data_deltatime.RData")

# define number of Fourier basis functions and number of clusters to consider
seq_basis <- seq(5,30, by=1)
num_clusters <- seq(3,6, by=1)

# Create funFEM models with all combinations
models_BIC_short <- c()
models_ICL_short <- c()

for(num_basis in seq_basis){
  
  basis = create.fourier.basis(range(SW_short[,1]) ,
                               nbasis = num_basis)
  fdobj = smooth.basis(SW_short[,1], SW_short[,-1], basis)
  
  set.seed(123)
  model <- funFEM(fdobj$fd,
                  K = num_clusters,
                  model = "AB",
                  init = "random",
                  maxit = 300)
  
  models_BIC_short <- cbind(models_BIC_short, model$allCriterions$bic)
  models_ICL_short <- cbind(models_ICL_short, model$allCriterions$icl)
}

rownames(models_BIC_short) <- paste0('C=', num_clusters)
colnames(models_BIC_short) <- paste0(seq_basis,rep(" basis",length(seq_basis)))
rownames(models_ICL_short) <- paste0('C=', num_clusters)
colnames(models_ICL_short) <- paste0(seq_basis,rep(" basis",length(seq_basis)))

short_AB <- list("models_BIC_short"=models_BIC_short,"models_ICL_short"=models_ICL_short)
save(short_AB, file = "short_AB.RData")


# Model selection ----
# to perform on the total output of the previous section performed on all models

#the section reports only the procedure for short whistles, but the exact same 
#code can be used for medium and long whistles by replacing the appropriate terms

lapply(paste0("funFEM_models/",list.files("funFEM_models/")), load, environment())

models <- c("AB", "ABk", "AjB", "AjBk", "AkB", "AkBk", "AkjB", "AkjBk", "DB", 
            "DBk", "DkB", "DkBk")

criteria <- c("BIC", "ICL")

for (m in 1:length(models)) {
  for (c in 1:length(criteria)) {
    assign(paste0("models_",criteria[c],"_",models[m]), as.data.frame(get(paste0("short_",models[m]))[[paste0("models_",criteria[c],"_short")]]))
  }
}


# Select best model based on BIC and ICL values
i <- 1

vector_BIC <- c()

for (m in 1:length(models)) {
  BIC_m <- sort(as.matrix(get(paste0("models_BIC_",models[m])), na.rm=TRUE), decreasing = TRUE)[i]
  vector_BIC <- c(vector_BIC, BIC_m)
}
models_BIC <- data.frame(vector_BIC, row.names = models)


vector_ICL <- c()

for (m in 1:length(models)) {
  ICL_m <- sort(as.matrix(get(paste0("models_ICL_",models[m])), na.rm=TRUE), decreasing = TRUE)[i]
  vector_ICL <- c(vector_ICL, ICL_m)
}
models_ICL <- data.frame(vector_ICL, row.names = models)


for (m in 1:length(models)) {
  for (c in 1:length(criteria)) {
    assign(paste0("ind_",criteria[c],"_",models[m]), 
           which(get(paste0("models_",criteria[c],"_",models[m])) == sort(as.matrix(get(paste0("models_",criteria[c],"_",models[m])), na.rm=TRUE), 
                                       decreasing = TRUE)[i], arr.ind=TRUE))
  }
}

clusters_bic <- c()
basis_bic <- c()
for (m in 1:length(models)) {
  clust_m <- rownames(models_BIC_AB)[get(paste0("ind_BIC_", models[m]))[1]]
  clusters_bic <- c(clusters_bic, clust_m)
  
  basis_m <- colnames(models_BIC_AB)[get(paste0("ind_BIC_", models[m]))[2]]
  basis_bic <- c(basis_bic, basis_m)
}


clusters_icl <- c()
basis_icl <- c()
for (m in 1:length(models)) {
  clust_m <- rownames(models_ICL_AB)[get(paste0("ind_ICL_", models[m]))[1]]
  clusters_icl <- c(clusters_icl, clust_m)
  
  basis_m <- colnames(models_ICL_AB)[get(paste0("ind_ICL_", models[m]))[2]]
  basis_icl <- c(basis_icl, basis_m)
}


models_BIC <- data.frame(models_BIC, clusters_bic, basis_bic)
models_ICL <- data.frame(models_ICL, clusters_icl, basis_icl)

models_BIC
models_ICL

best_model_BIC <- models_BIC[which(models_BIC$vector_BIC==max(models_BIC$vector_BIC)),]
best_model_ICL <- models_ICL[which(models_ICL$vector_ICL==max(models_ICL$vector_ICL)),]

best_model_BIC
best_model_ICL

# Extract parameters to use in the clustering (BIC)
mod <- row.names(best_model_BIC) 

num_clusters <- strsplit(best_model_BIC[[2]], "=")
num_clusters <- as.numeric(num_clusters[[1]][2])

num_basis <- strsplit(best_model_BIC[[3]], " ")
num_basis <- as.numeric(num_basis[[1]][1])


# Clustering ----

load("sw_data_deltatime.RData")

# choose the group (short, medium or long)
rvn.dat <- rvn.dat_short
SW <- SW_short

basis = create.fourier.basis(range(SW[,1]), nbasis = num_basis)
fdobj = smooth.basis(SW[,1], SW[,-1], basis)

set.seed(123)
model <- funFEM(fdobj$fd,
                K = num_clusters,
                model = "DkBk",
                init = "random",
                maxit = 300)


# Plot sounds in the first two dimensions of the discriminative functional space
# colours correspond to the assigned cluster

plot(t(fdobj$fd$coefs) %*% model$U, col=model$cls, pch=19, 
     xlab='First dimension', ylab='Second dimension', lwd=4)

legend("topright",legend=as.character(1:num_clusters), 
       col=1:max(model$cls), lwd=3, pch=19)


# Plot median contours of each cluster

nclust <- max(model$cls)
ind_clust <- model$cls

time <- SW[,1]
rvn.dat_clust_list <- list()


for (i in 1:nclust) {
  rvn.dat_clust_list[[i]] <- rvn.dat[ind_clust == i,]
}

sw_eval = eval.fd(SW[, 1], fdobj$fd)
sw_eval_list <- list()
df_long <- list()
sw_median <- list()

for (i in 1:nclust) {
  sw_eval_i <- as.data.frame(sw_eval[, ind_clust == i])
  sw_eval_i$time <- time
  sw_eval_list[[i]] <- sw_eval_i
  df_long[[i]] <- gather(sw_eval_list[[i]], key = "sw_variable", value = "sw_eval", -time)
  sw_median[[i]] <- data.frame(y = apply(sw_eval_list[[i]], 1, median), x = time)
}


for (i in 1:nclust) {
  assign(paste0("plot_clust_", i), ggplot(data = df_long[[i]], aes(x = time, y = sw_eval)) +
           geom_line(lwd = 1.2, col = "lightgrey", aes(group = sw_variable)) +
           labs(x = "Time (ms)", y = "Frequency (kHz)") +
           geom_line(data = sw_median[[i]], show.legend = FALSE, color = i, lwd = 3,aes(x=x, y=y)) +
           ggtitle(paste0("Cluster ", i)) +
           theme_minimal() +
           theme(axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15),
                 axis.title = element_text(size = 15),
                 plot.title = element_text(size = 20)))
}

grform <- c()
for (i in 1:nclust) {
  gr_i <- paste0("plot_clust_", i)
  grform <- c(grform,gr_i)
}

grform <- paste(grform, collapse= ",")
plot_list <- eval(parse(text = paste0("list(", grform, ")")))
do.call(grid.arrange, plot_list)


