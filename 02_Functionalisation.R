### FUNCTIONALISATION OF TONAL SOUNDS ###

# Load necessary packages ----
library(ggplot2) 
library(gridExtra)
library(fda)

# Load the data ----
SW <- readRDS("SW_00.RData")
time <- SW[, "time"]

# Select the optimal number of basis functions based on RMSE ----

min_basis = 4 #minimun n basis to consider
max_basis = 80 #maximum n basis to consider

# RMSE function
fun_rmse <- function(r) {sqrt(mean((r)^2))}

# Compute RMSE for each n basis
rmse <- c()
for (i in min_basis:max_basis) {
  basis = create.bspline.basis(range(SW[,1]), nbasis = i)
  fdobj = smooth.basis(SW[,1], SW[,-1], basis)
  sw_eval <- eval.fd(SW[, 1], fdobj$fd)
  res_mat <- sw_eval - SW[,-1]
  temp_rmse <- apply(res_mat, 2, fun_rmse)
  rmse <- c(rmse, mean(temp_rmse))
}

RMSE_df <- data.frame(nbasis = min_basis:max_basis, RMSE = rmse)

# RMSE plot 
# the vertical red line intersects the elbow of the curve, 
# corresponding to the optimal n basis
ggplot(RMSE_df, aes(x = nbasis, y = RMSE)) + 
  geom_point() +
  theme_minimal() +
  geom_vline(xintercept = 25, colour="red")

rm(fun_rmse, i, max_basis, min_basis, res_mat, rmse, RMSE_df, temp_rmse)

# Functionalisation ----

# define B-spline basis functions
n_basis = 25 
basis = create.bspline.basis(range(SW[,1]), nbasis = n_basis)

# create functional object through data smoothing
fdobj = smooth.basis(SW[,1], SW[,-1], basis)
fdobj.fd <- fdobj$fd
sw_eval = eval.fd(SW[, 1],fdobj.fd)

save.image("SW_fd.RData")


# plot functionalised sounds
sample <- c(170,190) #select sounds to plot
vcol <- c("darkcyan", "pink")

for (i in 1:length(sample)) {
  df_func <- data.frame(time = time, freq = SW[,sample[i]], func = sw_eval[,sample[i]-1])
  assign(paste0("p",i), ggplot(df_func, aes(x = time)) +
           geom_point(aes(y = freq, color = "functionalised"), size = 4) +
           geom_path(aes(y = func, color = "discretised"), size = 2) +
           labs(x = "Time (ms)", y = "Frequency (Hz)") +
           scale_color_manual(name = "Legend",
                              values = c("functionalised" = vcol[1],
                                         "discretised" = vcol[2])) +  
           theme_classic() +
           theme(axis.text = element_text(size = 16),    
                 axis.title = element_text(size = 20)) + 
           theme(legend.position = "none"))
  }

grform <- c()
for (i in 1:length(sample)) {
  gr_i <- paste0("p", i)
  grform <- c(grform,gr_i)
}

grform <- paste(grform, collapse= ",")
plot_list <- eval(parse(text = paste0("list(", grform, ")")))
do.call(grid.arrange, plot_list)


