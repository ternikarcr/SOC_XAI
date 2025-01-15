#%%%%%%%%%%%%%%%%%#
####Section 1####
#%%%%%%%%%%%%%%%%%#
packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", "moments", "caTools",
                        "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
                        "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
                        "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
                        "installr", "resemble", "prospectr", "magrittr", "doParallel", 
                        "parallel", "foreach", "ggplot2", "tidyr", "pls", "ChemoSpec", "MASS", "Johnson")
# Install and load packages
# install.packages(packages_to_install, dependencies = TRUE)
# Load the installed packages
sapply(packages_to_install, require, character.only = TRUE)

##FUNCTIONS
# Function to calculate the minimum indices
get_lowest_indices <- function(data, no_neighbor = 5) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  unique_min_indices <- unique(min_indices)
  return(unique_min_indices)
}

# Function to calculate R², RMSE, RPD, and RPIQ
calculate_statistics <- function(observed, predicted) {
  r_squared <- cor(observed, predicted)^2
  rmse <- sqrt(mean((observed - predicted)^2))
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  deviation <- observed - mean_observed
  rpd <- sd(observed) / rmse
  iqr_observed <- IQR(observed)
  rpiq <- iqr_observed / rmse
  bias <- mean_observed - mean_predicted
  mae <- mean(abs(observed - predicted))
  result <- list(R_squared = r_squared,RMSE = rmse,RPD = rpd,RPIQ = rpiq,Bias = bias,MAE = mae)
  return(result)
}

calculate_statistics_result <- function(observed_values1,Y2_hat1, ncomp1) {
  cal_statistics_result <- NULL
  for (i in seq(1, ncomp1, by = 1)) {
    predicted_values1 <- Y2_hat1[, 1, i] 
    cal_statistics_result[[i]] <- calculate_statistics(observed_values1, predicted_values1)
  }
  return(cal_statistics_result)
}

setwd("C:/Refined/Draft_3/")
#OSSL library
OSSL = read_csv("OSSL.csv")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,60:length(foo)]))
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(502, 2500, by = 2))
rm(OSSL,foo)
OSSL = foo1
rm(foo1)

#Data provision
IP = OSSL
#Component selection
property = "SOC"
summary(IP[,property])

IP_X1 = IP$spc[!is.na(IP[,property]),]
IP_Y1 = IP[,property][!is.na(IP[,property])]
summary(IP_Y1)
IP_X_O = IP_X1[order(IP_Y1), ]
IP_Y_O = IP_Y1[order(IP_Y1)]
#For SOC the values were made to selected >0.1
IP_Y1[IP_Y1 < 0.1] <- NA
IP_Y = IP_Y1[!is.na(IP_Y1)]
IP_X = IP_X1[!is.na(IP_Y1),]
summary(IP_Y)
IP_X_O = IP_X[order(IP_Y), ]
IP_Y_O = IP_Y[order(IP_Y)]


#One time for data analysis
length(IP_Y)
print(c(min(IP_Y), mean(IP_Y), max(IP_Y), sd(IP_Y), sd(IP_Y)/mean(IP_Y)), skewness(IP_Y), kurtosis(IP_Y))
pc = prcomp(IP_X)
principal_components = pc$x
explained_variance = pc$sdev^2
explained_variance_per = explained_variance / sum(explained_variance)
write.csv(explained_variance_per, "explained_variance_per.csv", row.names = FALSE)
principal_components_10 = principal_components[,1:10]
write.csv(principal_components_10, "principal_components.csv", row.names = FALSE)

#For plotting representative spectra's
#Summary of all spectra
foo = summary(IP_X_O)
clean_summary <- function(x) {
  as.numeric(gsub(".*?:", "", x))  # Remove text before ":" and convert to numeric
}
cleaned_numeric_matrix <- apply(foo, 2, clean_summary)
rownames(cleaned_numeric_matrix) = c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max")
write.csv(cleaned_numeric_matrix, "spectra_summary.csv", row.names = TRUE)
#Mean of spectra based on SOC values
result <- IP_Y[IP_Y >= 4.43]
print(mean(result))
foo2 = IP_X_O[order(result),]
SOC = colMeans(foo2)
write.csv(SOC, "spectra_summary_SOC.csv", row.names = FALSE)
rm(foo, cleaned_numeric_matrix, foo, foo1, foo2, result)

#Split 8:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

####1NEAREST NEIGHBOUR####
##Training##
foo = X1
foo1 = Y1

summary(foo1)
foo100 = diff(foo1[-1])
summary(foo100)
print(c((mean(foo100)), max(diff(foo100))))

#Optional or just check once when doing it for the 1sst iteration
#optimal_sel =  list(method = "opc", value = min(length(foo1)-4,100))
##optimal_sel =  list(method = "manual", value = 27)
##no_components = 27
#pls_tr_opc = ortho_projection(Xr = foo, Yr = foo1, method = "pls", pc_selection = optimal_sel)
#pls_tr_opc
#pls_comp_num = as.matrix(as.numeric(pls_tr_opc$opc_evaluation[,1]))
#pls_comp_value = as.matrix(as.numeric(pls_tr_opc$opc_evaluation[,2]))
#old_par = par("mfrow")
#par(mfrow = c(1,1))
#matplot(x = pls_comp_num, y = pls_comp_value, xlab = "No of PLS components",ylab = "RMSD of Yr",type = "h", lty = 1, col = "#FF1A00CC")
#title("method=pls.opc")
##Finding which k has the minimum rmse
#no_components = which.min(as.numeric(pls_tr_opc$opc_evaluation[,2]))
#print(paste("Number of components:", no_components))

#Comparing all the dissimilarity measures#
#PC dissimilarity with default settings (variance-based no. of components)
pcad = dissimilarity(foo, diss_method = "pca", scale = TRUE, return_projection = TRUE)
#PLS dissimilarity with default settings (variance-based no of components)
plsd = dissimilarity(foo, diss_method = "pls", Yr = foo1,scale = TRUE, return_projection = TRUE)
#PC dissimilarity with optimal selection of components
opc_sel =  list(method = "opc", value = min(length(foo1)-4,100))
o_pcad = dissimilarity(foo,diss_method = "pca",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
cat("Number of components as selected by Minimum RMSE in the NN algorithm:\n")
print(dim(o_plsd[["projection"]][["scores"]])[2])
#Correlation dissimilarity 
cd = dissimilarity(foo, diss_method = "cor", scale = TRUE, return_projection = TRUE)
#Moving window correlation dissimilarity 
mcd = dissimilarity(foo, diss_method = "cor", ws = 51, scale = TRUE, return_projection = TRUE)
#Euclidean dissimilarity 
ed = dissimilarity(foo, diss_method = "euclid", scale = TRUE, return_projection = TRUE)
#Cosine dissimilarity 
cosd = dissimilarity(foo, diss_method = "cosine", scale = TRUE, return_projection = TRUE)
#Spectral information divergence/dissimilarity 
sinfd = dissimilarity(foo, diss_method = "sid", scale = TRUE, return_projection = TRUE)

#Evaluations
Y_matrix = as.matrix(foo1)
ev = NULL
ev[["pcad"]] = sim_eval(pcad$dissimilarity, side_info = Y_matrix)
ev[["plsd"]] = sim_eval(plsd$dissimilarity, side_info = Y_matrix)
ev[["o_pcad"]] = sim_eval(o_pcad$dissimilarity, side_info = Y_matrix)
ev[["o_plsd"]] = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)
ev[["cd"]] = sim_eval(cd$dissimilarity, side_info = Y_matrix)
ev[["mcd"]] = sim_eval(mcd$dissimilarity, side_info = Y_matrix)
ev[["ed"]] = sim_eval(ed$dissimilarity, side_info = Y_matrix)
ev[["cosd"]] = sim_eval(cosd$dissimilarity, side_info = Y_matrix)
ev[["sinfd"]] = sim_eval(sinfd$dissimilarity, side_info = Y_matrix)

#Tabulating, saving the prediction files and plotting
statistics_result = NULL
for (label in names(ev)) {
  observed_values = ev[[label]]$first_nn[, 1]
  predicted_values = ev[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_calib = bind_rows(statistics_result, .id = "Methods")
print(r_calib)
write.csv(r_calib, "r_calib.csv", row.names = FALSE)

Y2_hat_pcad = ev[["pcad"]]$first_nn[, 2]
Y2_hat_plsd = ev[["plsd"]]$first_nn[, 2]
Y2_hat_o_pcad = ev[["o_pcad"]]$first_nn[, 2]
Y2_hat_o_plsd = ev[["o_plsd"]]$first_nn[, 2]
Y2_hat_cd = ev[["cd"]]$first_nn[, 2]
Y2_hat_mcd = ev[["mcd"]]$first_nn[, 2]
Y2_hat_ed = ev[["ed"]]$first_nn[, 2]
Y2_hat_cosd = ev[["cosd"]]$first_nn[, 2]
Y2_hat_sinfd = ev[["sinfd"]]$first_nn[, 2]
Y2_foo = cbind(foo1,Y2_hat_pcad,Y2_hat_plsd,Y2_hat_o_pcad,Y2_hat_o_plsd,Y2_hat_cd,Y2_hat_mcd,Y2_hat_ed,Y2_hat_cosd,Y2_hat_sinfd)
write.csv(Y2_foo, "Y1_hat_NN.csv", row.names = FALSE)
rm(Y2_foo,Y2_hat_pcad,Y2_hat_plsd,Y2_hat_o_pcad,Y2_hat_o_plsd,Y2_hat_cd,Y2_hat_mcd,Y2_hat_ed,Y2_hat_cosd,Y2_hat_sinfd)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
#par(mfrow = c(3, 3))
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("SOC Observed, %", "SOC Predicted(1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))])
             title(label)
             grid()
             abline(0, 1)
             #text((max(xy[,1])-10), (max(xy[,2])-10),labels = paste("RMSD:", round(x[[label]]$eval[1],3), "\nR:", round(x[[label]]$eval[2],3), "\nR2:", round(x[[label]]$eval[2]^2,3)),pos = 1,col = "black", cex = 1)
           },
           x = ev)
par(mfrow = c(1, 1))

##Testing##
#PC dissimilarity with default settings (variance-based no. of components)
pcad1 = dissimilarity(Xr=foo, Xu=X2, diss_method = "pca", scale = TRUE, return_projection = TRUE)
#PLS dissimilarity with default settings (variance-based no of components)
plsd1 = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls", scale = TRUE, return_projection = TRUE)
#PC dissimilarity with optimal selection of components
opc_sel =  list(method = "opc", value = min(length(Y2)-4,100))
o_pcad1 = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pca",pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd1 = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
#Correlation dissimilarity 
cd1 = dissimilarity(Xr=foo, Xu=X2, diss_method = "cor", scale = TRUE, return_projection = TRUE)
#Moving window correlation dissimilarity 
mcd1 = dissimilarity(Xr=foo, Xu=X2, diss_method = "cor", ws = 51, scale = TRUE, return_projection = TRUE)
#Euclidean dissimilarity 
ed1 = dissimilarity(Xr=foo, Xu=X2, diss_method = "euclid", scale = TRUE, return_projection = TRUE)
#Cosine dissimilarity 
cosd1 = dissimilarity(Xr=foo, Xu=X2, diss_method = "cosine", scale = TRUE, return_projection = TRUE)
#Spectral information divergence/dissimilarity 
sinfd1 = dissimilarity(Xr=foo, Xu=X2, diss_method = "sid", scale = TRUE, return_projection = TRUE)

#Evaluations
ev1 = NULL
observed_values = Y2
ev1[["pcad1"]]$first_nn = cbind(observed_values, foo1[apply(pcad1$dissimilarity, 2, which.min)])
ev1[["plsd1"]]$first_nn = cbind(observed_values, foo1[apply(plsd1$dissimilarity, 2, which.min)])
ev1[["o_pcad1"]]$first_nn = cbind(observed_values, foo1[apply(o_pcad1$dissimilarity, 2, which.min)])
ev1[["o_plsd1"]]$first_nn = cbind(observed_values, foo1[apply(o_plsd1$dissimilarity, 2, which.min)])
ev1[["cd1"]]$first_nn = cbind(observed_values, foo1[apply(cd1$dissimilarity, 2, which.min)])
ev1[["mcd1"]]$first_nn = cbind(observed_values, foo1[apply(mcd1$dissimilarity, 2, which.min)])
ev1[["ed1"]]$first_nn = cbind(observed_values, foo1[apply(ed1$dissimilarity, 2, which.min)])
ev1[["cosd1"]]$first_nn = cbind(observed_values, foo1[apply(cosd1$dissimilarity, 2, which.min)])
ev1[["sinfd1"]]$first_nn = cbind(observed_values, foo1[apply(sinfd1$dissimilarity, 2, which.min)])

statistics_result = NULL
for (label in names(ev1)) {
  observed_values = ev1[[label]]$first_nn[, 1]
  predicted_values = ev1[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_valid = bind_rows(statistics_result, .id = "Methods")
print(r_valid)
write.csv(r_valid, "r_valid.csv", row.names = FALSE)

Y2_hat_pcad = foo1[apply(pcad1$dissimilarity, 2, which.min)]
Y2_hat_plsd = foo1[apply(plsd1$dissimilarity, 2, which.min)]
Y2_hat_o_pcad = foo1[apply(o_pcad1$dissimilarity, 2, which.min)]
Y2_hat_o_plsd = foo1[apply(o_plsd1$dissimilarity, 2, which.min)]
Y2_hat_cd = foo1[apply(cd1$dissimilarity, 2, which.min)]
Y2_hat_mcd = foo1[apply(mcd1$dissimilarity, 2, which.min)]
Y2_hat_ed = foo1[apply(ed1$dissimilarity, 2, which.min)]
Y2_hat_cosd = foo1[apply(cosd1$dissimilarity, 2, which.min)]
Y2_hat_sinfd = foo1[apply(sinfd1$dissimilarity, 2, which.min)]
Y2_foo = cbind(Y2,Y2_hat_pcad,Y2_hat_plsd,Y2_hat_o_pcad,Y2_hat_o_plsd,Y2_hat_cd,Y2_hat_mcd,Y2_hat_ed,Y2_hat_cosd,Y2_hat_sinfd)
write.csv(Y2_foo, "Y2_hat_NN.csv", row.names = FALSE)
rm(Y2_foo,Y2_hat_pcad,Y2_hat_plsd,Y2_hat_o_pcad,Y2_hat_o_plsd,Y2_hat_cd,Y2_hat_mcd,Y2_hat_ed,Y2_hat_cosd,Y2_hat_sinfd)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
#par(mfrow = c(3, 3))
p = sapply(names(ev1), 
           FUN = function(x, label, labs = c("SOC Observed, %", "SOC Predicted(1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev1))])
             title(label)
             grid()
             abline(0, 1)
             #text((max(xy[,1])-10), (max(xy[,2])-10),labels = paste("RMSD:", round(x[[label]]$eval[1],3), "\nR:", round(x[[label]]$eval[2],3), "\nR2:", round(x[[label]]$eval[2]^2,3)),pos = 1,col = "black", cex = 1)
           },
           x = ev1)
par(mfrow = c(1, 1))


####2PLSR####
##Training##
foo = X1
foo1 = Y1
model = plsr(foo1 ~ foo, ncomp = min(length(foo1)-4,100), validation = "CV", segments = 10)  
summary(model)
explvar(model) #This is output of Xvar/Xtotvar
plot(RMSEP(model), legendpos = "topright")
plot(model$validation$PRESS[1,], legendpos = "topright")
par(mfrow = c(1,1))
comp_value = RMSEP(model)$val[1,1,]
matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Y", type = "h", lty = 1, col = "#FF1A00CC")

##DECIDE on the method to be used for choosing the number of components##
no_components = which.min(comp_value)-1
no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
cat("Number of components as selected by Minimum, onesigma and permutation:", no_components,no_components_onesigma,no_components_permut,"\n")
##DECISION ENDS##


####VIP Calculations####
calc_vip_custom <- function(model) {
  if (!inherits(model, "mvr")) {
    stop("Input must be an object of class 'mvr'.")
  }
  W <- model$loading.weights
  T <- model$scores
  Q <- model$Yloadings
  n_pred <- ncol(W)  # Number of predictors
  # Variance explained by each component
  TTQ_squared <- diag(t(T) %*% T %*% t(Q) %*% Q)
  total_variance <- sum(TTQ_squared)
  # Compute VIP scores
  Wnorm2 <- colSums(W^2)
  vip <- sqrt(n_pred * rowSums((W^2 / Wnorm2) * TTQ_squared) / total_variance)
  return(vip)
}

# Calculate VIP scores using the model
vip_scores <- calc_vip_custom(model)

# Print VIP scores
vip_scores
plot(vip_scores)
write.csv(vip_scores, "VIP.csv", row.names = FALSE)


#Plots
colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
plot(foo1, model$fitted.values[,1,no_components], xlab = "SOC Observed, %", ylab = "SOC Predicted (PLSR), %", col = colours1[9])
title("PLSR calibration")
grid()
abline(0, 1)
#text(90, 30, labels = paste("N_comp:", no_components, "\nRMSD:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$RMSE,3), "\nR2:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$R_squared,3)),pos = 1,col = "black", cex = 1.2)

plot(model, plottype = "scores", comps = 1:3)
plot(model, plottype = "loadings", comps = 1:3, legendpos = "topright",
     labels = "numbers", xlab = "Wavelength (nm)", ylab= "Loading value")
abline(h = 0)
plot(model, plottype = "coef", comps = 1:3, legendpos = "topright",
     labels = "numbers", xlab = "Wavelength (nm)", ylab= "Regression coefficients")
abline(h = 0)

#Statistics calculation
statistics_result = NULL
observed_values = foo1
for (i in seq(1, min(length(foo1)-1,100), by = 1)) {
  predicted_values = model$fitted.values[,1,i]
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_calib = bind_rows(statistics_result, .id = "no_components")
print(r_model_calib)
print(r_model_calib[no_components_onesigma,])
write.csv(r_model_calib, "r_model_calib.csv", row.names = FALSE)

Y1_hat_PLSR_all = model$fitted.values[,1,]
write.csv(Y1_hat_PLSR_all, "Y1_hat_PLSR_all.csv", row.names = FALSE)
Y1_hat_PLSR = as.numeric(Y1_hat[,1,no_components_onesigma])
write.csv(Y1_hat_PLSR, "Y1_hat_PLSR.csv", row.names = FALSE)

##Testing##
statistics_result = NULL
observed_values = Y2
Y2_hat = predict(model, newdata = X2)
for (i in seq(1, min(length(foo1)-4,100), by = 1)) {
  predicted_values = Y2_hat[,1,i] 
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
print(r_model_valid[no_components_onesigma,])
write.csv(r_model_valid, "r_model_valid.csv", row.names = FALSE)

Y2_hat_PLSR_all = Y2_hat[,1,]
write.csv(Y2_hat_PLSR_all, "Y2_hat_PLSR_all.csv", row.names = FALSE)
Y2_hat_PLSR = as.numeric(Y2_hat[,1,no_components_onesigma])
write.csv(Y2_hat_PLSR, "Y2_hat_PLSR.csv", row.names = FALSE)

par(mfrow = c(2,2))
matplot(r_model_calib$no_components,r_model_calib$RMSE, xlab = "No of PLS components", ylab = paste0("RMSD of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_calib$no_components,r_model_calib$R_squared, xlab = "No of PLS components", ylab = paste0("R2 of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$RMSE, xlab = "No of PLS components", ylab = paste0("RMSD of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$R_squared, xlab = "No of PLS components", ylab = paste0("R2 of ", property), type = "h", lty = 1, col = "#FF1A00CC")

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
par(mfrow = c(1,1))
plot(Y2, Y2_hat[,1,no_components_onesigma], xlab = "SOC Observed, %", ylab = "SOC Predicted (PLSR), %", col = colours1[9])
title("PLSR validation")
grid()
abline(0, 1)

#%%%%%%%%%%%%%%%%%#
####Section 2####
#%%%%%%%%%%%%%%%%%#
packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", "moments", "caTools",
                        "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
                        "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
                        "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
                        "installr", "resemble", "prospectr", "magrittr", "doParallel", 
                        "parallel", "foreach", "ggplot2", "tidyr", "pls", "ChemoSpec", "MASS", "Johnson")
sapply(packages_to_install, require, character.only = TRUE)

##FUNCTIONS
# Function to calculate the minimum indices
get_lowest_indices <- function(data, no_neighbor = 5) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  unique_min_indices <- unique(min_indices)
  return(unique_min_indices)
}

# Function to calculate R², RMSE, RPD, and RPIQ
calculate_statistics <- function(observed, predicted) {
  r_squared <- cor(observed, predicted)^2
  rmse <- sqrt(mean((observed - predicted)^2))
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  deviation <- observed - mean_observed
  rpd <- sd(observed) / rmse
  iqr_observed <- IQR(observed)
  rpiq <- iqr_observed / rmse
  bias <- mean_observed - mean_predicted
  mae <- mean(abs(observed - predicted))
  result <- list(R_squared = r_squared,RMSE = rmse,RPD = rpd,RPIQ = rpiq,Bias = bias,MAE = mae)
  return(result)
}

calculate_statistics_result <- function(observed_values1,Y2_hat1, ncomp1) {
  cal_statistics_result <- NULL
  for (i in seq(1, ncomp1, by = 1)) {
    predicted_values1 <- Y2_hat1[, 1, i] 
    cal_statistics_result[[i]] <- calculate_statistics(observed_values1, predicted_values1)
  }
  return(cal_statistics_result)
}

setwd("C:/Refined/Draft_3/")
#OSSL library
OSSL = read_csv("OSSL.csv")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,60:length(foo)]))
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(502, 2500, by = 2))
rm(OSSL,foo)
OSSL = foo1
rm(foo1)

#Data provision
IP = OSSL
#Component selection
property = "SOC"
summary(IP[,property])

IP_X1 = IP$spc[!is.na(IP[,property]),]
IP_Y1 = IP[,property][!is.na(IP[,property])]
summary(IP_Y1)
IP_X_O = IP_X1[order(IP_Y1), ]
IP_Y_O = IP_Y1[order(IP_Y1)]
#For SOC the values were made to selected >0.1
IP_Y1[IP_Y1 < 0.1] <- NA
IP_Y = IP_Y1[!is.na(IP_Y1)]
IP_X = IP_X1[!is.na(IP_Y1),]
summary(IP_Y)
IP_X_O = IP_X[order(IP_Y), ]
IP_Y_O = IP_Y[order(IP_Y)]

#Split 8:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 7:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 6:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 5:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 4:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 3:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 2:3
X1 = IP_X_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

#Split 1:3
X1 = IP_X_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

####1NEAREST NEIGHBOUR####
##Training##
foo = X1
foo1 = Y1

summary(foo1)
foo100 = diff(foo1[-1])
summary(foo100)
print(c((mean(foo100)), max(diff(foo100))))

#PLS dissimilarity with optimal selection of components
opc_sel =  list(method = "opc", value = min(length(foo1)-4,100))
o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
cat("Number of components as selected by Minimum RMSE in the NN algorithm:",dim(o_plsd[["projection"]][["scores"]])[2] ,"\n")

#Evaluations
Y_matrix = as.matrix(foo1)
ev = NULL
ev[["o_plsd"]] = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)

#Tabulating, saving the prediction files and plotting
statistics_result = NULL
for (label in names(ev)) {
  observed_values = ev[[label]]$first_nn[, 1]
  predicted_values = ev[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_calib = bind_rows(statistics_result, .id = "Methods")
print(r_calib)
write.csv(r_calib, "r_calib.csv", row.names = FALSE)

Y2_foo = cbind(observed_values,predicted_values)
write.csv(Y2_foo, "Y1_hat_NN.csv", row.names = FALSE)
rm(Y2_foo)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("SOC Observed, %", "SOC Predicted(1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))])
             title(label)
             grid()
             abline(0, 1)
             #text((max(xy[,1])-10), (max(xy[,2])-10),labels = paste("RMSD:", round(x[[label]]$eval[1],3), "\nR:", round(x[[label]]$eval[2],3), "\nR2:", round(x[[label]]$eval[2]^2,3)),pos = 1,col = "black", cex = 1)
           },
           x = ev)
par(mfrow = c(1, 1))

##Testing##
#PLS dissimilarity with optimal selection of components
o_plsd1 = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)

#Evaluations
ev1 = NULL
observed_values = Y2
ev1[["o_plsd1"]]$first_nn = cbind(observed_values, foo1[apply(o_plsd1$dissimilarity, 2, which.min)])

statistics_result = NULL
for (label in names(ev1)) {
  observed_values = ev1[[label]]$first_nn[, 1]
  predicted_values = ev1[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_valid = bind_rows(statistics_result, .id = "Methods")
print(r_valid)
write.csv(r_valid, "r_valid.csv", row.names = FALSE)

Y2_foo = cbind(observed_values,predicted_values)
write.csv(Y2_foo, "Y2_hat_NN.csv", row.names = FALSE)
rm(Y2_foo)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
p = sapply(names(ev1), 
           FUN = function(x, label, labs = c("SOC Observed, %", "SOC Predicted(1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev1))])
             title(label)
             grid()
             abline(0, 1)
             #text((max(xy[,1])-10), (max(xy[,2])-10),labels = paste("RMSD:", round(x[[label]]$eval[1],3), "\nR:", round(x[[label]]$eval[2],3), "\nR2:", round(x[[label]]$eval[2]^2,3)),pos = 1,col = "black", cex = 1)
           },
           x = ev1)
par(mfrow = c(1, 1))

####2PLSR####
##Training##
foo = X1
foo1 = Y1
model = plsr(foo1 ~ foo, ncomp = min(length(foo1)-4,100), validation = "CV", segments = 10)  
summary(model)
explvar(model) #This is output of Xvar/Xtotvar
plot(RMSEP(model), legendpos = "topright")
plot(model$validation$PRESS[1,], legendpos = "topright")
par(mfrow = c(1,1))
comp_value = RMSEP(model)$val[1,1,]
matplot(comp_value, xlab = "No of PLS components", ylab = "RMSD of Y", type = "h", lty = 1, col = "#FF1A00CC")

##DECIDE on the method to be used for choosing the number of components##
no_components = which.min(comp_value)-1
no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
cat("Number of components as selected by Minimum, onesigma and permutation:", no_components,no_components_onesigma,no_components_permut,"\n")
##DECISION ENDS##

#Plots
colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
plot(foo1, model$fitted.values[,1,no_components_onesigma], xlab = "SOC Observed, %", ylab = "SOC Predicted (PLSR), %", col = colours1[9])
title("PLSR calibration")
grid()
abline(0, 1)
#text(90, 30, labels = paste("N_comp:", no_components, "\nRMSD:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$RMSE,3), "\nR2:", round(calculate_statistics(foo1,model$fitted.values[,1,no_components])$R_squared,3)),pos = 1,col = "black", cex = 1.2)

plot(model, plottype = "scores", comps = 1:3)
plot(model, plottype = "loadings", comps = 1:3, legendpos = "topright",
     labels = "numbers", xlab = "Wavelength (nm)", ylab= "Loading value")
abline(h = 0)
plot(model, plottype = "coef", comps = 1:3, legendpos = "topright",
     labels = "numbers", xlab = "Wavelength (nm)", ylab= "Regression coefficients")
abline(h = 0)

#Statistics calculation
statistics_result = NULL
observed_values = foo1
for (i in seq(1, min(length(foo1)-1,100), by = 1)) {
  predicted_values = model$fitted.values[,1,i]
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_calib = bind_rows(statistics_result, .id = "no_components")
print(r_model_calib)
print(r_model_calib[no_components_onesigma,])
write.csv(r_model_calib, "r_model_calib.csv", row.names = FALSE)

Y1_hat_PLSR_all = model$fitted.values[,1,]
write.csv(Y1_hat_PLSR_all, "Y1_hat_PLSR_all.csv", row.names = FALSE)
Y1_hat_PLSR = as.numeric(model$fitted.values[,1,no_components_onesigma])
write.csv(Y1_hat_PLSR, "Y1_hat_PLSR.csv", row.names = FALSE)

##Testing##
statistics_result = NULL
observed_values = Y2
Y2_hat = predict(model, newdata = X2)
for (i in seq(1, min(length(foo1)-4,100), by = 1)) {
  predicted_values = Y2_hat[,1,i] 
  statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
}
r_model_valid = bind_rows(statistics_result, .id = "no_components")
print(r_model_valid)
print(r_model_valid[no_components_onesigma,])
write.csv(r_model_valid, "r_model_valid.csv", row.names = FALSE)

Y2_hat_PLSR_all = Y2_hat[,1,]
write.csv(Y2_hat_PLSR_all, "Y2_hat_PLSR_all.csv", row.names = FALSE)
Y2_hat_PLSR = as.numeric(Y2_hat[,1,no_components_onesigma])
write.csv(Y2_hat_PLSR, "Y2_hat_PLSR.csv", row.names = FALSE)

par(mfrow = c(2,2))
matplot(r_model_calib$no_components,r_model_calib$RMSE, xlab = "No of PLS components", ylab = paste0("RMSD of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_calib$no_components,r_model_calib$R_squared, xlab = "No of PLS components", ylab = paste0("R2 of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$RMSE, xlab = "No of PLS components", ylab = paste0("RMSD of ", property), type = "h", lty = 1, col = "#FF1A00CC")
matplot(r_model_valid$no_components,r_model_valid$R_squared, xlab = "No of PLS components", ylab = paste0("R2 of ", property), type = "h", lty = 1, col = "#FF1A00CC")

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
par(mfrow = c(1,1))
plot(Y2, Y2_hat[,1,no_components_onesigma], xlab = "SOC Observed, %", ylab = "SOC Predicted (PLSR), %", col = colours1[9])
title("PLSR validation")
grid()
abline(0, 1)

#%%%%%%%%%%%%%%%%%#
####Section 3####
#%%%%%%%%%%%%%%%%%#
##3ITERATIONS##
##Here we check the variations using the randomly bootstraping 8:3 split.
packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", "moments", "caTools",
                        "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
                        "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
                        "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
                        "installr", "resemble", "prospectr", "magrittr", "doParallel", 
                        "parallel", "foreach", "ggplot2", "tidyr", "pls", "ChemoSpec", "MASS", "Johnson")
sapply(packages_to_install, require, character.only = TRUE)

##FUNCTIONS
# Function to calculate the minimum indices
get_lowest_indices <- function(data, no_neighbor = 5) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  unique_min_indices <- unique(min_indices)
  return(unique_min_indices)
}

# Function to calculate R², RMSE, RPD, and RPIQ
calculate_statistics <- function(observed, predicted) {
  r_squared <- cor(observed, predicted)^2
  rmse <- sqrt(mean((observed - predicted)^2))
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  deviation <- observed - mean_observed
  rpd <- sd(observed) / rmse
  iqr_observed <- IQR(observed)
  rpiq <- iqr_observed / rmse
  bias <- mean_observed - mean_predicted
  mae <- mean(abs(observed - predicted))
  result <- list(R_squared = r_squared,RMSE = rmse,RPD = rpd,RPIQ = rpiq,Bias = bias,MAE = mae)
  return(result)
}

calculate_statistics_result <- function(observed_values1,Y2_hat1, ncomp1) {
  cal_statistics_result <- NULL
  for (i in seq(1, ncomp1, by = 1)) {
    predicted_values1 <- Y2_hat1[, 1, i] 
    cal_statistics_result[[i]] <- calculate_statistics(observed_values1, predicted_values1)
  }
  return(cal_statistics_result)
}

setwd("C:/Refined/Draft_3/")
#OSSL library
OSSL = read_csv("OSSL.csv")
setwd("C:/Refined/Draft_3/NEW1/")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,60:length(foo)]))
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(502, 2500, by = 2))
rm(OSSL,foo)
OSSL = foo1
rm(foo1)

#Data provision
IP = OSSL
#Component selection
property = "SOC"
summary(IP[,property])

IP_X1 = IP$spc[!is.na(IP[,property]),]
IP_Y1 = IP[,property][!is.na(IP[,property])]
summary(IP_Y1)
IP_X_O = IP_X1[order(IP_Y1), ]
IP_Y_O = IP_Y1[order(IP_Y1)]
#For SOC the values were made to selected >0.1
IP_Y1[IP_Y1 < 0.1] <- NA
IP_Y = IP_Y1[!is.na(IP_Y1)]
IP_X = IP_X1[!is.na(IP_Y1),]
summary(IP_Y)
IP_X_O = IP_X[order(IP_Y), ]
IP_Y_O = IP_Y[order(IP_Y)]

##3Random cal val split at 8:3 ##
set.seed(123)
results_NN = matrix(ncol = 14, nrow = 100)
colnames(results_NN) = c("Iteration","no_components","R_squared_cal","RMSE_cal","RPD_cal","RPIQ_cal","Bias_cal","MAE_cal","R_squared","RMSE","RPD","RPIQ","Bias","MAE")
results_PLSR = matrix(ncol = 16, nrow = 100)
colnames(results_PLSR) = c("Iteration","no_components","no_components_one_sigma","no_components_permutation","R_squared_cal","RMSE_cal","RPD_cal","RPIQ_cal","Bias_cal","MAE_cal","R_squared","RMSE","RPD","RPIQ","Bias","MAE")

start_time = proc.time()
for (m in seq(55)) {
  split = sample.split(IP_Y_O, SplitRatio = 8/11)
  X1 = subset(IP_X_O, split == TRUE)
  X2 = subset(IP_X_O, split == FALSE)
  Y1 = subset(IP_Y_O, split == TRUE)
  Y2 = subset(IP_Y_O, split == FALSE)
  if (m < 51) {
    print(m)
  } else {
    ##1NEAREST NEIGHBOUR##
    ##Training##
    foo = X1
    foo1 = Y1
    opc_sel =  list(method = "opc", value = min(length(foo1)-4,100))
    o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
    cat("Iteration:", m,"Number of components as selected by Minimum RMSE in the NN algorithm:", dim(o_plsd[["projection"]][["scores"]])[2], "\n")
    #Evaluations
    Y_matrix = as.matrix(foo1)
    ev = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)
    observed_values = ev$first_nn[, 1]
    predicted_values = ev$first_nn[, 2]
    statistics_cal = calculate_statistics(observed_values, predicted_values)
    no_components = dim(o_plsd[["projection"]][["scores"]])[2]
    ##Testing##
    #PLS dissimilarity with optimal selection of components
    o_plsd1 = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
    #Evaluations
    observed_values = Y2
    predicted_values = foo1[apply(o_plsd1$dissimilarity, 2, which.min)]
    statistics_val = calculate_statistics(observed_values, predicted_values)
    results_NN[m,] = as.numeric(c(m, no_components, statistics_cal, statistics_val))
    ##2PLSR##
    ##Training##
    model = plsr(foo1 ~ foo, ncomp = min(length(foo1)-4,100), validation = "CV", segments = 10)  
    ##DECIDE on the method to be used for choosing the number of components##
    comp_value = RMSEP(model)$val[1,1,]
    no_components = which.min(comp_value)-1
    no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
    no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
    cat("Number of components as selected by Minimum, onesigma and permutation:", no_components,no_components_onesigma,no_components_permut,"\n")
    ##DECISION ENDS##
    #Statistics calculation
    statistics_result = NULL
    observed_values = foo1
    for (i in seq(1, min(length(foo1)-1,100), by = 1)) {
      predicted_values = model$fitted.values[,1,i]
      statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
    }
    r_model_calib = bind_rows(statistics_result, .id = "no_components")
    statistics_cal = r_model_calib[no_components_onesigma,]
    ##Testing##
    statistics_result = NULL
    observed_values = Y2
    Y2_hat = predict(model, newdata = X2)
    for (i in seq(1, min(length(foo1)-4,100), by = 1)) {
      predicted_values = Y2_hat[,1,i] 
      statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
    }
    r_model_valid = bind_rows(statistics_result, .id = "no_components")
    statistics_val = r_model_valid[no_components_onesigma,]
    results_PLSR[m,] = as.numeric(c(m,no_components,no_components_onesigma,no_components_permut,statistics_cal[2:7],statistics_val[2:7]))
    print(m)
    flush.console()
  }
  }
end_time = proc.time()
user_time = end_time[1] - start_time[1]
cat("Elapsed time:", user_time/3600, "hours \n")

write.csv(results_NN, "results_NN_60.csv", row.names = FALSE)
write.csv(results_PLSR, "results_PLSR_60.csv", row.names = FALSE)

#%%%%%%%%%%%%%%%%%#
####Section 4####
#%%%%%%%%%%%%%%%%%#
##4Correlation between errors of NN and PLSR methods##
##Getting the errors fromm the two methods: Y2_hat_NN, Y2_hat_PLSR and Y2
#Getting predictions for 8:3 splts
setwd("C:/Refined/Draft_3/Results/3Errors/")
e_Y2_hat_pcad = Y2 - Y2_hat_pcad
e_Y2_hat_plsd = Y2 - Y2_hat_plsd
e_Y2_hat_o_pcad = Y2 - Y2_hat_o_pcad
e_Y2_hat_o_plsd = Y2 - Y2_hat_o_plsd
e_Y2_hat_cd = Y2 - Y2_hat_cd
e_Y2_hat_mcd = Y2 - Y2_hat_mcd
e_Y2_hat_ed = Y2 - Y2_hat_ed
e_Y2_hat_cosd = Y2 - Y2_hat_cosd
e_Y2_hat_sinfd = Y2 - Y2_hat_sinfd
e_Y2_hat_PLSR = Y2 - Y2_hat_PLSR 

data = data.frame(value = c(e_Y2_hat_pcad, e_Y2_hat_plsd, e_Y2_hat_o_pcad, e_Y2_hat_o_plsd, e_Y2_hat_cd, e_Y2_hat_mcd, e_Y2_hat_ed, e_Y2_hat_cosd, e_Y2_hat_sinfd, e_Y2_hat_PLSR),
                  Method = rep(c("pcad", "plsd", "o_pcad", "o_plsd", "cd", "mcd", "ed", "cosd", "sinfd", "PLSR"), each = length(e_Y2_hat_pcad)))
error_data = cbind(e_Y2_hat_pcad, e_Y2_hat_plsd, e_Y2_hat_o_pcad, e_Y2_hat_o_plsd, e_Y2_hat_cd, e_Y2_hat_mcd, e_Y2_hat_ed, e_Y2_hat_cosd, e_Y2_hat_sinfd, e_Y2_hat_PLSR)
error_data = as.data.frame(error_data)
colnames(error_data) = c("pcad", "plsd", "o_pcad", "o_plsd", "cd", "mcd", "ed", "cosd", "sinfd", "PLSR")
ggpairs(error_data)


#For 8_3 PLSR v/s NN
error_NN = Y2 - Y2_hat_NN
error_PLSR = Y2 - Y2_hat_PLSR
error_cor = cor(error_NN, error_PLSR)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
par(mfrow = c(1,1))
plot(error_NN, error_PLSR, xlab = "error_NN, %", ylab = "error_PLSR, %", col = colours1[9])
title("Error correlation analysis")
grid()
abline(0, 1)

Y2_foo = cbind(Y2, Y2_hat_NN, Y2_hat_PLSR)
write.csv(Y2_foo, "Y2_hat.csv", row.names = FALSE)
rm(Y2_foo)

error_dataframe = as.data.frame(cbind(error_NN,error_PLSR))
#Plot the errors
ggplot(error_dataframe) +
  geom_density(alpha = 0.5, aes(x = error_NN, fill = "NN")) +
  geom_density(alpha = 0.5, aes(x = error_PLSR, fill = "PLSR")) +
  labs(title = "PDF of Errors from two methods",x = "Error, %",y = "Density",fill = "Method") +
  theme(panel.grid.major = element_line(),panel.grid.minor = element_line())

ggplot(error_dataframe) +
  geom_density(alpha = 0.5, aes(x = error_NN, fill = "NN")) +
  geom_density(alpha = 0.5, aes(x = error_PLSR, fill = "PLSR")) +
  xlim(-15, 15) +
  labs(title = "PDF of Errors from two methods",x = "Error, %",y = "Density",fill = "Method") +
  theme(panel.grid.major = element_line(),panel.grid.minor = element_line())

# Create a data frame
data = data.frame(value = c(error_NN,error_PLSR),Method = rep(c("NN", "PLSR"), each = length(error_NN)))
ggplot(data, aes(x = Method, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  ylim(-6, 6) +
  scale_fill_manual(values = c("NN" = "#F8766D", "PLSR" = "#00BFC4")) +
  labs(title = "Boxplot of Two Vectors",x = "Method",y = "Error, %") +
  theme(panel.grid.major = element_line(),panel.grid.minor = element_line())

rm(error_dataframe,data)


#%%%%%%%%%%%%%%%%%#
####Section 5####
#%%%%%%%%%%%%%%%%%#
##5ITERATIONS##
##Here we check the variations using the stratified 8:3 split and then randomly bootstraping from them to get 3:3 split.
packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", "moments", "caTools",
                        "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
                        "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
                        "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
                        "installr", "resemble", "prospectr", "magrittr", "doParallel", 
                        "parallel", "foreach", "ggplot2", "tidyr", "pls", "ChemoSpec", "MASS", "Johnson")
sapply(packages_to_install, require, character.only = TRUE)

##FUNCTIONS
# Function to calculate the minimum indices
get_lowest_indices <- function(data, no_neighbor = 5) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  unique_min_indices <- unique(min_indices)
  return(unique_min_indices)
}

# Function to calculate R², RMSE, RPD, and RPIQ
calculate_statistics <- function(observed, predicted) {
  r_squared <- cor(observed, predicted)^2
  rmse <- sqrt(mean((observed - predicted)^2))
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  deviation <- observed - mean_observed
  rpd <- sd(observed) / rmse
  iqr_observed <- IQR(observed)
  rpiq <- iqr_observed / rmse
  bias <- mean_observed - mean_predicted
  mae <- mean(abs(observed - predicted))
  result <- list(R_squared = r_squared,RMSE = rmse,RPD = rpd,RPIQ = rpiq,Bias = bias,MAE = mae)
  return(result)
}

calculate_statistics_result <- function(observed_values1,Y2_hat1, ncomp1) {
  cal_statistics_result <- NULL
  for (i in seq(1, ncomp1, by = 1)) {
    predicted_values1 <- Y2_hat1[, 1, i] 
    cal_statistics_result[[i]] <- calculate_statistics(observed_values1, predicted_values1)
  }
  return(cal_statistics_result)
}

setwd("C:/Refined/Draft_3/")
#OSSL library
OSSL = read_csv("OSSL.csv")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,60:length(foo)]))
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(502, 2500, by = 2))
rm(OSSL,foo)
OSSL = foo1
rm(foo1)

#Data provision
IP = OSSL
#Component selection
property = "SOC"
summary(IP[,property])

IP_X1 = IP$spc[!is.na(IP[,property]),]
IP_Y1 = IP[,property][!is.na(IP[,property])]
summary(IP_Y1)
IP_X_O = IP_X1[order(IP_Y1), ]
IP_Y_O = IP_Y1[order(IP_Y1)]
#For SOC the values were made to selected >0.1
IP_Y1[IP_Y1 < 0.1] <- NA
IP_Y = IP_Y1[!is.na(IP_Y1)]
IP_X = IP_X1[!is.na(IP_Y1),]
summary(IP_Y)
IP_X_O = IP_X[order(IP_Y), ]
IP_Y_O = IP_Y[order(IP_Y)]

##Stratified cal val split at 8:3 ##
#Split 8:3
X0 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), ]
Y0 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

set.seed(123)
results_NN = matrix(ncol = 14, nrow = 100)
colnames(results_NN) = c("Iteration","no_components","R_squared_cal","RMSE_cal","RPD_cal","RPIQ_cal","Bias_cal","MAE_cal","R_squared","RMSE","RPD","RPIQ","Bias","MAE")
results_PLSR = matrix(ncol = 16, nrow = 100)
colnames(results_PLSR) = c("Iteration","no_components","no_components_one_sigma","no_components_permutation","R_squared_cal","RMSE_cal","RPD_cal","RPIQ_cal","Bias_cal","MAE_cal","R_squared","RMSE","RPD","RPIQ","Bias","MAE")
Pred_NN = matrix(ncol = 100, nrow = 17268)
Pred_PLSR = matrix(ncol = 100, nrow = 17268)

start_time = proc.time()
for (m in seq(20)) {
  split = sample.split(Y0, SplitRatio = 17271)
  X1 = subset(X0, split == TRUE)
  Y1 = subset(Y0, split == TRUE)
  if (m>0) {
    ##1NEAREST NEIGHBOUR##
    ##Training##
    foo = X1
    foo1 = Y1
    opc_sel =  list(method = "opc", value = min(length(foo1)-4,100))
    o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
    cat("Iteration:", m,"Number of components as selected by Minimum RMSE in the NN algorithm:", dim(o_plsd[["projection"]][["scores"]])[2], "\n")
    #Evaluations
    Y_matrix = as.matrix(foo1)
    ev = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)
    observed_values = ev$first_nn[, 1]
    predicted_values = ev$first_nn[, 2]
    statistics_cal = calculate_statistics(observed_values, predicted_values)
    no_components = dim(o_plsd[["projection"]][["scores"]])[2]
    ##Testing##
    #PLS dissimilarity with optimal selection of components
    o_plsd1 = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
    #Evaluations
    observed_values = Y2
    predicted_values = foo1[apply(o_plsd1$dissimilarity, 2, which.min)]
    statistics_val = calculate_statistics(observed_values, predicted_values)
    results_NN[m,] = as.numeric(c(m, no_components, statistics_cal, statistics_val))
    Pred_NN[,m] = predicted_values
    ##2PLSR##
    ##Training##
    model = plsr(foo1 ~ foo, ncomp = min(length(foo1)-4,100), validation = "CV", segments = 10)  
    ##DECIDE on the method to be used for choosing the number of components##
    comp_value = RMSEP(model)$val[1,1,]
    no_components = which.min(comp_value)-1
    no_components_onesigma = selectNcomp(model, method = "onesigma", plot = TRUE)
    no_components_permut = selectNcomp(model, method = "randomization", plot = TRUE)
    cat("Number of components as selected by Minimum, onesigma and permutation:", no_components,no_components_onesigma,no_components_permut,"\n")
    ##DECISION ENDS##
    #Statistics calculation
    statistics_result = NULL
    observed_values = foo1
    for (i in seq(1, min(length(foo1)-1,100), by = 1)) {
      predicted_values = model$fitted.values[,1,i]
      statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
    }
    r_model_calib = bind_rows(statistics_result, .id = "no_components")
    statistics_cal = r_model_calib[no_components_onesigma,]
    ##Testing##
    statistics_result = NULL
    observed_values = Y2
    Y2_hat = predict(model, newdata = X2)
    for (i in seq(1, min(length(foo1)-4,100), by = 1)) {
      predicted_values = Y2_hat[,1,i] 
      statistics_result[[i]] = calculate_statistics(observed_values, predicted_values)
    }
    r_model_valid = bind_rows(statistics_result, .id = "no_components")
    statistics_val = r_model_valid[no_components_onesigma,]
    results_PLSR[m,] = as.numeric(c(m,no_components,no_components_onesigma,no_components_permut,statistics_cal[2:7],statistics_val[2:7]))
    Pred_PLSR[,m] = Y2_hat[,1,no_components_onesigma]
    print(m)
    flush.console()
  }
  else print(m)
}
end_time = proc.time()
user_time = end_time[1] - start_time[1]
cat("Elapsed time:", user_time/3600, "hours \n")

write.csv(results_NN, "results_NN_20.csv", row.names = FALSE)
write.csv(results_PLSR, "results_PLSR_20.csv", row.names = FALSE)
write.csv(Pred_NN, "Pred_NN_20.csv", row.names = FALSE)
write.csv(Pred_PLSR, "Pred_PLSR_20.csv", row.names = FALSE)
#write.csv(Y2, "Y2.csv", row.names = FALSE)

#%%%%%%%%%%%%%%%%%#
####Section 6####
#%%%%%%%%%%%%%%%%%#
##6ROUGH##
optimal_sel =  list(method = "opc", value = min(length(foo1)-4,40))
pls_tr_opc2 = ortho_projection(Xr = foo, Yr = foo1, method = "pls", pc_selection = optimal_sel)
pls_tr_opc1 = ortho_projection(Xr = foo, Yr = foo1, Xu = X2, Yu = Y2, method = "pls",pc_selection = optimal_sel, scale = TRUE)
length(pls_tr_opc[["scores"]][1:10,1])
length(pls_tr_opc1[["scores"]][1:10,1])

pls_final_model <- plsr(foo1 ~ foo, ncomp = 22)
summary(pls_final_model)
foo_X2_plsr = predict(pls_final_model, newdata = X2)

foo100 = pls_tr_opc[["scores"]]
foo_proj_matrix = pls_tr_opc[["projection_mat"]]
foo_X2_projected = X2 %*% foo_proj_matrix

get_all_indices <- function(data, no_neighbor = 1) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  return(min_indices)
}
all_indices = get_all_indices(o_plsd$dissimilarity, 1)
all_indices[1]
o_plsd$dissimilarity[153,1]

f1 = as.matrix(foo100[153,])
f2 = as.matrix(foo_X2_plsr[1,1,])
md_foo100_f2 = f_diss(foo100,f2, "mahalanobis")
rm(foo100,f1,f2)

# Final model with the optimal number of components
pls_final_model <- plsr(y ~ X, ncomp = ncomp_opt)
# Summary of the final model
summary(pls_final_model)
#Check return projection ?
opc_sel =  list(method = "opc", value = min(length(foo1)-4,100))
#PLS dissimilarity with optimal selection of components
o_plsd_foo = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)
o_plsd_foo = dissimilarity(Xr=foo[1:1000,], Yr = foo1[1:1000], Xu = X2[1:100,], Yu = Y2[1:100], diss_method = "pls",pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)

data(NIRsoil)
Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
# Euclidean distances between all the observations in Xr
ed100 <- f_diss(Xr = Xr, diss_method = "euclid")
# Euclidean distances between observations in Xr and observations in Xu
ed101 <- f_diss(Xr, Xu)
# Mahalanobis distance computed on the first 20 spectral variables
md_xr_xu <- f_diss(Xr[, 1:20], Xu[, 1:20], "mahalanobis")
dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE, return_projection = TRUE)

#Plots
colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan", "brown")
colours1 = c("black", "black","black","black","black","black","black","black","black")
par(mfrow = c(4, 3))
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("SOC, %", "SOC (NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))], cex.lab = 1.2, cex.axis = 1.2)
             title(label, cex.main = 1.5)
             grid()
             abline(0, 1)
             text(
               (max(xy[,1])-10), (max(xy[,2])-10),
               labels = paste("RMSE:", round(statistics_result[[label]]$RMSE,3), "\nR2:", round(statistics_result[[label]]$R_squared,3)),
               pos = 1,
               col = "black", cex = 1.2
             )
           },
           x = ev)
plot(Y2, Y2_hat[,1,no_components_onesigma], xlab = "SOC, %", ylab = "SOC (PLSR), %", col = colours1[10], cex.lab = 1.2, cex.axis = 1.2)
title("plsr", cex.main = 1.5)
grid()
abline(0, 1)
text((max(Y2)-10), (max(Y2)-10), labels = paste("N_comp:", no_components_onesigma, "\nRMSE:", round(calculate_statistics(Y2,Y2_hat[,1,no_components_onesigma])$RMSE,3), "\nR2:", round(calculate_statistics(Y2,Y2_hat[,1,no_components_onesigma])$R_squared,3)),
     pos = 1,col = "black", cex = 1.2)
par(mfrow = c(1, 1))
rm(foo,foo1)



#%%%%%%%%%%%%%%%%%#
#Errors
#%%%%%%%%%%%%%%%%%#
##Correlation between NN and PLSR methods##
##Getting the errors fromm the two methods: Y2_hat_NN, Y2_hat_PLSR and Y2
#For 8_3 PLSR v/s NN
error_NN = Y2 - Y2_hat_NN
error_PLSR = Y2 - Y2_hat_PLSR
error_cor = cor(error_NN, error_PLSR)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
colours1 = c("black", "black","black","black","black","black","black","black","black")
par(mfrow = c(1,1))
plot(error_NN, error_PLSR, xlab = "error_NN, %", ylab = "error_PLSR, %", col = colours1[9])
title("Error correlation analysis")
grid()
abline(0, 1)

Y2_foo = cbind(Y2, Y2_hat_NN, Y2_hat_PLSR)
write.csv(Y2_foo, "Y2_hat.csv", row.names = FALSE)
rm(Y2_foo)

error_dataframe = as.data.frame(cbind(error_NN,error_PLSR))
#Plot the errors
ggplot(error_dataframe) +
  geom_density(alpha = 0.5, aes(x = error_NN, fill = "NN")) +
  geom_density(alpha = 0.5, aes(x = error_PLSR, fill = "PLSR")) +
  labs(title = "PDF of Errors from two methods",x = "Error, %",y = "Density",fill = "Method") +
  theme(panel.grid.major = element_line(),panel.grid.minor = element_line())

ggplot(error_dataframe) +
  geom_density(alpha = 0.5, aes(x = error_NN, fill = "NN")) +
  geom_density(alpha = 0.5, aes(x = error_PLSR, fill = "PLSR")) +
  xlim(-15, 15) +
  labs(title = "PDF of Errors from two methods",x = "Error, %",y = "Density",fill = "Method") +
  theme(panel.grid.major = element_line(),panel.grid.minor = element_line())

# Create a data frame
data = data.frame(value = c(error_NN,error_PLSR),Method = rep(c("NN", "PLSR"), each = length(error_NN)))
ggplot(data, aes(x = Method, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  ylim(-6, 6) +
  scale_fill_manual(values = c("NN" = "#F8766D", "PLSR" = "#00BFC4")) +
  labs(title = "Boxplot of Two Vectors",x = "Method",y = "Error, %") +
  theme(panel.grid.major = element_line(),panel.grid.minor = element_line())

rm(error_dataframe,data)