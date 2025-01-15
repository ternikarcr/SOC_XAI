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

setwd("C:/Refined/Draft_6_Ex_AI_Book_Chapter")
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

train = data.frame(Y = Y1, X = X1)
test = data.frame(Y = Y2, X = X2)
write.csv(train, "train_8_3.csv", row.names = FALSE)
write.csv(test, "test_8_3.csv", row.names = FALSE)

#Split 7:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

train = data.frame(Y = Y1, X = X1)
test = data.frame(Y = Y2, X = X2)
write.csv(train, "train_7_3.csv", row.names = FALSE)
write.csv(test, "test_7_3.csv", row.names = FALSE)

#Split 6:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

train = data.frame(Y = Y1, X = X1)
test = data.frame(Y = Y2, X = X2)
write.csv(train, "train_6_3.csv", row.names = FALSE)
write.csv(test, "test_6_3.csv", row.names = FALSE)

#Split 5:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

train = data.frame(Y = Y1, X = X1)
test = data.frame(Y = Y2, X = X2)
write.csv(train, "train_5_3.csv", row.names = FALSE)
write.csv(test, "test_5_3.csv", row.names = FALSE)

#Split 4:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

train = data.frame(Y = Y1, X = X1)
test = data.frame(Y = Y2, X = X2)
write.csv(train, "train_4_3.csv", row.names = FALSE)
write.csv(test, "test_4_3.csv", row.names = FALSE)

#Split 3:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
#rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

train = data.frame(Y = Y1, X = X1)
test = data.frame(Y = Y2, X = X2)
write.csv(train, "train_3_3.csv", row.names = FALSE)
write.csv(test, "test_3_3.csv", row.names = FALSE)

