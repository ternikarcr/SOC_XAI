# SOC_XAI
# Soil Organic Carbon Analysis
This repository contains a collection of scripts and models for analyzing Soil Organic Carbon (SOC) data using Traditional method and Explainable AI methods for SOC prediction using VNIR lab spectroscopy (OSSL). The project explores SOC prediction through explainable AI frameworks.

## Contents

### Scripts

1. **1Dataset\_prep.R**

   - Prepares and preprocesses the SOC dataset for analysis.
   - Features:
     - Reads and reformats data from CSV files.
     - Performs spectral transformation (log10).
     - Splits data into training and testing sets using various ratios.
     - Outputs preprocessed datasets into CSV files.
   - Key Outputs:
     - Train-test datasets in multiple splits (e.g., 8:3, 7:3).

2. **1SOC\_NN.R**

   - Implements nearest neighbor and partial least squares regression (PLSR) for SOC prediction.
   - Features:
     - Calculates evaluation metrics like RMSE, RPD, RPIQ.
     - Compares various distance metrics.
     - Validates the models with calibration and validation statistics.
   - Key Outputs:
     - Model evaluation metrics (e.g., RMSE, R-squared).
     - Predictions for validation datasets.

## Prerequisites

- **Languages:** R
- **Dependencies:**
  - R Packages: `readxl`, `writexl`, `dplyr`, `ggplot2`, `pls`, `reshape2`, `heatmap3`, `foreach`, `parallel`, `ChemoSpec`
  
## How to Use

1. Clone this repository:

   git clone <repository-url>
   
2. Install dependencies:

     install.packages(c("readxl", "writexl", "dplyr", "ggplot2", "pls", "reshape2", "heatmap3", "foreach", "parallel", "ChemoSpec"))
   
3. Run individual scripts as needed:

   - R scripts can be executed in RStudio.
   
## File Structure

- **1Dataset\_prep.R**: Data preprocessing script.
- **1SOC\_NN.R**: Nearest neighbor and PLSR script.


## Acknowledgments

This project leverages data and methods from the field of soil spectroscopy and explainable AI. Special thanks to the open-source libraries and contributors.

---

For any issues or questions, please feel free to raise a GitHub issue or contact the repository maintainer.
