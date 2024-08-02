![INFORMS Journal on Computing Logo](https://camo.githubusercontent.com/1b8f04b8ff248ffd132c13343858d070c4805406bbd4c4651f9b27e9c2f01a58/68747470733a2f2f494e464f524d534a6f432e6769746875622e696f2f6c6f676f732f494e464f524d535f4a6f75726e616c5f6f6e5f436f6d707574696e675f4865616465722e6a7067)
# Fraud detection by integrating multisource heterogeneous presence-only data
This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](https://github.com/INFORMSJoC/2023.0257/blob/master/LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Fraud detection by integrating multisource heterogeneous presence-only data](https://doi.org/10.1287/ijoc.2023.0366) by Y. Qiu, Y. Chen, K. Fang, L. Yu, and K. Fang.
## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0366

https://doi.org/10.1287/ijoc.2023.0366.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{QiuIPU2024,
  author =        {Y. Qiu, Y. Chen, K. Fang, L. Yu, and K. Fang},
  publisher =     {INFORMS Journal on Computing},
  title =         {Fraud detection by integrating multisource heterogeneous presence-only data},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0366.cd},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0366},
} 
```

## Abstract

In credit fraud detection practice, certain fraudulent transactions often evade detection due to the hidden nature of fraudulent behavior. To address this issue, an increasing number of Positive Unlabelled (PU) learning techniques have been employed by more and more financial institutions. However, most of these methods are designed for single datasets and do not take into account the heterogeneity of data when they are collected from different sources. In this paper, we propose an integrative PU learning method (IPU) for pooling information from multiple heterogeneous PU datasets. A novel approach that penalizes group differences is developed to explicitly and automatically identify the cluster structures of coefficients across different datasets, thus offering a plausible interpretation of heterogeneity. Furthermore, we apply a bi-level selection method to detect the sparse structure at both the group level and within-group level. Theoretically, we show that our proposed estimator has the oracle property. Computationally, we design an expectation-maximization (EM) algorithm framework and propose an alternating direction method of multipliers (ADMM) algorithm to solve it. Simulation results show that our proposed method has better numerical performance in terms of variable selection, parameter estimation, and prediction ability. Finally, a real-world application showcases the effectiveness of our method in identifying distinct coefficient clusters and its superior prediction performance compared to direct data merging or separate modeling. This result also offers valuable insights for financial institutions in developing targeted fraud detection systems.  

## Description

This repository provides data for the problem and code for the method. The main folders are `data`and `code`.

- `data`: This folder includes the data used in the paper.
- `code`: This folder contains the source code and the code for experimental comparison.

## Data

- The raw data used in our study can be downloaded from [Kaggle](https://www.kaggle.com/competitions/ieee-fraud-detection/data). The process of feature engineering is documented on [GitHub](https://github.com/xiaoluoyfy/IEEE-CIS-Fraud-Detection/) (note that variables derived from regions have been removed). The cleaned dataset is available in the data folder.
- A demo dataset is also provided for quick start.

## Replicating

The `code` folder contains all the code which implements the framework of this paper. 

* `I-PU.R` contains functions for building and optimizing the I-PU model.
* Comparison Methods
  - `single model.R` includes several single models, and their estimation results can be selected as initial values for a warm start. To implement the code, you need to install the R packages **glmnet**, **PUlasso**, **ncvreg**, and **grpreg**.
  - `Oracle.R` includes functions for calculating the Oracle model.
  - `I-LR.R` contains functions for calculating the I-LR model.
  - `PYMTL.py` contains functions for calculating the Clustered-LR, Lowrank-LR, and AR-LR models. The code implements the method from Duan, Y. and Wang, K (2023)<sup>[1]</sup>. To ensure the code runs correctly, you need to download the **ARMUL.py** and **MTL.py** files from the GitHub repository at https://github.com/kw2934/armul and install the required Python packages.
  - `RMTL.R` contains the function for calculating the Sparse-LR model<sup>[2]</sup>. To implement the code, you need to install the R package **RMTL**.
* `criteria.R` contains functions for calculating evaluation criteria.
* `data_generate.R` contains the data generators for the three examples used in the simulation.
* `demo.html` contains the code and results of the demo.

## Demo 

To demonstrate the application of our proposed method, we have created a demo. Readers can refer to **demo.html** to review the code and results.

```{r message=FALSE}
library(PUlasso)
source('data_generate.R')
source('I-PU.R')
source('criteria.R')
```

### Data generation

```{r}
# p0: Coefficient dimension without an intercept term.
# example: 1, 2, 3 correspond to the three examples in the simulation.
# l: Random number seed.

p0 = 200
example = 1
l=12345
# Use `data_generate` to create data
Data = data_generate(l, example = example, p0=p0)
# or you can directly read the `demo_data.Rdata`: `load("demo_data.Rdata")`.
beta_true = Data$beta_true
M  = Data$M; n  = Data$n; nu = Data$nu; nl = Data$nl; pi = Data$pi
sample_size = Data$sample_size; N = Data$N; p = Data$p
X_m  = Data$X_m; 
Y = Data$Y; Y_test = Data$Y_test;  Y_valid = Data$Y_valid
X_test = Data$X_test; X_valid = Data$X_valid
Z = Data$Z; Z_test = Data$Z_test; Z_valid = Data$Z_valid
group = Data$group
important_coef = Data$important_coef
```

### Obtain the initial values

Here, we use the PUlasso package to estimate initial values for demonstration; you may also choose other suitable initial values.

```{r warning=FALSE}
beta_pu_grplasso = matrix(0, M, p)
for(m in 1 : M){
  CV_pu_grplasso = cv.grpPUlasso(X = X_m[[m]], z = Z[[m]], py1 = pi[[m]], 
                                 lambda = seq(0, 0.025, 0.001),
                                 group = group[-1])
  lambda_index = which(CV_pu_grplasso$lambda == CV_pu_grplasso$lambda.min)
  beta_pu_grplasso[m, ] = coef(CV_pu_grplasso)[ ,lambda_index]
}
init_beta = beta_pu_grplasso
```

### Model training and parameter selection

```{r}
# Set the sequence for lambda1 and lambda2.
lambda1_seq = c(seq(0.15,0.1,-0.01))
lambda2_seq = c(seq(0.08,0.07,-0.002))
# Train the model, and the result will return the estimated coefficients corresponding to each set of parameters, as well as the optimal parameters selected using the validation set and their corresponding indices.
fit_IPU = I_PU_train(init_beta,X_m,X,Z,X_valid,Y_valid,Z_valid,sample_size,M,pi,N,p,group,lambda1_seq,lambda2_seq,loop=50,a = 3,rho1=1,rho2=1,beta_true = NULL)

```

### Evaluate the selected optimal estimated parameters

```{r}
best_ind = fit_IPU$ind_min
beta_IPU = fit_IPU$beta_IPU_list[[best_ind]]
print(cal_criteria(beta_IPU, beta_true, M, p0, X_test, Y_test, Z_test, pi))
print(cal_criteria(beta_pu_grplasso, beta_true, M, p0, X_test, Y_test, Z_test, pi))
```

## Contact & Support

Please contact Qiu Yongqin via qiuyongqin@ustc.edu.cn for support in using this software or questions about the paper.

## Reference

[1] Duan Y, Wang K (2023) Adaptive and robust multi-task learning. The Annals of Statistics 51(5):2015–2039.  
[2] Cao H, Zhou J, Schwarz E (2019) RMTL: an R library for multi-task learning. Bioinformatics 35(10):1797–1798.  
