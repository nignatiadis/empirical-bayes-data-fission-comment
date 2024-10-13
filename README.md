# empirical-bayes-data-fission-comment

The main scripts are the following: 

* `gaussian_aurora.R`: Simulations for EB estimation with a normal likelihood. Running this script generates two files:
  * "double_gaussian_plot.png" which is panel (a) of Figure 1
  * "double_gaussian_MSE.csv", which contains the data for the column "Gaussian MSE" of Table 1

* `poisson_aurora.R`: Simulations for EB estimation with a Poisson likelihood. Running this script generates two files:
  * "double_poisson_plot.png" which is panel (b) of Figure 1
  * "double_poisson_MSE.csv", which contains the data for the column "Poisson MSE" of Table 1


## Further instructions

Running the above two scripts presupposes a working installation of R and several packages. The following file

* `R_session_info.txt`

records the session information about the version of R used as well as of all required packages.

For the implementation of the nonparametric maximum likelihood estimator (NPMLE), we used the CRAN R package [REBayes](https://cran.r-project.org/web/packages/REBayes/index.html). 
REBayes requires a working installation of the commercial optimization solver [Mosek](https://www.mosek.com/) (we used Mosek version 10.1). Mosek provides [free licenses for academic use](https://www.mosek.com/products/academic-licenses/).
