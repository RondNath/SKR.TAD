# SKR_TraitAbundanceDistribution

Based on an analysis of the Skewness-Kurtosis Relationship (SKR), this framework promotes the study of the shapes of the Trait Abundance Distributions (TADs) to better understand community assembly processes, and predict community dynamics under environmental changes.

Installing the package
---------------------------------
The version 1.0.0 is available on CRAN, https://cran.r-project.org/web/packages/TAD/index.html, can be installed using:

install.packages("TAD")


The development version on Github can be installed using the devtools package:

devtools::install_github("RondNath/SKR.TAD")


Step 1 - Function: Randomize the abundances
---------------------------------
A “randomization” consists in swapping the lines between them (i.e. for community ecology, thus keeping the abundance structure and the number of species constant). This randomization must be carried out n times.

The randomization can be done by respecting a series of categorical factors (i.e. by year, fertility level, …)


Step 2 - Function: SKR analysis of the TADs & Graphical representation
----------------------------------------------------
-   Calculate the different moments (mean, variance, skewness, kurtosis) for all abundance matrices (observed and n randomized)
-   Calculate deviations from the n randomization (null model) via a Confidence Interval (default: CI = 0.95%)
-   It is allowed to change the significance thresholds
-   Calculate the SKR for the observations and the null model
-   Extract the parameters
    -   R²
    -   TADstab (RMSE)
    -   intercept (alpha)
    -   slope (beta)
    -   distance to a reference distribution family (default TADeve, skew-uniform family (slope = 1; intercept = 1.86))
    -   conditional probabilities
    -   Null envelope
-   Get the confidence interval for each parameter and compare it to the null envelope

  
Step 3 - Function: Graphical representation of the SKR and parameters
----------------------------------------------------
- Graphical representation of the observed moments relative to random moments
- Graphical representation of the SKR (Kurtosis ~ Skewness²)
- Graphical representation of the observed SKR parameters (linear regression parameters) compared to the randomized parameters
