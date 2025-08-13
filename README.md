# SKR_TraitAbundanceDistribution

Based on an analysis of the Skewness-Kurtosis Relationship (SKR), this framework promotes the study of the shapes of the Trait Abundance Distributions (TADs) to better understand community assembly processes, and predict community dynamics under environmental changes.

<!-- badges: start -->
[![Pipeline mainbranch](https://forgemia.inra.fr/urep/data_processing/tad/badges/main/pipeline.svg?key_text=Pipeline+main+branch&key_width=130)](https://forgemia.inra.fr/urep/data_processing/tad/pipelines/main/latest)
[![LatestRelease](https://forgemia.inra.fr/urep/data_processing/tad/-/badges/release.svg)](https://forgemia.inra.fr/urep/data_processing/tad/-/releases)
[![Coverage](https://forgemia.inra.fr/urep/data_processing/tad/badges/main/coverage.svg?key_text=Coverage)](https://urep.pages.mia.inra.fr/data_processing/tad/coverage/report.html)
[![Usermanual](https://urep.pages.mia.inra.fr/data_processing/tad/manual.svg)](https://urep.pages.mia.inra.fr/dev_utils/r_utils/r4urep/index.html)
<!-- badges: end -->

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
