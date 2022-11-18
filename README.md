# StemAnalysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/forestscientist/StemAnalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/forestscientist/StemAnalysis/actions/workflows/R-CMD-check.yaml)
[![Travis build status](https://travis-ci.com/forestscientist/StemAnalysis.svg?branch=main)](https://travis-ci.com/forestscientist/StemAnalysis)
<!-- badges: end -->

## Abstract
#### StemAnalysis R package is a tool for designed to reconstruct stem growth profiles, construct height-diameter relationships, and consequently to compute growth trends in terms of diameter at breast height (DBH), tree height, stem volume, tree biomass and carbon storage for an individual tree. This vignette provides an overview of this package functions and options. We provide a working examples that demonstrates the basic functionality and use of the package.

## Purpose
#### Accurate information about age dynamics of timber production and carbon storage in forest ecosystems is frequently required by scientists, stakeholders, and policymakers. Stem analysis is a technique for measuring tree growth (Salas-Eljatib, 2021). The computational burden of reconstructing  temporal, radial, and longitudinal patterns of tree growth, fitting height-diameter relationships, and calculating diameter with bark from radial annual-ring increment sequences measured on multiple cross-sectional discs, may present a hindrance to application of stem analysis methodology in forest research investigations and operational forest multifunctional management (Newton, 2019). Therefore, a standardized tool, StemAnalysis R package, is developed to calculate tree growth dynamics and then make the stem analysis technique more conveniently applied to forest multifunctional investigation.


# 1. Installation
#### To install the current (development) version from the repository, run the following command:

```{r, eval=FALSE}
if(!require(devtools)){install.packages(devtools)}
devtools::install_github(repo = "forestscientist/StemAnalysis")
```

# 2. Load the package

```{r, eval=FALSE}
library(StemAnalysis)
```

# 3. Load the stem analysis data stored in the package
#### d_stem is a dataset containing the input data of stem analysis. Note: If a user uses the StemAnalysis package to analysis a very big tree, the number of inner growth rings that diameter measured for some cross-sectional discs may be more than 11, the Dnobark12, Dnobark13, and much more variables can be added, which also could successfully run.

```{r, eval=TRUE, cache=TRUE}
data(d_stem)
str(d_stem)
```
### The list of input variables and their description in d_stem dataset
#### No: The line number
#### Treeno: The tree number for the sampled tree, the same number represents the same tree
#### TreeTH: Tree total height (m)
#### stemheight: The stem height that the cross-sectional discs were obtained (m). The detail information see Figure 1
#### stemage: The age at the stem height, namely the number of growth rings of the cross-sectional disc (year)
#### Dwithbark: The maximum diameter of the cross-sectional disc with bark (cm)
#### Dnobark0: The maximum diameter of the cross-sectional disc without bark (cm)
#### Dnobark1: The diameter for the 1th age class inner growth ring of the cross-sectional disc (cm). The detail information see Figure 2
#### Dnobark2: The diameter for the 2th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark3: The diameter for the 3th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark4: The diameter for the 4th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark5: The diameter for the 5th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark6: The diameter for the 6th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark7: The diameter for the 7th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark8: The diameter for the 8th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark9: The diameter for the 9th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark10: The diameter for the 10th age class inner growth ring of the cross-sectional disc (cm)
#### Dnobark11: The diameter for the 11th age class inner growth ring of the cross-sectional disc (cm)

![image](https://github.com/forestscientist/StemAnalysis/blob/main/man/Figures/stemheight.png)

#### Figure 1 Two sampling scenarios of cross-sectional discs obtained from individual trees. Scenario 1, for the tree height lower than 10 m, the cross-sectional discs were obtained at the height of 0 m (ground level), 0.5 m, 1.3 m, 1.5 m, and at 1-m intervals thereafter; scenario 2, for the tree height higher than 10 m, were collected at the height of 0 m, 1.3 m, 3.6 m, and at 2-m intervals thereafter. Tree age of the example tree for scenario 1 is 12-years-old and for scenario 2 is 31-years-old.

![image](https://github.com/forestscientist/StemaAnalysis/blob/main/man/Figures/measuremethod.png)

#### Figure 2 A schematic diagram for the measurement of a given cross-sectional disc. The age class defined as the tolerance of the arithmetic sequence of the maximum stem age. In fact, the age class can be 1 year, 2 years, or even 5 years, depending on the diameter measured for the inner growth rings.

# 4. Load the parameter data stored in the package
#### Tree biomass estimated using allometric models [ln(Bi)=a+b×ln(DBH^2H)] (Xiang et al., 2021). The parameter data, including parameters a, b, and carbon concentration for each tissues, are optional inputs.

```{r, eval=TRUE, cache=TRUE}
data(d_parameters)
str(d_parameters)
```
### The list of input variables and their description in d_parameters dataset
#### tissues: The tree tissues including stem, branch, leaf, root, and total tree
#### a: The parameter a in the allometric model ln(Bi)=a+b×ln(DBH2H)
#### b: The parameter b in the allometric model ln(Bi)=a+b×ln(DBH2H)
#### Cconcentration: The tree carbon concentration (kg C kg-1)

# 5. Load the biomass expansion factor data stored in the package
#### Total tree biomass estimated using volume model [V*WD*BEF*(1+R)] (IPCC, 2003). The biomass expansion factor data, including wood density (WD), biomass expansion factor (BEF), root-shoot ratio (R), and C concentration, are optional inputs.

```{r, eval=TRUE, cache=TRUE}
data(d_BEF)
str(d_BEF)
```
### The list of input variables and their description in d_BEF dataset
#### WD: The wood density (kg m-3)
#### BEF: The Biomass Expasion Factor
#### R: The Root-to-Shoot Ratio in terms of biomass
#### Cconcentration: The tree carbon concentration (kg C kg-1)

# 6. Application of StemAnalysis package
## 6.1 Stem growth analysis
#### Reconstructed stem growth patterns and calculated DBH and tree height growth trends, and stem volume with bark increment trends using stem analysis data.   

```{r, eval=TRUE, cache=TRUE}
stemanalysism(xtree = 8, stemgrowth = TRUE, stemdata = d_stem)
```

![image](https://github.com/forestscientist/StemaAnalysis/blob/main/man/Figures/StemGrowth.png)

#### Figure 3 Stem growth patterns of an individual tree. (a) shows the stem growth pattern; (b), (c), and (d) are the cumulative growth and (e), (f), and (g) are the mean annual increment (red dotted line) and current annual increment (blue dashed line) of DBH, tree height, and wood volume, respectively.

## 6.2 Estimation Of tree carbon accumulation
### 6.2.1 Tree biomass and carbon accumulation estimated by allometric models
#### If set 'treecarbon = TRUE' and provide parameter data, tree biomass and carbon accumulation estimated by allometric models.

```{r, eval=TRUE, cache=TRUE}
stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = d_stem, parameterdata = d_parameters)
```
![image](https://github.com/forestscientist/StemaAnalysis/blob/main/man/Figures/TreeCarbon_allometric.png)

#### Figure 4 The age dynamics of total tree biomass (a) and carbon storage (b) for the 20-year-old tree estimated using allometric models.

### 6.2.2 Tree biomass and carbon accumulation estimated by volume model
#### If set 'treecarbon = TRUE' and provide biomass expansion factor data, tree biomass and carbon accumulation estimated by volume model.

```{r, eval=TRUE, cache=TRUE}
stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = d_stem, BEFdata = d_BEF)
```
![image](https://github.com/forestscientist/StemaAnalysis/blob/main/man/Figures/TreeCarbon_volume.png)

#### Figure 5 The age dynamics of total tree biomass (a) and carbon storage (b) for the 20-year-old tree estimated using volume model.

## 6.3 Construction of height-diameter relationship
#### If set 'HDmodel = TRUE', tree height-diameter relationship will be constructed by nonlinear models, and the fitted statistics are showed in a graph.

```{r, eval=FALSE, cache=TRUE}
stemanalysism(xtree = 8, HDmodel = TRUE, stemdata = d_stem)
```
![image](https://github.com/forestscientist/StemaAnalysis/blob/main/man/Figures/TreeCarbon_volume.png)

#### Figure 6 Tree height-diameter relationships for the 20-year-old tree. The fitted curves of the Chapman-Richards model (red line), Logistic model (blue line), Weibull model (green line), and Gomperz model (yellow line) as well as their fitted statistics. a, b and C are the parameters of the nonlinear models; R2 is the coefficient of determination; RMSE is the root mean square error; RSS is the residual sum of squares.

# References
#### IPCC. (2003) Good Practice Guidance for Land Use, Land-Use Change and Forestry; IPCC/IGES: Hayama, Japan.

#### Newton, P.F. (2019) Examining naturogenic processes and anthropogenic influences on tree growth and development via stem analysis: data processing and computational analytics. Forests 10, 1058.

#### Salas-Eljatib, C. (2021) A new algorithm for reconstructing the height growth with stem analysis data. Methods Ecol. Evol. 12, 2008–2016.

#### Xiang, W.H., Li, L.H., Ouyang, S., Xiao, W.F., Zeng, L.X., Chen, L., Lei, P.F., Deng, X.W., Zeng, Y.L., Fang, J.P. & Forrester, D.I. (2021) Effects of stand age on tree biomass partitioning and allometric equations in Chinese fir (Cunninghamia lanceolata) plantations. Eur. J. For. Res. 140, 317–332.
