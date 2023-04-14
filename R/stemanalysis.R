#' Reconstructing Tree Growth and Carbon Accumulation with Stem Analysis Data
#'
#' @param xtree Xtree is the tree number (Treeno), which is used to choose a target tree to be analyzed
#' @param stemgrowth If stemgrowth is 'TRUE', stem growth profile and growth trends in terms of diameter at breast height (DBH), tree height, and stem volume will be showed in a graph. A example graph is man/Figures/StemGrowth.png
#' @param treecarbon If treecarbon is 'TRUE', total tree biomass and carbon storage will be estimated by general allometric models (National Forestry and Grassland Administration, 2014) and volume model (Fang et al., 2001). The example graphs are man/Figures/TreeCarbon_allometric.png and TreeCarbon_volume. In addition, although treecarbon is 'TRUE', the estimation of tree biomass and carbon storage by allometric models will skip if 'allompardata' is missing, and the same is true for the estimation by volume model if 'volumepardata' is missing.
#' @param HDmodel If HDmodel is 'TRUE', height-diameter relationship will be fitted with nonlinear models (Mehtatalo, 2017) and showed the fitted results in a graph. A example graph is man/Figures/HDmodel.png
#' @param stemdata table as described in \code{\link{stemdata}} containing the information about stem analysis data.
#' @param allompardata table as described in \code{\link{allomPardata}} containing the information about the parameter data of allometric models that can be optionally inputted by users
#' @param volumepardata table as described in \code{\link{volumePardata}} containing the information about the biomass conversion factor data for volume model that can be optionally inputted by users
#'
#' @note The \code{stemanalysis} was performed on individual trees
#'
#' @return A list with class "output" containing three data.frame.
#'    - `StemGrowth`: the estimated stem growth trends data for a target tree, including the tree age class and the corresponding growth data of diameter at breast height (DBH), stem height, and stem volume. More details on the output is \code{\link{StemGrowth}}
#'    - `allomCarbon`: the estimated tree biomass and carbon storage data by using allometric models for a target tree, including tree biomass and carbon storage for aboveground, belowground, and total tree. More details on the output is \code{\link{allomCarbon}}
#'    - `volumeCarbon`: the estimated tree biomass and carbon storage data by using volume model for a target tree, including tree biomass and carbon storage for aboveground, belowground, and total tree. More details on the output is \code{\link{volumeCarbon}}
#'
#' @export
#'
#' @references Fang, J., Chen, A., Peng, C., et al. (2001)
#' Changes in forest biomass carbon storage in China between 1949 and 1998.
#' \emph{Science}
#' \bold{292}, 2320-2322. {doi:10.1126/science.1058629}
#'
#' Mehtatalo, L. (2017)
#' Lmfor: Functions for forest biometrics.
#' {https://CRAN.R-project.org/package=lmfor}
#'
#' National Forestry and Grassland Administration. (2014) Tree biomass models and
#' related parameters to carbon accounting for Cunninghamria lanceolata.
#' \emph{Forestry industry standards of the People's Republic of China}
#' {Beijing, LY/T 2264—2014}
#'
#' @examples
#'
#' library(StemAnalysis)
#'
#' # Load the data sets
#' data(stemdata)
#' data(volumePardata)
#' data(allomPardata)
#'
#' # To calculating stem growth trends for an individual tree is needed
#' stemanalysism(xtree = 8, stemgrowth = TRUE, stemdata = stemdata)
#'
#' # To calculating tree carbon storage by allometric models is needed
#' stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = stemdata,
#'              allompardata = allomPardata)
#'
#' # To calculating tree carbon storage by volume model is needed
#' stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = stemdata,
#'              volumepardata = volumePardata)
#'
#' # To fitting the height-diameter relationships
#' stemanalysism(xtree = 8, HDmodel = TRUE, stemdata = stemdata)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics axis layout layout.show legend lines mtext par
#' @importFrom stats anova lm nls nls.control predict residuals
#' @importFrom utils write.table
#' @importFrom stats AIC logLik

stemanalysism <- function(xtree, stemgrowth = FALSE,
                          treecarbon = FALSE, HDmodel = FALSE,
                          stemdata, allompardata, volumepardata) {
  Treeno <- output <- NULL

  # list estimation data
  output <- list()

  # Extract data of a given tree by selecting the appropriate Treeno
  stemdata1 <- subset(stemdata, Treeno == xtree)

  # Check if there is enough data to compute the stem analysis
  treeDBH <- stemdata1[which(stemdata1[,4] == 1.3), 6]
  if (treeDBH < 4.1) {
    stop("The tree DBH was smaller than 4.1 cm")
  }


  treeTH <- with(stemdata1, max(stemheight))
  if (treeTH < 1.5) {
    stop("The total tree height was lower than 1.5 m")
  }


  if (xtree != max(stemdata1$Treeno)) {
    stop("The data has not included in the stemdata")
  }


  if (!is.logical(stemgrowth)) {
    stop("stemgrowth argument must be a boolean")
  }


  if (!is.logical(treecarbon)) {
    stop("treecarbon argument must be a boolean")
  }


  if (!is.logical(HDmodel)) {
    stop("HDmodel argument must be a boolean")
  }


   # Calculating tree growth and carbon accumulation with stem analysis data

  # Set the age class according to the stem analysis
  # Obtain the starting column for the maximum age of discs along the stem
  stemj <- grep("Dnobark0", names(stemdata1))
  sstart <- stemj
  ys <- stemdata1[1, stemj]
  while (ys != 0) {
    sstart <- sstart + 1
    ys <- stemdata1[1, sstart]
  }
  Agemax <- with(stemdata1, max(stemage)) # Tree age
  nn <- sstart - stemj
  if (Agemax %% nn == 0) {
    Ageclass <- Agemax %/% nn
  } else {
    Ageclass <- Agemax %/% (nn - 1)
  } # set age class

  # Extract tree basic parameters
  hmax <- with(stemdata1, max(stemheight)) # Tree height
  Dmax <- with(stemdata1, max(Dwithbark)) # Diameter with bark at ground


  # Calculate the number of discs
  stemi <- seq_len(length.out = length(stemdata1$stemheight))



  # Calculate tree height

  # Obtain the starting column for the maximum age of discs along the stem
  # stemj<-grep("Dnobark0",names(stemdata1))
  if (Agemax %% Ageclass == 0) {
    # The ages of cross-sectional discs by which the Agemax is divisible
    stemdj <- seq(0, Agemax, Ageclass)
  } else {
    # The ages of cross-sectional discs by which Agemax is not divisible
    stemdj <- c(seq(0, Agemax, Ageclass), Agemax)
  }


  # Set the initial value for calculating the stem height of each disc
  TreeH <- seq(1, length(stemdj) - 2, 1)

  # Find the position of the disc that does not reach the height of next disc
  x <- ceiling(max(stemdj / Ageclass))
  for (k in seq(1, x - 1, 1)) {
    ind <- which(stemdata1[, stemj + k] == 0.00)[1]

    # Calculate tree heights at a given age using Ratio method
    Hitemp <- stemdata1$stemheight[ind - 1] +
      ((stemdata1[ind - 1, stemj + k]) /
         (stemdata1[ind - 1, stemj + k - 1])) *
      (stemdata1$stemheight[ind] - stemdata1$stemheight[ind - 1])
    TreeH[k] <- Hitemp
  }
  TreeH <- as.numeric(TreeH)
  # Calculate the diameter with bark of a given tree age

  # Construct the linear relationship between Dwithbark and Dnobark0 of discs
  DBmodel <- lm(Dwithbark ~ Dnobark0, stemdata1)
  #summary(DBmodel) # Check the model fit effect
  stemdata2 <- stemdata1

  # Predict diameter with bark of a given tree age via the above model
  for (i in seq_len(length.out = nrow(stemdata1))) {
    for (j in stemj:(sstart - 1)) {
      if (stemdata2[i, j] != 0) {
        colnames(stemdata2)[j] <- "Dnobark0"
        stemdata2 <- as.data.frame(stemdata2)
        stemdata1[i, j] <- predict(DBmodel, list(Dnobark0 = stemdata2[i, j]))
      } else {
        stemdata1[i, j] <- stemdata1[i, j]
      }
    }
  }



  # Calculate the stem volume with bark of a given tree age

  # Set the initial value for the stem volume of each individual log
  Vage <- seq(1, length(stemdj) - 2, 1)

  # Set the initial value for the stem volume of the “hidden tip” of the disc
  VageH <- seq(1, length(stemdj) - 2, 1)

  # Find the position of the disc that does not reach to the height of next disc
  x <- ceiling(max(stemdj / Ageclass))
  for (j in seq(1, x - 1, 1)) {
    ind <- which(stemdata1[, stemj + j] == 0.00)[1]

    # Set the initial value for calculating the stem volume of each disc
    VageTemp <- 0
    for (i in seq(1, ind - 2)) {
      # Calculate the stem volume without “hidden tip”
      if (i > 0) {
        Vagej <- VageTemp + 3.14 * ((stemdata1[i, (stemj + j)] / 200)^2 +
                                      ((stemdata1[i, (stemj + j)] / 200)+
                                         (stemdata1[i + 1, (stemj + j)] / 200))^2+
                                      (stemdata1[i + 1, (stemj + j)] / 200)^2) *
          (stemdata1$stemheight[i + 1] - stemdata1$stemheight[i]) / 6
        VageTemp <- Vagej

        # Calculate the stem volume of “hidden tip”
        VHTemp <- 3.14 * ((stemdata1[i, (stemj + j)] / 200)^2) *
          (stemdata1$stemheight[i + 1] - stemdata1$stemheight[i]) / 3
        VageH[j] <- VHTemp
      } else {
        Vagej <- 0
      }
    }
    Vage[j] <- Vagej
  }


  # Calculate the stem volume with bark of a given tree age
  Vage <- as.numeric(Vage)
  VageH <- as.numeric(VageH)
  Vlage <- Vage + VageH


  # Set the initial value for calculating the total stem volume with bark
  VbarkTemp <- 0

  # Calculate the total stem volume without “hidden tip”
  for (i in seq(1, length(stemi) - 2)) {
    Vbark <- VbarkTemp + 3.14 * ((stemdata1$Dnobark0[i] / 200)^2 +
                                   (stemdata1$Dnobark0[i + 1] / 200)^2) *
      (stemdata1$stemheight[i + 1] - stemdata1$stemheight[i]) / 2
    VbarkTemp <- Vbark
  }

  # Calculate the total stem volume with “hidden tip”
  Vbark <- VbarkTemp +
    3.14 * ((stemdata1$Dnobark0[length(stemi) - 1] / 200)^2) *
    (stemdata1$stemheight[length(stemi)] -
       stemdata1$stemheight[length(stemi) - 1]) / 3


  # Extract DBH data for each age class
  DBH1 <- stemdata1[which(stemdata1$stemheight == 1.3), stemj:(stemj + x - 1)]
  DBH2 <- cbind(0, rev(DBH1))
  DBH3 <- as.data.frame(DBH2)
  DBH <- t(DBH3)
  colnames(DBH) <- c("DBHt")
  DBH

  # Construct tree height data for each age class
  Height <- c(0, rev(TreeH), hmax)

  # Construct stem volume with bark data for each age class
  Volume <- c(0, rev(Vlage), Vbark)

  # Construct data for each age class
  Mydata <- data.frame(stemdj, DBH, Height, Volume)


  # Calculate the mean annual increment of stem growth
  AvincreD <- with(Mydata, DBHt / stemdj)
  AvincreH <- with(Mydata, Height / stemdj)
  AvincreV <- with(Mydata, Volume / stemdj)
  AvincreD[is.na(AvincreD)] <- 0
  AvincreH[is.na(AvincreH)] <- 0
  AvincreV[is.na(AvincreV)] <- 0

  # Calculate the current annual increment of stem growth
  AnincreD1 <- seq(1, x, 1)
  AnincreH1 <- seq(1, x, 1)
  AnincreV1 <- seq(1, x, 1)
  for (i in seq(1, x, 1)) {
    AnincreD2 <- with(Mydata, (DBHt[i + 1] - DBHt[i]) /
                        (stemdj[i + 1] - stemdj[i]))
    AnincreD1[i] <- AnincreD2
    AnincreH2 <- with(Mydata, (Height[i + 1] - Height[i]) /
                        (stemdj[i + 1] - stemdj[i]))
    AnincreH1[i] <- AnincreH2
    AnincreV2 <- with(Mydata, (Volume[i + 1] - Volume[i]) /
                        (stemdj[i + 1] - stemdj[i]))
    AnincreV1[i] <- AnincreV2
  }
  AnincreD <- c(0, AnincreD1)
  AnincreH <- c(0, AnincreH1)
  AnincreV <- c(0, AnincreV1)

  # Construct data frame containing DBH, height, and volume growth data
  StemGrowth <- data.frame(Mydata, AnincreD, AvincreD,
                           AnincreH, AvincreH, AnincreV, AvincreV)
  output$StemGrowth <- round(StemGrowth,3)



  # Estimation of tree biomass and carbon storage
  # Check if there a allompardata exists
  if(!missing(allompardata)){

    # Calculate tree aboveground biomass using allometric models
    abovegrounddata <- subset(allompardata, allompardata$tissues=="aboveground")
    abovegrounddata0 <- subset(abovegrounddata, abovegrounddata$DBH=="<5")
    abovegrounddata5 <- subset(abovegrounddata, abovegrounddata$DBH==">5")

    abovegroundB <- ifelse(DBH<5,
                           exp(log(abovegrounddata0$a)+
                                 abovegrounddata0$b*log(DBH)+
                                 abovegrounddata0$c*log(Height)),
                           exp(log(abovegrounddata5$a)+
                                 abovegrounddata5$b*log(DBH)+
                                 abovegrounddata5$c*log(Height)))
    abovegroundB <- as.numeric(abovegroundB)
    abovegroundB[is.na(abovegroundB)]=0

    # Calculate tree belowground biomass using allometric models
    belowgrounddata <- subset(allompardata, allompardata$tissues=="belowground")
    belowgrounddata0 <- subset(belowgrounddata, belowgrounddata$DBH=="<5")
    belowgrounddata5 <- subset(belowgrounddata, belowgrounddata$DBH==">5")

    belowgroundB <- ifelse(DBH<5,
                           exp(log(belowgrounddata0$a)+
                                 belowgrounddata0$b*log(DBH)+
                                 belowgrounddata0$c*log(Height)),
                           exp(log(belowgrounddata5$a)+
                                 belowgrounddata5$b*log(DBH)+
                                 belowgrounddata5$c*log(Height)))
    belowgroundB <- as.numeric(belowgroundB)
    belowgroundB[is.na(belowgroundB)]=0

    totalB <- abovegroundB+belowgroundB
    totalB <- as.numeric(totalB)

    # Calculate tree carbon storage
    abovegroundC <-  abovegroundB*abovegrounddata5$Cconcentration
    abovegroundC <- as.numeric(abovegroundC)

    belowgroundC <-  belowgroundB*belowgrounddata5$Cconcentration
    belowgroundC <- as.numeric(belowgroundC)

    totalC <- abovegroundC + belowgroundC
    totalC <- as.numeric(totalC)

    # Construct data frame that contains tree biomass and carbon storage
    treeage <- as.numeric(stemdj)
    allomCarbon <- cbind(treeage, abovegroundB, belowgroundB, totalB,
                         abovegroundC, belowgroundC, totalC)
    allomCarbon <- as.data.frame(allomCarbon)
    output$allomCarbon <- round(allomCarbon,3)
  }

  # Check if there a volumepardata exists
  if(!missing(volumepardata)){

    # Calculate BCF
    BCFdata <- subset(volumepardata, volumepardata$factors=="BCF")
    BCFdata0 <- subset(BCFdata, BCFdata$DBH=="<5")
    BCFdata5 <- subset(BCFdata, BCFdata$DBH==">5")

    BCF <- ifelse(DBH<5,
                  exp(log(BCFdata0$a)+BCFdata0$b*log(DBH)+BCFdata0$c*log(Height)),
                  exp(log(BCFdata5$a)+BCFdata5$b*log(DBH)+BCFdata5$c*log(Height)))
    BCF <- as.numeric(BCF)
    BCF[is.na(BCF)]=0


    # Calculate RSR
    Rdata <- subset(volumepardata, volumepardata$factors=="RSR")
    Rdata0 <- subset(Rdata, Rdata$DBH=="<5")
    Rdata5 <- subset(Rdata, Rdata$DBH==">5")

    RSR <- ifelse(DBH<5,
                  exp(log(Rdata0$a)+Rdata0$b*log(DBH)+Rdata0$c*log(Height)),
                  exp(log(Rdata5$a)+Rdata5$b*log(DBH)+Rdata5$c*log(Height)))
    RSR <- as.numeric(RSR)
    RSR[is.na(RSR)]=0

    # Calculate tree biomass using volume model
    abovegroundB <- ifelse(Volume!=0,Volume*BCF*1000,Volume)
    abovegroundB <- as.numeric(abovegroundB)

    belowgroundB <- ifelse(Volume!=0,Volume*BCF*1000*RSR,Volume)
    belowgroundB <- as.numeric(belowgroundB)

    totalB <- ifelse(Volume!=0,Volume*BCF*(1+RSR)*1000,Volume)
    totalB <- as.numeric(totalB)

    # Calculate tree carbon storage
    abovegroundC <- abovegroundB*BCFdata5$Cconcentration
    abovegroundC <- as.numeric(abovegroundC)

    belowgroundC <- belowgroundB*BCFdata5$Cconcentration
    belowgroundC <- as.numeric(belowgroundC)

    totalC <- abovegroundC+belowgroundC
    totalC <- as.numeric(totalC)


    # Construct data frame that contains tree biomass and tree carbon storage
    treeage <- as.numeric(stemdj)

    volumeCarbon <- cbind(treeage, BCF, RSR, abovegroundB, belowgroundB, totalB,
                          abovegroundC, belowgroundC, totalC)
    volumeCarbon <- as.data.frame(volumeCarbon)
    output$volumeCarbon <- round(volumeCarbon,3)
  }




  #######################################################
  # define drawStemgrowth() function which show the stem growth trend
  drawStemgrowth <- function() {

    # Open a new graph frame which names as “Stem growth trend”
    dev.new(title = "Stem growth trend", width = 4500,
            height = 4500, noRStudioGD = TRUE)


    # Define the parameters of the graph
    mat <- matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), 3, 3, byrow = FALSE)

    # mat
    layout(mat)
    layout.show(7)
    par(oma = c(0, 0, 0, 0), mar = c(5, 5, 1, 1))

    # Draw the base coordinates
    plot(c(0, 0), c(0, hmax), type = "l", lty = "dotted", col = "red",
         lwd = 2, xlim = c(-ceiling(Dmax / 20), ceiling(Dmax / 20)),
         ylim = c(0, ceiling(hmax)), ann = FALSE, axes = FALSE)
    legend("topleft", "(a)", cex = 1.1, box.lty = 0, bg = NULL)
    axis(1, at = seq(-ceiling(Dmax / 20), ceiling(Dmax / 20), 0.2),
         lwd = 1.5, font = 1, tck = 0.02, cex = 0.4)

    # Calculate the number of discs
    stemi <- seq_len(length.out = length(stemdata1$stemheight))

    # Draw the position lines for discs at different stem heights
    for (i in stemi) {
      lines(c(-stemdata1$Dwithbark[i] / 20, stemdata1$Dwithbark[i] / 20),
            c(stemdata1$stemheight[i], stemdata1$stemheight[i]),
            lty = "dotted", col = "red", lwd = 2)
    }
    # Draw the lateral surface with and without bark
    for (i in stemi - 1) {
      lines(c(stemdata1$Dwithbark[i] / 20, stemdata1$Dwithbark[i + 1] / 20),
            c(stemdata1$stemheight[i], stemdata1$stemheight[i + 1]),
            lty = "solid", col = "black", lwd = 2)
      lines(c(-stemdata1$Dwithbark[i] / 20, -stemdata1$Dwithbark[i + 1] / 20),
            c(stemdata1$stemheight[i], stemdata1$stemheight[i + 1]),
            lty = "solid", col = "black", lwd = 2)
      lines(c(stemdata1$ Dnobark0[i] / 20, stemdata1$ Dnobark0[i + 1] / 20),
            c(stemdata1$stemheight[i], stemdata1$stemheight[i + 1]),
            lty = "solid", col = "black", lwd = 1)
      lines(c(-stemdata1$ Dnobark0[i] / 20, -stemdata1$ Dnobark0[i + 1] / 20),
            c(stemdata1$stemheight[i], stemdata1$stemheight[i + 1]),
            lty = "solid", col = "black", lwd = 1)
    }

    # Find the disc that does not reach the height of next cross-sectional disc
    # Set the initial value for calculating the stem height of each disc
    TreeH <- seq(1, length(stemdj) - 2, 1)
    x <- ceiling(max(stemdj / Ageclass))

    for (k in seq(1, x - 1, 1)) {
      ind <- which(stemdata1[, stemj + k] == 0.00)[1]

      # Calculate tree heights at a given age using Ratio method
      Hitemp <- stemdata1$stemheight[ind - 1] +
        ((stemdata1[ind - 1, stemj + k]) /
           (stemdata1[ind - 1, stemj + k - 1])) *
        (stemdata1$stemheight[ind] - stemdata1$stemheight[ind - 1])
      TreeH[k] <- Hitemp

      # Draw the lateral surface of the individual logs
      if(ind-2==0){
      lines(c(stemdata1[1, (stemj + k)] / 20,
              stemdata1[2, (stemj + k)] / 20),
            c(0, stemdata1$stemheight[ind]),
            lty = "solid", col = "black", lwd = 1)
      lines(c(-stemdata1[1, (stemj + k)] / 20,
              -stemdata1[2, (stemj + k)] / 20),
            c(0, stemdata1$stemheight[ind]),
            lty = "solid", col = "black", lwd = 1)
      }else if(ind-2==1){
      lines(c(stemdata1[ind-2, (stemj + k)] / 20,
              stemdata1[ind-2+1, (stemj + k)] / 20),
            c(stemdata1$stemheight[ind-2], stemdata1$stemheight[ind-2+1]),
            lty = "solid", col = "black", lwd = 1)
      lines(c(-stemdata1[ind-2, (stemj + k)] / 20,
              -stemdata1[ind-2+1, (stemj + k)] / 20),
            c(stemdata1$stemheight[ind-2], stemdata1$stemheight[ind-2+1]),
            lty = "solid", col = "black", lwd = 1)


    }else if(ind-2>1){
            for (i in seq(1, ind - 2)) {
                  lines(c(stemdata1[i, (stemj + k)] / 20,
                stemdata1[i+1, (stemj + k)] / 20),
              c(stemdata1$stemheight[i], stemdata1$stemheight[i+1]),
              lty = "solid", col = "black", lwd = 1)
        lines(c(-stemdata1[i, (stemj + k)] / 20,
                -stemdata1[i+1, (stemj + k)] / 20),
              c(stemdata1$stemheight[i], stemdata1$stemheight[i+1]),
              lty = "solid", col = "black", lwd = 1)
            }
    }


      # Draw the lateral surface of the “hidden tip” of the disc
      lines(c(stemdata1[ind - 1, (stemj + k)] / 20, 0),
            c(stemdata1$stemheight[ind - 1], Hitemp),
            lty = "solid", col = "black", lwd = 1)
      lines(c(-stemdata1[ind - 1, (stemj + k)] / 20, 0),
            c(stemdata1$stemheight[ind - 1], Hitemp),
            lty = "solid", col = "black", lwd = 1)
      }



    # Draw graphics for the cumulative growth
    with(output$StemGrowth, plot(stemdj, DBHt, type = "b", pch = 16,
                                 col = "forestgreen", lwd = 2, cex = 1.5,
                                 xlim = c(0, max(stemdj) + Ageclass),
                                 ylim = c(0, 1.2 * max(DBHt)), xlab = "Age (years)",
                                 ylab = "DBH (cm)", las = 1, cex.lab = 1.1))
    legend("topleft", "(b)", cex = 1.1, box.lty = 0, bg = NULL)

    with(output$StemGrowth, plot(stemdj, Height, type = "b", pch = 16,
                                 col = "forestgreen", lwd = 2, cex = 1.5,
                                 xlim = c(0, max(stemdj) + Ageclass),
                                 ylim = c(0, 1.2 * max(Height)), xlab = "Age (years)",
                                 ylab = "Height (m)", las = 1, cex.lab = 1.1))
    legend("topleft", "(c)", cex = 1.1, box.lty = 0, bg = NULL)

    with(output$StemGrowth, plot(stemdj, Volume, type = "b", pch = 16,
                                 col = "forestgreen", lwd = 2, cex = 1.5,
                                 xlim = c(0, max(stemdj) + Ageclass),
                                 ylim = c(0, 1.2 * max(Volume)), xlab = "Age (years)",
                                 ylab = expression(paste("Volume (m"^"3", ")")),
                                 las = 1, cex.lab = 1.1))
    legend("topleft", "(d)", cex = 1.1, box.lty = 0, bg = NULL)

    # Draw graphics for the mean annual increment and current annual increment
    with(output$StemGrowth, plot(stemdj, AnincreD, type = "l", lty = "dashed", pch = 16,
                                 col = "blue", lwd = 2, cex = 1.2,
                                 xlim = c(0, max(stemdj) + Ageclass),
                                 ylim = c(0, 1.2 * (max(max(AnincreD), max(AvincreD)))),
                                 xlab = "Age (years)",
                                 ylab = expression(paste
                                                   ("DBH increment (cm", " year"^"-1", ")")),
                                 las = 1, cex.lab = 1.1))
    with(output$StemGrowth, lines(stemdj, AvincreD, lty = "dotted", pch = 16,
                                  col = "red", lwd = 2, cex = 1.2))
    legend("topleft", "(e)", cex = 1.1, box.lty = 0, bg = NULL)
    legend("topright",
           col = c("blue", "red"), cex = 1, box.lty = 0, lwd = 2,
           lty = c("dashed", "dotted"), bg = NULL,
           c("Current annual increment", "Mean annual increment")
    )

    with(output$StemGrowth, plot(stemdj, AnincreH, type = "l", lty = "dashed", pch = 1,
                                 col = "blue", lwd = 2, cex = 1.2,
                                 xlim = c(0, max(stemdj) + Ageclass),
                                 ylim = c(0, 1.2 * (max(max(AnincreH), max(AvincreH)))),
                                 xlab = "Age (years)",
                                 ylab = expression(paste
                                                   ("Height increment (m", " year"^"-1", ")")),
                                 las = 1, cex.lab = 1.1))
    with(output$StemGrowth, lines(stemdj, AvincreH, lty = "dotted", pch = 16,
                                  col = "red", lwd = 2, cex = 1.2))
    legend("topleft", "(f)", cex = 1.1, box.lty = 0, bg = NULL)

    with(output$StemGrowth, plot(stemdj, AnincreV, type = "l", lty = "dashed", pch = 1,
                                 col = "blue", lwd = 2, cex = 1.2,
                                 xlim = c(0, max(stemdj) + Ageclass),
                                 ylim = c(0, 1.2 * (max(max(AnincreV), max(AvincreV)))),
                                 xlab = "Age (years)",
                                 ylab = expression(paste
                                                   ("Volume increment (m"^"3", " year"^"-1", ")"),
                                                   las = 1, cex.lab = 1.1)))
    with(output$StemGrowth, lines(stemdj, AvincreV, lty = "dotted", pch = 16,
                                  col = "red", lwd = 2, cex = 1.2))
    legend("topleft", "(g)", cex = 1.1, box.lty = 0, bg = NULL)
  }

  # Whether to call the stemgrowth graph
  if (stemgrowth == TRUE) {
    drawStemgrowth()
  }



  #####################################################################
  #define calCarbon() function
  calCarbon <- function(allompardata, volumepardata){

    # Check if there a allompardata exists
    if(!missing(allompardata)){

      # Combined tree growth trend data and carbon storage data
      Mydata2 <- data.frame(output$StemGrowth, output$allomCarbon)

      # Open a new graph frame
      dev.new(title =
                "tree biomass and carbon storage estimated by allometric models",
              width = 4500, height = 2250, noRStudioGD = TRUE)

      # Define the parameters of graph
      mat <- matrix(c(1,2),1,2,byrow=FALSE)
      #mat
      layout(mat)
      layout.show(2)
      par(oma=c(0,0,0,0),mar=c(5,5,1,1))

      # Draw a graph for tree biomass and carbon storage across tree age
      # The changes in tree biomass across tree age
      with(Mydata2, plot(treeage, totalB, type="b", lty=1, lwd=2, pch = 19,
                         col = "blue", cex.axis=1.3, cex.lab=1.5,cex = 2,
                         xlim = c(0, max(treeage)+Ageclass),
                         ylim = c(0, 1.2*(max(totalB))),
                         xlab = "Age (years)",
                         ylab = "Total tree biomass (kg)",las=1))
      legend("topleft", "(a)", cex = 1.5, box.lty = 0, bg=NULL)
      legend("topright", cex = 1.5, box.lty = 0, bg=NULL,
             c("Estimated by allometric models"))

      # The changes in carbon storage across tree age
      with(Mydata2, plot(treeage, totalC, type="b", lty=1, lwd=2,pch = 19,
                         col = "blue", cex.axis=1.3,cex.lab=1.5,cex = 2,
                         xlim = c(0, max(treeage)+Ageclass),
                         ylim = c(0, 1.2*(max(max(totalC)))),
                         xlab = "Age (years)",
                         ylab = "Total tree C storage (kg)", las=1))
      legend("topleft", "(b)", cex = 1.5, box.lty = 0, bg=NULL)
      legend("topright", cex = 1.5, box.lty = 0, bg=NULL,
             c("Estimated by allometric models"))
    }


    # Check if there a volumepardata exists
    if(!missing(volumepardata)){

      # Combined tree growth trend data and carbon storage data
      Mydata3 <- data.frame(output$StemGrowth, output$volumeCarbon)

      # Open a new graph frame
      dev.new(title =
                "tree biomass and carbon storage estimated by volume model",
              width = 4500, height = 2250, noRStudioGD = TRUE)

      # Define the parameters of graph
      mat <- matrix(c(1,2),1,2,byrow=FALSE)
      # mat
      layout(mat)
      layout.show(2)
      par(oma=c(0,0,0,0),mar=c(5,5,1,1))

      # Draw a graph for tree biomass and carbon storage across tree age
      # The changes in tree biomass across tree age
      with(Mydata3, plot(treeage, totalB, type="b", lty=1, lwd=2, pch = 19,
                         col = "blue", cex.axis=1.3, cex.lab=1.5,cex = 2,
                         xlim = c(0, max(treeage)+Ageclass),
                         ylim = c(0, 1.2*(max(totalB))),
                         xlab = "Age (years)",
                         ylab = "Total tree biomass (kg)",las=1))
      legend("topleft", "(a)", cex = 1.5, box.lty = 0, bg=NULL)
      legend("topright", cex = 1.5, box.lty = 0, bg=NULL,
             c("Estimated by volume model"))

      # The changes in carbon storage across tree age
      with(Mydata3, plot(treeage, totalC, type="b", lty=1, lwd=2,pch = 19,
                         col = "blue", cex.axis=1.3,cex.lab=1.5,cex = 2,
                         xlim = c(0, max(treeage)+Ageclass),
                         ylim = c(0, 1.2*(max(max(totalC)))),
                         xlab = "Age (years)",
                         ylab = "Total tree C storage (kg)", las=1))
      legend("topleft", "(b)", cex = 1.5, box.lty = 0, bg=NULL)
      legend("topright", cex = 1.5, box.lty = 0, bg=NULL,
             c("Estimated by volume model"))
    }
  }

  # Whether to call the treecarbon graph
  if (treecarbon == TRUE) {
    calCarbon(allompardata, volumepardata)
  }




  ###############################################################
  # Define allomHD() function
  allomHD <- function() {
    # Open a new graph frame
    dev.new(title = "height-diameter relationships",
            width = 5000, height = 5500, noRStudioGD = TRUE)

    # Define the parameters of graph
    mat <- matrix(c(1, 1, 2), 3, 1, byrow = FALSE)
    # mat
    layout(mat)
    layout.show(2)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5, 5, 1, 1))

    # Extract sub-datasets with tree height greater than 1.3 m
    Mydata <- subset(output$StemGrowth, Height >= 1.3)

    # Draw a scatter plot
    plot(Mydata$DBHt, Mydata$Height,
         pch = 21, bg = "purple",
         cex.axis = 1.5, cex.lab = 1.8, cex = 3,
         xlim = c(0, 1.1 * max(Mydata$DBHt)),
         ylim = c(0, 1.1 * max(Mydata$Height)),
         xlab = "Tree DBH (cm)", ylab = "Tree height (m)", las = 1
    )

    # Develop Chapman-Richards model
    theta_richards <- lmfor::startHDrichards(d = Mydata$DBHt, h = Mydata$Height)
    nlc <- nls.control(maxiter = 1000)
    HDrichards <- nls(Height ~ 1.3 + a * (1 - exp(-b * DBHt))^c,
                      control = nlc,
                      start = list(a = theta_richards[1],
                                   b = theta_richards[2],
                                   c = theta_richards[3]), data = Mydata)
    summary(HDrichards)
    xv1_richards <- seq(0.3 * min(Mydata$DBHt), 1.1 * max(Mydata$DBHt), 0.0001)
    yv1_richards <- predict(HDrichards, list(DBHt = xv1_richards))
    lines(xv1_richards, yv1_richards, col = "orangered", lwd = 3)
    SSre_richards <- sum(residuals(HDrichards)^2)
    SStot_richards <- sum((Mydata$Height - mean(Mydata$Height))^2)
    R2_richards <- 1 - SSre_richards / SStot_richards
    R2_richards
    N <- nrow(Mydata)
    RMSE_richards <- sqrt(sum(residuals(HDrichards)^2) / (N - 2))
    RMSE_richards

    # Develop Logistic model
    theta_logistic <- lmfor::startHDlogistic(d = Mydata$DBHt, h = Mydata$Height)
    nlc <- nls.control(maxiter = 1000)
    HDlogistic <- nls(Height ~ 1.3 + a / (1 + b * exp(-c * DBHt)),
                      control = nlc,
                      start = list(a = theta_logistic[1],
                                   b = theta_logistic[2],
                                   c = theta_logistic[3]), data = Mydata)
    summary(HDlogistic)
    xv1_logistic <- seq(0.3 * min(Mydata$DBHt), 1.1 * max(Mydata$DBHt), 0.0001)
    yv1_logistic <- predict(HDlogistic, list(DBHt = xv1_logistic))
    lines(xv1_logistic, yv1_logistic, col = "blue", lwd = 3)
    SSre_logistic <- sum(residuals(HDlogistic)^2)
    SStot_logistic <- sum((Mydata$Height - mean(Mydata$Height))^2)
    R2_logistic <- 1 - SSre_logistic / SStot_logistic
    R2_logistic
    N <- nrow(Mydata)
    RMSE_logistic <- sqrt(sum(residuals(HDlogistic)^2) / (N - 2))
    RMSE_logistic

    # Develop Weibull model
    theta_weibull <- lmfor::startHDweibull(d = Mydata$DBHt, h = Mydata$Height)
    nlc <- nls.control(maxiter = 1000)
    HDweibull <- nls(Height ~ 1.3 + a * (1 - exp(-b * DBHt^c)), control = nlc,
                     start = list(a = theta_weibull[1],
                                  b = theta_weibull[2],
                                  c = theta_weibull[3]), data = Mydata)
    summary(HDweibull)
    xv1_weibull <- seq(0.3 * min(Mydata$DBHt), 1.1 * max(Mydata$DBHt), 0.0001)
    yv1_weibull <- predict(HDweibull, list(DBHt = xv1_weibull))
    lines(xv1_weibull, yv1_weibull, col = "forestgreen", lwd = 3)
    SSre_weibull <- sum(residuals(HDweibull)^2)
    SStot_weibull <- sum((Mydata$Height - mean(Mydata$Height))^2)
    R2_weibull <- 1 - SSre_weibull / SStot_weibull
    R2_weibull
    N <- nrow(Mydata)
    RMSE_weibull <- sqrt(sum(residuals(HDweibull)^2) / (N - 2))
    RMSE_weibull

    # Develop gomperz model
    theta_gomperz <- lmfor::startHDgomperz(d = Mydata$DBHt, h = Mydata$Height)
    nlc <- nls.control(maxiter = 1000)
    HDgomperz <- nls(Height ~ 1.3 + a * exp(-b * exp(-c * DBHt)),
                     control = nlc,
                     start = list(a = theta_gomperz[1],
                                  b = theta_gomperz[2],
                                  c = theta_gomperz[3]), data = Mydata)
    summary(HDgomperz)
    xv1_gomperz <- seq(0.3 * min(Mydata$DBHt), 1.1 * max(Mydata$DBHt), 0.0001)
    yv1_gomperz <- predict(HDgomperz, list(DBHt = xv1_gomperz))
    lines(xv1_gomperz, yv1_gomperz, col = "purple", lwd = 3)
    SSre_gomperz <- sum(residuals(HDgomperz)^2)
    SStot_gomperz <- sum((Mydata$Height - mean(Mydata$Height))^2)
    R2_gomperz <- 1 - SSre_gomperz / SStot_gomperz
    R2_gomperz
    N <- nrow(Mydata)
    RMSE_gomperz <- sqrt(sum(residuals(HDgomperz)^2) / (N - 2))
    RMSE_gomperz

    legend("topleft",
           cex = 1.5, lwd = 2.5, box.lty = 0, bg = NULL,
           text.col = c("orangered", "blue", "forestgreen", "purple"),
           legend = c(
             expression(paste
                        ("HDrichards  H = 1.3 + a(1 - e"^"-bDBH","",")"^"c")),
             expression(paste
                        ("HDlogistic  H = 1.3 + a/(1 + be"^"-cDBH","",")")),
             expression(paste
                        ("HDweibull  H = 1.3 + a(1 - e"^"-bDBH"^"c","",")")),
             expression(paste
                        ("HDgomperz  H = 1.3 + ae"^("-be"^"-cDBH")))),
           lty = 1, col = c("orangered", "blue", "forestgreen", "purple"),
           title = "Models"
    )

    # Show the fitted statistics
    # Draw a scatter plot
    par(new=T)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(0, 0, 0, 0))
    plot(Mydata$DBHt, Mydata$Height,
         type = "n",xaxt = "n", yaxt ="n", bty = "n", ann = FALSE,
         xlim = c(0, 1.1 * max(Mydata$DBHt)),
         ylim = c(0, 1.1 * max(Mydata$Height))
    )

    # Filter the model providing the best fit by anova
    min_RSS <- min(anova(HDrichards, HDlogistic, HDgomperz, HDweibull)[2])

    HDrichards_RSS <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[1, 2]
    HDlogistic_RSS <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[2, 2]
    HDweibull_RSS <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[3, 2]
    HDgomperz_RSS <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[4, 2]

    # Filter the model providing the best fit by AIC
    HDrichards_AIC <- AIC(HDrichards, HDlogistic, HDweibull, HDgomperz)[1,2]
    HDlogistic_AIC <- AIC(HDrichards, HDlogistic, HDweibull, HDgomperz)[2,2]
    HDweibull_AIC <- AIC(HDrichards, HDlogistic, HDweibull, HDgomperz)[3,2]
    HDgomperz_AIC <- AIC(HDrichards, HDlogistic, HDweibull, HDgomperz)[4,2]

    legend("left",
           cex = 1.5, bty = "n",
           legend = c(
             "Equation", "HDrichards", "HDlogistic", "HDweibull", "HDgomperz",
             "a",round(summary(HDrichards)$parameters[1], 3),
             round(summary(HDlogistic)$parameters[1], 3),
             round(summary(HDweibull)$parameters[1], 3),
             round(summary(HDgomperz)$parameters[1], 3),
             "b", round(summary(HDrichards)$parameters[2], 3),
             round(summary(HDlogistic)$parameters[2], 3),
             round(summary(HDweibull)$parameters[2], 3),
             round(summary(HDgomperz)$parameters[2], 3),
             "c",round(summary(HDrichards)$parameters[3], 3),
             round(summary(HDlogistic)$parameters[3], 3),
             round(summary(HDweibull)$parameters[3], 3),
             round(summary(HDgomperz)$parameters[3], 3),
             "R2", round(R2_richards, 3), round(R2_logistic, 3),
             round(R2_weibull, 3), round(R2_gomperz, 3),
             "RSS",round(HDrichards_RSS, 3), round(HDlogistic_RSS, 3),
             round(HDweibull_RSS, 3), round(HDgomperz_RSS, 3),
             "AIC",round(HDrichards_AIC, 3), round(HDlogistic_AIC, 3),
             round(HDweibull_AIC, 3), round(HDgomperz_AIC, 3),
             "logLik",round(logLik(HDrichards), 3), round(logLik(HDlogistic), 3),
             round(logLik(HDweibull), 3), round(logLik(HDgomperz), 3)),

           title = "Fitted statistics of HDmodel", title.font = 2,
           ncol = 8, text.width = c(1.6,0.75,0.75,0.75,0.75,0.75,0.75,0.75))
  }

  # Whether to call the HDmodel graph
  if (HDmodel == TRUE) {
    allomHD()
  }

  return(output)
}


