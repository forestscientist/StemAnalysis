#' Reconstructing tree growth and carbon accumulation with stem analysis data
#'
#' @param xtree Xtree is the tree number (Treeno), which isused to choose the
#'     target tree to be analyzed
#' @param stemgrowth If stemgrowth is 'TRUE', stem growth profile and growth
#'     trends in terms of diameter at breast height (DBH), tree height, and
#'     stem volume will be showed in a graph
#' @param treecarbon If treecarbon is 'TRUE', total tree biomass and carbon
#'     storage will be estimated by allometric models and volume model. In
#'     addition, although treecarbon is 'TRUE', the estimation of tree biomass
#'     and carbon storage by allometric models will skip if data 'parameterdata'
#'     is missing, and the same is true for the estimation by volume model if
#'     data 'BEFdata' is missing
#' @param HDmodel If HDmodel is 'TRUE', height-diameter relationship will be
#'     fitted with nonlinear models and showed the fitted results in a graph
#' @param stemdata Stemdata is the stem analysis data that has been inputted
#' @param parameterdata Parameterdata is the parameter data of allometric
#'     models that can be optionally inputted
#' @param BEFdata BEFdata is the biomass estimation factor data of volume model
#'     that can be optionally inputted by users
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' library(StemAnalysis)
#'
#' # Load the data sets
#' data(stemdata)
#' data(parameterdata)
#'
#' # To calculating tree growth and carbon accumulation with input data sets
#' stemanalysism(xtree = 8, stemdata = stemdata)
#'
#' # If the graph of stem growth profile and growth trends is needed
#' stemanalysism(xtree = 8, stemgrowth = TRUE, treecarbon = TRUE,
#'     stemdata = stemdata, parameterdata = parameterdata)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics axis layout layout.show legend lines mtext par
#' @importFrom stats anova lm nls nls.control predict residuals
#' @importFrom utils write.table



stemanalysism <- function(xtree, stemgrowth = FALSE,
                          treecarbon = FALSE, HDmodel = FALSE, stemdata,
                          parameterdata, BEFdata) {
  Treeno <- NULL

  # Extract data of a given tree by selecting the appropriate Treeno
  stemdata1 <- subset(stemdata, Treeno == xtree)



  # Check if there is enough data to compute the stem analysis
  treeDBH <- max(stemdata1[which(stemdata1$stemheight == 1.3),
                           stemdata1$Dwithbark])
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



  # Calculate the diameter with bark of a given tree age

  # Construct the linear relationship between Dwithbark and Dnobark0 of discs
  DBmodel <- lm(Dwithbark ~ Dnobark0, stemdata1)
  summary(DBmodel) # Check the model fit effect

  # Predict diameter with bark of a given tree age via the above model
  for (i in seq_len(length.out = nrow(stemdata1))) {
    for (j in stemj:(sstart - 1)) {
      if (stemdata1[i, j] != 0) {
        stemdata1[i, j] <- predict(DBmodel, list(Dnobark0 = stemdata1[i, j]))
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
                                      (stemdata1[i + 1, (stemj + j)] / 200)^2) *
          (stemdata1$stemheight[i + 1] - stemdata1$stemheight[i]) / 2
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

  # Construct data frame containing DBH, height, and volume and growth data
  StemGrowth <- data.frame(Mydata, AnincreD, AvincreD,
                        AnincreH, AvincreH, AnincreV, AvincreV)
  write.table(StemGrowth, file ="../StemGrowth.csv", sep =",")

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
    legend("topleft", "(a)", cex = 1.1, box.lty = 0)
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
      for (i in seq(1, ind - 2)) {
        lines(c(stemdata1[i, (stemj + k)] / 20,
                stemdata1[i + 1, (stemj + k)] / 20),
              c(stemdata1$stemheight[i], stemdata1$stemheight[i + 1]),
              lty = "solid", col = "black", lwd = 1)
        lines(c(-stemdata1[i, (stemj + k)] / 20,
                -stemdata1[i + 1, (stemj + k)] / 20),
              c(stemdata1$stemheight[i], stemdata1$stemheight[i + 1]),
              lty = "solid", col = "black", lwd = 1)
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
    with(StemGrowth, plot(stemdj, DBHt, type = "b", pch = 16,
                       col = "forestgreen", lwd = 2, cex = 1.5,
                       xlim = c(0, max(stemdj) + Ageclass),
                       ylim = c(0, 1.2 * max(DBHt)), xlab = "Age (years)",
                       ylab = "DBH (cm)", las = 1, cex.lab = 1.1))
    legend("topleft", "(b)", cex = 1.1, box.lty = 0)

    with(StemGrowth, plot(stemdj, Height, type = "b", pch = 16,
                       col = "forestgreen", lwd = 2, cex = 1.5,
                       xlim = c(0, max(stemdj) + Ageclass),
                       ylim = c(0, 1.2 * max(Height)), xlab = "Age (years)",
                       ylab = "Height (m)", las = 1, cex.lab = 1.1))
    legend("topleft", "(c)", cex = 1.1, box.lty = 0)

    with(StemGrowth, plot(stemdj, Volume, type = "b", pch = 16,
                       col = "forestgreen", lwd = 2, cex = 1.5,
                       xlim = c(0, max(stemdj) + Ageclass),
                       ylim = c(0, 1.2 * max(Volume)), xlab = "Age (years)",
                       ylab = expression(paste("Volume (m"^"3", ")")),
                       las = 1, cex.lab = 1.1))
    legend("topleft", "(d)", cex = 1.1, box.lty = 0)

    # Draw graphics for the mean annual increment and current annual increment
    with(StemGrowth, plot(stemdj, AnincreD, type = "l", lty = "dashed", pch = 16,
                       col = "blue", lwd = 2, cex = 1.2,
                       xlim = c(0, max(stemdj) + Ageclass),
                       ylim = c(0, 1.2 * (max(max(AnincreD), max(AvincreD)))),
                       xlab = "Age (years)",
                       ylab = expression(paste
                                  ("DBH increment (cm", " year"^"-1", ")")),
                       las = 1, cex.lab = 1.1))
    with(StemGrowth, lines(stemdj, AvincreD, lty = "dotted", pch = 16,
                        col = "red", lwd = 2, cex = 1.2))
    legend("topleft", "(e)", cex = 1.1, box.lty = 0)
    legend("topright",
           col = c("blue", "red"), cex = 1, box.lty = 0, lwd = 2,
           lty = c("dashed", "dotted"), bg = NULL,
           c("Current annual increment", "Mean annual increment")
    )

    with(StemGrowth, plot(stemdj, AnincreH, type = "l", lty = "dashed", pch = 1,
                       col = "blue", lwd = 2, cex = 1.2,
                       xlim = c(0, max(stemdj) + Ageclass),
                       ylim = c(0, 1.2 * (max(max(AnincreH), max(AvincreH)))),
                       xlab = "Age (years)",
                       ylab = expression(paste
                                  ("Height increment (m", " year"^"-1", ")")),
                       las = 1, cex.lab = 1.1))
    with(StemGrowth, lines(stemdj, AvincreH, lty = "dotted", pch = 16,
                        col = "red", lwd = 2, cex = 1.2))
    legend("topleft", "(f)", cex = 1.1, box.lty = 0)

    with(StemGrowth, plot(stemdj, AnincreV, type = "l", lty = "dashed", pch = 1,
                       col = "blue", lwd = 2, cex = 1.2,
                       xlim = c(0, max(stemdj) + Ageclass),
                       ylim = c(0, 1.2 * (max(max(AnincreV), max(AvincreV)))),
                       xlab = "Age (years)",
                       ylab = expression(paste
                                ("Volume increment (m"^"3", " year"^"-1", ")"),
                                         las = 1, cex.lab = 1.1)))
    with(StemGrowth, lines(stemdj, AvincreV, lty = "dotted", pch = 16,
                        col = "red", lwd = 2, cex = 1.2))
    legend("topleft", "(g)", cex = 1.1, box.lty = 0)
  }

  # Whether to call the stemgrowth graph
  if (stemgrowth == TRUE) {
    drawStemgrowth()
  } else {
    print("The user does not call the stemgrowth graph")
  }



  #####################################################################
  #define calCarbon() function
  calCarbon <- function(parameterdata, BEFdata){

    # Check if there a parameterdata exists
    # Estimation of tree biomass and carbon storage from allometric models
    if(!missing(parameterdata)){

      # Calculate tree biomass using allometric models
      stem_pa <- subset(parameterdata,parameterdata$tissues=="stem")
      stem_biomass <- ifelse(DBH!=0,
                             exp(stem_pa$a+stem_pa$b*log(DBH^2*Height)), DBH)
      stem_biomass <- as.numeric(stem_biomass)

      branch_pa <- subset(parameterdata,parameterdata$tissues=="branch")
      branch_biomass <- ifelse(DBH!=0,
                               exp(branch_pa$a+branch_pa$b*log(DBH^2*Height)),
                               DBH)
      branch_biomass <- as.numeric(branch_biomass)

      leaf_pa <- subset(parameterdata,parameterdata$tissues=="leaf")
      leaf_biomass <- ifelse(DBH!=0,
                             exp(leaf_pa$a+leaf_pa$b*log(DBH^2*Height)), DBH)
      leaf_biomass <- as.numeric(leaf_biomass)

      root_pa <- subset(parameterdata,parameterdata$tissues=="root")
      root_biomass <- ifelse(DBH!=0,
                             exp(root_pa$a+root_pa$b*log(DBH^2*Height)), DBH)
      root_biomass <- as.numeric(root_biomass)

      total_pa <- subset(parameterdata,parameterdata$tissues=="total")
      total_biomass <- stem_biomass+branch_biomass+root_biomass+leaf_biomass
      total_biomass <- as.numeric(total_biomass)

      # Calculate tree carbon storage
      stem_C <-  stem_biomass*stem_pa$Cconcentration
      stem_C <- as.numeric(stem_C)

      branch_C <- branch_biomass*branch_pa$Cconcentration
      branch_C <- as.numeric(branch_C)

      leaf_C <- leaf_biomass*leaf_pa$Cconcentration
      leaf_C <- as.numeric(leaf_C)

      root_C <- root_biomass*root_pa$Cconcentration
      root_C <- as.numeric(root_C)

      total_C <- total_biomass*total_pa$Cconcentration
      total_C <- as.numeric(total_C)

      # Construct data frame that contains tree biomass and tree carbon storage
      stemdj <- as.numeric(stemdj)
      allomCarbon <- cbind(stemdj, stem_biomass, branch_biomass, leaf_biomass,
                           root_biomass, total_biomass, stem_C, branch_C,
                           leaf_C, root_C, total_C)
      allomCarbon <- as.data.frame(allomCarbon)

      write.table(allomCarbon, file ="../allomCarbon.csv", sep =",")

      # Combined tree growth trend data and carbon storage data
      Mydata2 <- data.frame(StemGrowth, allomCarbon)

      # Open a new graph frame
      dev.new(title =
                "tree biomass and carbon storage estimated by allometric models",
              width = 4500, height = 2250, noRStudioGD = TRUE)

      # Define the parameters of graph
      mat <- matrix(c(1,2),1,2,byrow=FALSE)
      #mat
      layout(mat)
      layout.show(2)
      par(oma=c(0,0,0,0),mar=c(5,5,5,1))

      # Draw a graph for tree biomass and carbon storage across tree age
      # The changes in tree biomass across tree age
      with(Mydata2, plot(stemdj, total_biomass, type="b", lty=1, lwd=2, pch = 19,
                         col = "blue", cex.axis=1.5, cex.lab=1.8,cex = 2,
                         xlim = c(0, max(stemdj)+Ageclass),
                         ylim = c(0, 1.2*(max(total_biomass))),
                         xlab = "Age (years)",
                         ylab = "Total tree biomass (kg)",las=1))
      legend("topleft", "(a)", cex = 1.8, box.lty = 0)
      legend("topright", cex = 1.8, box.lty = 0, bg=NULL,
             c("Estimated by allometric models"))

      # The changes in carbon storage across tree age
      with(Mydata2, plot(stemdj, total_C, type="b", lty=1, lwd=2,pch = 19,
                         col = "blue", cex.axis=1.5,cex.lab=1.8,cex = 2,
                         xlim = c(0, max(stemdj)+Ageclass),
                         ylim = c(0, 1.2*(max(max(total_C)))),
                         xlab = "Age (years)",
                         ylab = "Total tree C storage (kg)", las=1))
      legend("topleft", "(b)", cex = 1.8, box.lty = 0)
      legend("topright", cex = 1.8, box.lty = 0, bg=NULL,
             c("Estimated by allometric models"))
    }



    # Check if there a BEFdata exists
    if(!missing(BEFdata)){

      # Calculate tree biomass using volume model
      stem_bio <- ifelse(Volume!=0,Volume*BEFdata$WD*1000,Volume)
      stem_bio <- as.numeric(stem_bio)

      aboveground_bio <- ifelse(Volume!=0,Volume*BEFdata$WD*BEFdata$BEF*1000,Volume)
      aboveground_bio <- as.numeric(aboveground_bio)

      belowground_bio <- ifelse(Volume!=0,Volume*BEFdata$WD*BEFdata$BEF*BEFdata$R*1000,Volume)
      belowground_bio <- as.numeric(belowground_bio)

      total_bio <- ifelse(Volume!=0,Volume*BEFdata$WD*BEFdata$BEF*(1+BEFdata$R)*1000,Volume)
      total_bio <- as.numeric(total_bio)

      # Calculate tree carbon storage
      stem_Carbon <-  stem_bio*BEFdata$Cconcentration
      stem_Carbon <- as.numeric(stem_Carbon)

      aboveground_Carbon <- aboveground_bio*BEFdata$Cconcentration
      aboveground_Carbon <- as.numeric(aboveground_Carbon)

      belowground_Carbon <- belowground_bio*BEFdata$Cconcentration
      belowground_Carbon <- as.numeric(belowground_Carbon)

      total_Carbon <- total_bio*BEFdata$Cconcentration
      total_Carbon <- as.numeric(total_Carbon)


      # Construct data frame that contains tree biomass and tree carbon storage
      stemdj <- as.numeric(stemdj)

      volumeCarbon <- cbind(stemdj, stem_bio, aboveground_bio,
                            belowground_bio, total_bio, stem_Carbon,
                            aboveground_Carbon, belowground_Carbon,
                            total_Carbon)
      volumeCarbon <- as.data.frame(volumeCarbon)

      write.table(volumeCarbon, file ="../volumeCarbon.csv", sep =",")

      # Combined tree growth trend data and carbon storage data
      Mydata3 <- data.frame(StemGrowth, volumeCarbon)


      # Open a new graph frame
      dev.new(title =
                "tree biomass and carbon storage estimated by volume model",
              width = 4500, height = 2250, noRStudioGD = TRUE)

      # Define the parameters of graph
      mat <- matrix(c(1,2),1,2,byrow=FALSE)
      # mat
      layout(mat)
      layout.show(2)
      par(oma=c(0,0,0,0),mar=c(5,5,5,1))

      # Draw a graph for tree biomass and carbon storage across tree age
      # The changes in tree biomass across tree age
      with(Mydata3, plot(stemdj, total_bio, type="b", lty=1, lwd=2, pch = 19,
                         col = "blue", cex.axis=1.5, cex.lab=1.8,cex = 2,
                         xlim = c(0, max(stemdj)+Ageclass),
                         ylim = c(0, 1.2*(max(total_bio))),
                         xlab = "Age (years)",
                         ylab = "Total tree biomass (kg)",las=1))
      legend("topleft", "(a)", cex = 1.8, box.lty = 0)
      legend("topright", cex = 1.8, box.lty = 0, bg=NULL,
             c("Estimated by volume model"))

      # The changes in carbon storage across tree age
      with(Mydata3, plot(stemdj, total_Carbon, type="b", lty=1, lwd=2,pch = 19,
                         col = "blue", cex.axis=1.5,cex.lab=1.8,cex = 2,
                         xlim = c(0, max(stemdj)+Ageclass),
                         ylim = c(0, 1.2*(max(max(total_Carbon)))),
                         xlab = "Age (years)",
                         ylab = "Total tree C storage (kg)", las=1))
      legend("topleft", "(b)", cex = 1.8, box.lty = 0)
      legend("topright", cex = 1.8, box.lty = 0, bg=NULL,
             c("Estimated by volume model"))
    }
  }

  # Whether to call the treecarbon graph
  if (treecarbon == TRUE) {
    calCarbon(parameterdata, BEFdata)
  } else {
    print("The user does not call the treecarbon graph")
  }




  ###############################################################
  # Define allomHD() function
  allomHD <- function() {
    # Open a new graph frame
    dev.new(title = "height-diameter relationships",
            width = 4500, height = 2250, noRStudioGD = TRUE)

    # Define the parameters of graph
    mat <- matrix(c(1, 2), 1, 2, byrow = FALSE)
    # mat
    layout(mat)
    layout.show(2)
    par(oma = c(0, 0, 0, 0), mar = c(5, 5, 1, 1))

    # Extract sub-datasets with tree height greater than 1.3 m
    Mydata <- subset(StemGrowth, StemGrowth$Height >= 1.3)

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
    lines(xv1_gomperz, yv1_gomperz, col = "gold", lwd = 3)
    SSre_gomperz <- sum(residuals(HDgomperz)^2)
    SStot_gomperz <- sum((Mydata$Height - mean(Mydata$Height))^2)
    R2_gomperz <- 1 - SSre_gomperz / SStot_gomperz
    R2_gomperz
    N <- nrow(Mydata)
    RMSE_gomperz <- sqrt(sum(residuals(HDgomperz)^2) / (N - 2))
    RMSE_gomperz

    legend("topright", "(a)", cex = 1.2, box.lty = 0)
    legend("topleft",
           cex = 1.5, lwd = 2.5,
           legend = c("HDrichards", "HDlogistic", "HDgomperz", "HDweibull"),
           lty = 1, col = c("orangered", "blue", "forestgreen", "gold"),
           title = "Models"
    )

    # Filter the model providing the best fit
    min_Res <- min(anova(HDrichards, HDlogistic, HDgomperz, HDweibull)[2])

    HDrichards_Res <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[1, 2]
    HDlogistic_Res <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[2, 2]
    HDweibull_Res <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[3, 2]
    HDgomperz_Res <- anova(HDrichards, HDlogistic, HDweibull, HDgomperz)[4, 2]

    if (min_Res == HDrichards_Res) {
      plot(Mydata$DBHt, Mydata$Height,
           pch = 21, bg = "purple",
           cex.axis = 1.5, cex.lab = 1.8, cex = 3,
           xlim = c(0, 1.1 * max(Mydata$DBHt)),
           ylim = c(0, 1.1 * max(Mydata$Height)),
           xlab = "Tree DBH (cm)", ylab = "Tree height (m)", las = 1
      )
      lines(xv1_richards, yv1_richards, col = "orangered", lwd = 3)
      legend("topright", "(b)", cex = 1.2, box.lty = 0)
      legend("topleft", cex = 1.5,
             legend = c(expression(paste(italic(R)^2 == "")),
                        expression(paste(italic(RMSE) == "")),
                        round(R2_richards, 3), round(RMSE_richards, 3)),
             title = "Associated statistics", ncol = 2)
      legend("bottomright",
             cex = 1.5,
             legend = c(
               "a =", "b =", "c =",
               round(summary(HDrichards)$parameters[1], 3),
               round(summary(HDrichards)$parameters[2], 3),
               round(summary(HDrichards)$parameters[3], 3)
             ),
             title = "HDrichards Parameters",
             ncol = 2
      )
    } else if (min_Res == HDlogistic_Res) {
      plot(Mydata$DBHt, Mydata$Height,
           pch = 21, bg = "purple",
           cex.axis = 1.5, cex.lab = 1.8, cex = 3,
           xlim = c(0, 1.1 * max(Mydata$DBHt)),
           ylim = c(0, 1.1 * max(Mydata$Height)),
           xlab = "Tree DBH (cm)", ylab = "Tree height (m)", las = 1
      )
      lines(xv1_logistic, yv1_logistic, col = "blue", lwd = 3)
      legend("topright", "(b)", cex = 1.2, box.lty = 0)
      legend("topleft", cex = 1.5,
             legend = c(expression(paste(italic(R)^2 == "")),
                        expression(paste(italic(RMSE) == "")),
                        round(R2_logistic, 3), round(RMSE_logistic, 3)),
             title = "Associated statistics", ncol = 2)
      legend("bottomright",
             cex = 1.5,
             legend = c(
               "a =", "b =", "c =",
               round(summary(HDlogistic)$parameters[1], 3),
               round(summary(HDlogistic)$parameters[2], 3),
               round(summary(HDlogistic)$parameters[3], 3)
             ),
             title = "HDlogistic Parameters",
             ncol = 2
      )
    } else if (min_Res == HDgomperz_Res) {
      plot(Mydata$DBHt, Mydata$Height,
           pch = 21, bg = "purple",
           cex.axis = 1.5, cex.lab = 1.8, cex = 3,
           xlim = c(0, 1.1 * max(Mydata$DBHt)),
           ylim = c(0, 1.1 * max(Mydata$Height)),
           xlab = "Tree DBH (cm)", ylab = "Tree height (m)", las = 1
      )
      lines(xv1_gomperz, yv1_gomperz, col = "gold", lwd = 3)
      legend("topright", "(b)", cex = 1.2, box.lty = 0)
      legend("topleft",
             cex = 1.5,
             legend = c(expression(paste(italic(R)^2 == "")),
                        expression(paste(italic(RMSE) == "")),
                        round(R2_gomperz, 3), round(RMSE_gomperz, 3)),
             title = "Associated statistics", ncol = 2
      )
      legend("bottomright",
             cex = 1.5,
             legend = c(
               "a =", "b =", "c =",
               round(summary(HDgomperz)$parameters[1], 3),
               round(summary(HDgomperz)$parameters[2], 3),
               round(summary(HDgomperz)$parameters[3], 3)
             ),
             title = "HDgomperz Parameters",
             ncol = 2
      )
    } else if (min_Res == HDweibull_Res) {
      plot(Mydata$DBHt, Mydata$Height,
           pch = 21, bg = "purple",
           cex.axis = 1.5, cex.lab = 1.8, cex = 3,
           xlim = c(0, 1.1 * max(Mydata$DBHt)),
           ylim = c(0, 1.1 * max(Mydata$Height)),
           xlab = "Tree DBH (cm)", ylab = "Tree height (m)", las = 1
      )
      lines(xv1_weibull, yv1_weibull, col = "gold", lwd = 3)
      legend("topright", "(b)", cex = 1.2, box.lty = 0)
      legend("topleft",
             cex = 1.5,
             legend = c(expression(paste(italic(R)^2 == "")),
                        expression(paste(italic(RMSE) == "")),
                        round(R2_weibull, 3), round(RMSE_weibull, 3)),
             title = "Associated statistics", ncol = 2
      )
      legend("bottomright",
             cex = 1.5,
             legend = c(
               "a =", "b =", "c =",
               round(summary(HDweibull)$parameters[1], 3),
               round(summary(HDweibull)$parameters[2], 3),
               round(summary(HDweibull)$parameters[3], 3)
             ),
             title = "HDweibull parameters",
             ncol = 2
      )
    }
  }

  # Whether to call the HDmodel graph
  if (HDmodel == TRUE) {
    allomHD()
  } else {
    print("The user does not call the HDmodel graph")
  }
}
