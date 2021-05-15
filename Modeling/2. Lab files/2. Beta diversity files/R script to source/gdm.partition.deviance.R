#' @title Perform deviance partitioning for up to two variable sets
#' (e.g., soils and climate), plus (optionally) space.
#' @rdname gdm.partition.deviance
#' @name gdm.partition.deviance
#' @description Partitions deviance...
#' @param sitePairTable A correctly formatted site-pair table from
#' \code{\link[gdm]{formatsitepair}.
#' @param varSets A list in which each element is a vector of variable names
#' across which deviance partitioning is to be performed, excluding
#' geographic distance (which is set by the partSpace argument). Variable names
#' must match those used to build the site-pair table. See example.
#' @param partSpace Whether or not to perform the partitioning using
#' geographic space.
#' @return A dataframe summarizing partitioning results.
#' @export
#' @author Matt Fitzpatrick and Karel Mokany
#' 
#' @examples
#' # set up site-pair table
#' load(system.file("./data/gdm.RData", package="gdm"))
#' sppData <- gdmExpData[c(1,2,13,14)]
#' envTab <- gdmExpData[c(2:ncol(gdmExpData))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", 
#' siteColumn="site", predData=envTab) 
#' 
#' # Make list of variable sets for partitioning
#' varSet <- vector("list", 2)
#' varSet[[1]] <- c("awcA", "phTotal", "sandA", "shcA", "solumDepth")
#' varSet[[2]] <- c("bio5", "bio6", "bio15", "bio18", "bio19")
#' 
#' # name sets
#' names(varSet) <- c("soil", "climate")
#' 
#' # partition soils, climate, and space
#' scgPart <- gdm.partition(sitePairTab, varSet, partSpace=T)
#' 
#' # partition soils and climate only
#' scPart <- gdm.partition(sitePairTab, varSet, partSpace=F)
#' 
gdm.partition.deviance <- function(sitePairTable, varSets=list(), partSpace=T){
  
  #require(VennDiagram)
  #sitePairTable <- sitePairTab
  #varSets <- lll
  #partSpace <- T
  
  if(class(varSets)!="list"){
    stop("Variables names must be given in a list.")
  }
  
  if(!"gdmData" %in% class(sitePairTable)){
    stop("Site-pair table must be class=gdmData.")
  }
  
  # number of variables sets to partition across,
  # excluding geo
  partLength <- length(varSets)
  
  # make a third 'total' variable set by combining
  # the two provided (assuming two are provided, otherwise, skip)
  if(partLength==2){
    varSets[[length(varSets)+1]] <- as.character(unlist(varSets))
    names(varSets)[partLength+1] <- paste(names(varSets)[1:partLength],
                                          collapse = " and ")}
  
  # if partitioning using geo, setup to fit models w/ w/o geo
  if(partSpace){
    partGeo <- c(FALSE, TRUE)
  } else {partGeo <- FALSE}
  
  # objects to hold results
  dev <- NULL
  variableSets <- NULL
  
  # Work through each partition
  # Start with geo=T, as needed
  # then select variable sets from site-pair table & fit GDM to each
  if(partSpace){
    #print("geo = ")
    dev <- gdm(sitePairTable[,1:6], geo=T)$explained
    variableSets <- ("geo") 
    # loop through all combinations of geo & other variable sets
    for(g in partGeo){
      for(v in 1:length(varSets)){
        #print(paste("geo = ", g, " and ", names(varSets)[v], sep=""))
        # select variable sets
        greppers <- c(paste("s1.", varSets[[v]], sep=""), 
                      paste("s2.", varSets[[v]], sep=""))
        modSPT <- sitePairTable[, c(1:6,which(names(sitePairTable) %in% greppers))]
        # fit GDM, get deviance
        dev <- c(dev, gdm(modSPT, geo=g)$explained)
        if(g){
          ggg <- "geo and "
        } else {
          ggg <- ""
        }
        variableSets <- c(variableSets, paste(ggg, names(varSets)[v], sep=""))
      }
    }
    
    # data frame with results
    devTable <- data.frame(VARIABLE_SET = variableSets, DEVIANCE = dev,
                           stringsAsFactors = F)
    devTable[7,1] <- paste("ALL VARIABLES (", devTable[7,1], ")", sep="")
    
    # partition deviance
    if(length(varSets)==1){
      devTable[4,1] <- paste(names(varSets), "alone", sep=" ")
      devTable[4,2] <- devTable[3,2]-devTable[1,2]
      devTable[5,1] <- "geo alone"
      devTable[5,2] <- devTable[3,2]-devTable[2,2]
      devTable[6,1] <- "Unexplained"
      devTable[6,2] <- 100-devTable[3,2]
    }
    
    if(length(varSets)==3){
      devTable <- devTable[c(3,2,1,4:7),]
      devTable[8,1] <- paste(names(varSets)[2], "alone", sep=" ")
      devTable[8,2] <- devTable[7,2]-devTable[5,2]
      devTable[9,1] <- paste(names(varSets)[1], "alone", sep=" ")
      devTable[9,2] <- devTable[7,2]-devTable[6,2]
      devTable[10,1] <- "geo alone"
      devTable[10,2] <- devTable[7,2]-devTable[4,2]
      devTable[11,1] <- paste(names(varSets)[1], " union ", names(varSets)[2],
                              ", exclude geo", sep="")
      devTable[11,2] <- devTable[1,2]-devTable[8,2]-((devTable[1,2]+devTable[3,2])-devTable[6,2])
      devTable[12,1] <- paste(names(varSets)[1], " union geo",
                              ", exclude ", names(varSets)[2], sep="")
      devTable[12,2] <- devTable[2,2]-devTable[9,2]-((devTable[1,2]+devTable[2,2])-devTable[4,2])
      devTable[13,1] <- paste("geo union ", names(varSets)[2],
                              ", exclude ", names(varSets)[1], sep="")
      devTable[13,2] <- devTable[3,2]-devTable[10,2]-((devTable[2,2]+devTable[3,2])-devTable[5,2])
      devTable[14,1] <- paste(names(varSets)[2], " union ", names(varSets)[1],
                              " union geo", sep="")
      devTable[14,2] <- devTable[7,2]-devTable[8,2]-devTable[9,2]-devTable[10,2]-devTable[11,2]-devTable[12,2]-devTable[13,2]
      devTable[15,1] <- "Unexplained"
      devTable[15,2] <- 100-devTable[7,2]
    }
  }
  
  # same as above, but excluding geo if requested
  if(partSpace==F){
    for(g in partGeo){
      for(v in 1:length(varSets)){
        #print(paste("geo = ", g, " and ", names(varSets)[v], sep=""))
        greppers <- c(paste("s1.", varSets[[v]], sep=""), 
                      paste("s2.", varSets[[v]], sep=""))
        modSPT <- sitePairTable[, c(1:6,which(names(sitePairTable) %in% greppers))]
        
        dev <- c(dev, gdm(modSPT, geo=g)$explained)
        
        if(g){
          ggg <- "geo and "
        } else {
          ggg <- ""
        }
        variableSets <- c(variableSets, paste(ggg, names(varSets)[v], sep=""))
      }
    }
    
    devTable <- data.frame(variableSet = variableSets, deviance = dev,
                           stringsAsFactors = F)
    devTable[3,1] <- paste("ALL VARIABLES (", devTable[3,1], ")", sep="")
    
    devTable[4,1] <- paste(names(varSets)[1], "alone", sep=" ")
    devTable[4,2] <- devTable[3,2]-devTable[2,2]
    devTable[5,1] <- paste(names(varSets)[2], "alone", sep=" ")
    devTable[5,2] <- devTable[3,2]-devTable[1,2]
    devTable[6,1] <- paste(names(varSets)[1], " union ", names(varSets)[2],sep="")
    devTable[6,2] <- devTable[3,2]-devTable[4,2]-devTable[5,2]
    devTable[7,1] <- "Unexplained"
    devTable[7,2] <- 100-devTable[3,2]
    
    # vennPlot <- draw.pairwise.venn(area1=round(devTable[1,2], 2), 
    #                                    area2=round(devTable[2,2], 2), 
    #                                    cross.area = round(devTable[6,2], 2),
    #                                    category=c(names(varSets)[1],
    #                                               names(varSets)[2]),
    #                                    euler.d=T, cat.col="darkred",
    #                                scaled=T)
    # grid.newpage()
    #     grid.draw(vennPlot)
    
  }
  devTable[,2] <- round(devTable[,2],3)
  return(devTable)
}