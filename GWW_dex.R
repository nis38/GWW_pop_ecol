#### dex ####
#PACKAGES
if (!require(XML)) install.packages("XML")
if (!require(methods)) install.packages("methods")
if (!require(plyr)) install.packages("plyr")
library(XML)
library(methods)
library(plyr)

dex <- function(filepath) {
  #read in SVG
  p <- xmlParse(filepath) #reads SVG
  d <- xmlRoot(p) 
  g <- d[-c(1,2)] # extract layers
  
  h <- as.numeric(gsub("mm", "", xmlGetAttr(d, "height")))
  w <- as.numeric(gsub("mm", "", xmlGetAttr(d, "width")))
  
  # Extract all Specimen IDs ---------------------------------------------------
  
  #get path ids
  f <- lapply(g, function(x) {x["ellipse"]}) # each element(layer) is a list of ellipses
  m <- lapply(f, function(x){lapply(x, xmlGetAttr, "id")})
  ids <- unname(unlist(m)) # vector of path id for every ellipse

  # Assign Taxon Names to  Specimen IDs --------------------------------------
  
  # get taxon name from the layer names
  m <- lapply(g, xmlGetAttr, "inkscape:label") # list of taxon names
  
  # number of ellipses (specimens) for each taxa (layer)
  n <- lapply(f, length) # list of specimen numbers

  # replicate taxa for number of specimens
  taxa <- unname(unlist(mapply(rep, m, times = n)))# vector of taxa corresponding to each path id

  # Extract Coordinates and Ellipse dimensions --------------------------------------
  
  # extract x, y, rx, ry
  get_discxy <- function(X) {
    x <- as.numeric(xmlGetAttr(X, "cx"))
    y <- as.numeric(xmlGetAttr(X, "cy"))
    ifelse(xmlName(X) == "ellipse", { 
      rx <- as.numeric(xmlGetAttr(X, "rx"))
      ry <- as.numeric(xmlGetAttr(X, "ry"))
    }, {
      rx <- as.numeric(xmlGetAttr(X, "r"))
      ry <- as.numeric(xmlGetAttr(X, "r"))
    })
    A <- pi * rx * ry
    return(cbind(x, y, rx, ry, A))# may be slower to covert to dataframe
  }
  
  e <- lapply(unlist(f), get_discxy) # list of x, y, rx, ry for each layer (taxon)
  coords <- lapply(e, function(x){x[,1:2]})
  coords <- ldply(coords)[,-1]
  
  # Transform coordinate values --------------------------------------
  
  # extract transformation values 
  get_transform <- function(x) {
    t <- xmlGetAttr(x, "transform")
    if (!is.null(t)) {
      t <- gsub("[\\(\\)]", "", regmatches(t, gregexpr("\\(.*?\\)", t))[[1]])
      t <- strsplit(t, ",")
      t <- as.numeric(t[[1]])
      
    }
    return(t)
  }
  
  get_transform_type <- function(x) {
    t <- xmlGetAttr(x, "transform")
    if (!is.null(t)) {
      t <- sub("[^[:alpha:]]+", "", t)
    }
    return(t)
  }
  
  #SVG transforms
  matrix_transform <- function(x, y, t){
    x1 <- t[1] * x + t[3] * y + t[5]
    y <- t[2] * x + t[4] * y + t[6]
    return(c(x1, y))
  }
  
  rotate_transform <- function(x, y, t){
    l <- length(t)
    
    ifelse(l == 1, {
      t <- t * (pi / 180)
      x1 <- cos(t) * x - sin(t) * y
      y <- sin(t) * x + cos(t) * y
    },
    {
      ta <- t[1] * (pi / 180)
      x1 <- cos(ta) * t[2] - sin(ta) * -t[3]
      y <- sin(ta) * t[2] + cos(ta) * -t[3]
    })
    
    return(c(x1, y))
  }
  
  translate_transform <- function(x, y, t){
    x1 <- x + t[1]
    y <- y + t[2]
    
    return(c(x1, y))
  }
  
  
  #apply transformation - will add more options/warnings
  apply_transform <- function(x, y, values, type) {
    
    if (!is.null(type)) {
      ifelse(type == "translate", {
        res <- translate_transform(x, y, values)
      },
      {ifelse(type == "matrix", {
        res <- matrix_transform(x, y, values)
      },
      {ifelse(type == "rotate", {
        res <- rotate_transform(x, y, values)
      },
      {res <- c(x,y)})
      })
      })
    }
    else{
      res <- c(x,y)
    }
    
    return(res)
  }
  
  
  # Ellipse transform - rightmost/lowermost transform first
  te <- lapply(unlist(f), get_transform) # list of ellipse transform values
  type_e <- lapply(unlist(f), get_transform_type) # list of ellipse transform types, use table(type) to see all used

  xy_trans <- as.data.frame(t(mapply(apply_transform, coords$x, coords$y, te, type_e)))
  
  # Layer transformation
  t <- lapply(unlist(g), get_transform) # List of transform values for each layer
  type <- lapply(unlist(g), get_transform_type) # List of transform types for each layer

  # duplicate so transform specified for each ellipse
  typed <- unname(unlist(mapply(function(x, times){
    ifelse(!is.null(x), {
      return(rep(x,times))
    }, {return(vector(mode = "character", length = times))})
  }, type, times = n)))

  td <- list()
  for(i in 1:length(t)){
    td[[i]] <- rep(list(t[[i]]), n[[i]])
  }


  td <- do.call(c, td)

  #apply transform
  xy <- as.data.frame(t(mapply(apply_transform, xy_trans[,1], xy_trans[,2], td, typed)))
  names(xy) <- c("x", "y")
  
  #account for inkscape window axes
  xy$y <- - xy$y + h
  
  # Combine to final data frame  --------------------------------------
  dims <- ldply(e)[,-1]
  df <- cbind.data.frame(taxa, ids, xy, rx = dims$rx, ry = dims$ry, area = dims$A)
  
  plot(y ~ x, df)
  
  return(list(df, h, w))
  
}
