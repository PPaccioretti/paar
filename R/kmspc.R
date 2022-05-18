#' MULTISPATI-PCA clustering
#'
#'
#' @inheritParams remove_inlier
#' @param data sf object
#' @param variables variables to use for clustering, if missing, all numeric
#' variables will be used
#' @param number_cluster \code{numeric} vector with number of final clusters
#' @param explainedVariance \code{numeric} number in percentage of explained variance
#' from PCA analysis to keep and make cluster process
#' @param fuzzyness A number greater than 1 giving the degree of fuzzification.
#' @param center a logical or numeric value, centring option
#' if TRUE, centring by the mean
#' if FALSE no centring
#' if a numeric vector, its length must be equal to the number of
#' columns of the data frame df and gives the decentring
#' @param distance \code{character} Must be one of the following:
#' If "euclidean", the mean square error, if "manhattan", the mean
#' absolute error is computed. Abbreviations are also accepted.
#' @return a list with classification results and indices to select bes number of
#' clusters.
#' @export
#'
#'
kmspc <- function(data,
                  variables,
                  number_cluster = 3:5,
                  explainedVariance = 70,
                  ldist = 0,
                  udist = 40,
                  center = TRUE,
                  fuzzyness = 1.2,
                  distance = "euclidean",
                  zero.policy = FALSE) {

  if (missing(variables)) {
    warning("All numeric Variables will be used to make clusters",
            call. = FALSE)
    myNumVars <-
      unlist(lapply(sf::st_drop_geometry(data), is.numeric))

    if (sum(myNumVars) == 0) {
      stop('Non numeric variables were found in data')
    }
    variables <- names(data)[myNumVars]
  }

  if (!inherits(data, "sf") & (length(variables) > 1)) {
    stop('data must be an sf object')
  }

  if (between(explainedVariance, 0, 1)) {
    explainedVariance <- explainedVariance * 100
  }

  if (!between(explainedVariance, 0, 100)) {
    stop("explainedVariance must bu a value between 0 and 100")
  }

  data <- data[, variables]
  data <- stats::na.omit(data)
  data_clust <- data
  if (ncol(sf::st_drop_geometry(data)) > 1 ) {
  lw <- spatial_weights(data, ldist, udist, zero.policy = zero.policy)
  pca <-
    dudy_pca(
    # ade4::dudi.pca(
      sf::st_drop_geometry(data),
      center = center,
      scale = TRUE,
      scannf = FALSE,
      nf = length(variables)
    )
  suppressWarnings(
    ms <- multispati(
      # adespatial::multispati(
        pca,
        lw,
        scannf = FALSE,
        nfnega = length(variables),
        nfposi = length(variables)
      )
  )

  invisible(utils::capture.output(resms <- summary(ms)))
  var_ms <- resms[, 2, drop = F]
  nfila_ms <- length(ms$eig)
  propvar_ms <- var_ms / nfila_ms
  propvaracum_ms <- cumsum(propvar_ms) * 100

  eje_ms <- c(1:nfila_ms)
  resultado_ms <-
    data.frame(eje_ms, resms$eig, resms$var, propvar_ms, propvaracum_ms)
  names(resultado_ms) <-
    c("Axis",
      "Eigenvalue",
      "Spatial Variance",
      "Prop",
      "Acum. Prop.")


  num_sPC <-
    seq_len(Position(function(x) {
      x > explainedVariance
    }, resultado_ms[, 5]))

  data_clust <- ms$li[num_sPC]
  }
  if (inherits(data_clust, "sf")) {
    data_clust <-  sf::st_drop_geometry(data_clust)
  }
  make_clasification(data_clust,
                     number_cluster,
                     fuzzyness = fuzzyness,
                     distance = distance)

}






#' Auxiliary function for kmspc
#' @noRd
cmeans_vectorized <-
  function(data,
           nclusters,
           ...,
           index = c(
             "xie.beni",
             # "fukuyama.sugeno",
             "partition.coefficient",
             "partition.entropy"
           )) {
    lapply(nclusters, function(center, x, index, ...) {
      myClusters <- e1071::cmeans(
        x = data,
        centers = center,
        method = "cmeans",
        iter.max = 100,
        ...
      )

      myIndices <- fclustIndex(myClusters,
                               x = data,
                               index = index)

      list("cluster" = myClusters,
           "indices" = myIndices)
    },
    x = data,
    index = index,
    ...)


  }




#' Auxiliary function for kmspc
#' @noRd
spatial_weights <- function(data, ldist, udist, zero.policy = FALSE) {
  oldzp <- spdep::get.ZeroPolicyOption()
  spdep::set.ZeroPolicyOption(zero.policy)

  gri <-
    spdep::dnearneigh(data,
                      ldist,
                      udist)
  lw <- tryCatch(
    {spdep::nb2listw(gri, style = "W")},
    error = function(e) {
      if (agrepl("Empty neighbour sets found", e)) {
        stop("Empty neighbour sets found", call. = FALSE)
      }
    },
    silent = TRUE
  )

  spdep::set.ZeroPolicyOption(oldzp)
  lw
}


#' Auxiliary function for kmspc
#' @noRd
between <- function(x, min, max) {
  x > min & x < max
}

#' Auxiliary function for kmspc
#' @noRd
normalize <- function(x) {
  if (length(x) > 1) {
    return(x / max(x))
  } else {x}

}




#' Auxiliary function for kmspc
#' @noRd
summarize_indices <- function(indices, number_cluster) {

  indices <- do.call("rbind", indices)

  IndN <- indices
  if (nrow(indices) > 1) {
    IndN <- apply(indices, 2, normalize)
    IndN <- apply(IndN, 1, function(xx) {
      sqrt(sum(xx ^ 2))
    })
  }

  indicesresults <- data.frame(number_cluster, indices, IndN)
  names(indicesresults) = c(
    "Num. Cluster",
    "Xie Beni",
    # "Fukuyama Sugeno",
    "Partition Coefficient",
    "Entropy of Partition",
    "Summary Index"
  )
  indicesresults
}

#' Auxiliary function for kmspc
#' @noRd
summarize_clusters_metrics <- function(clasifications, number_cluster) {
  res_iter <- lapply(clasifications, '[[', "iter")
  res_scdd <- lapply(clasifications, '[[', "withinerror")

  data.frame("Clusters" = number_cluster,
             "Iterations" = unlist(res_iter),
             "SSDW" = unlist(res_scdd))
}


#' Auxiliary function for kmspc
#' @noRd
summarize_clusters <- function(clasifications, number_cluster) {
  res_clas <-
    lapply(clasifications, '[[', "cluster")
  res_clas <- data.frame(res_clas)
  names(res_clas) <-
    paste("Cluster", "_",
          number_cluster,
          sep = "")
  res_clas
}



#' Auxiliary function for kmspc
#' @noRd
#'
make_clasification <- function(data, number_cluster, fuzzyness, distance) {
  clasifications <- cmeans_vectorized(data,
                                       number_cluster,
                                       dist = distance,
                                       m = fuzzyness)
  indices <- lapply(clasifications, '[[', "indices")
  indices <- summarize_indices(indices, number_cluster)


  clasifications <- lapply(clasifications, '[[', "cluster")




  clasifResults <-
    summarize_clusters_metrics(clasifications, number_cluster)
  FinalCluster <-
    summarize_clusters(clasifications, number_cluster)

  list("summaryResults" = clasifResults,
       "indices" = indices,
       "cluster" = FinalCluster)
}






















#' @noRd
dudy_pca <- function(df, row.w = rep(1, nrow(df))/nrow(df), col.w = rep(1,
                                                                         ncol(df)), center = TRUE, scale = TRUE, scannf = TRUE, nf = 2)
{
  df <- as.data.frame(df)
  nc <- ncol(df)
  if (any(is.na(df)))
    stop("na entries in table")
  f1 <- function(v) sum(v * row.w)/sum(row.w)
  f2 <- function(v) sqrt(sum(v * v * row.w)/sum(row.w))
  if (is.logical(center)) {
    if (center) {
      center <- apply(df, 2, f1)
      df <- sweep(df, 2, center)
    }
    else center <- rep(0, nc)
  }
  else if (is.numeric(center) && (length(center) == nc))
    df <- sweep(df, 2, center)
  else stop("Non convenient selection for center")
  if (scale) {
    norm <- apply(df, 2, f2)
    norm[norm < 1e-08] <- 1
    df <- sweep(df, 2, norm, "/")
  }
  else norm <- rep(1, nc)
  X <- as_dudi(df, col.w, row.w, scannf = scannf, nf = nf,
               call = match.call(), type = "pca")
  X$cent <- center
  X$norm <- norm
  X
}


#' @noRd
as_dudi <-
  function(df,
           col.w,
           row.w,
           scannf,
           nf,
           call,
           type,
           tol = 1e-07,
           full = FALSE)
  {
    if (!is.data.frame(df))
      stop("data.frame expected")
    lig <- nrow(df)
    col <- ncol(df)
    if (length(col.w) != col)
      stop("Non convenient col weights")
    if (length(row.w) != lig)
      stop("Non convenient row weights")
    if (any(col.w < 0))
      stop("col weight < 0")
    if (any(row.w < 0))
      stop("row weight < 0")
    if (full)
      scannf <- FALSE
    transpose <- FALSE
    if (lig < col)
      transpose <- TRUE
    res <- list(tab = df, cw = col.w, lw = row.w)
    df <- as.matrix(df)
    df.ori <- df
    df <- df * sqrt(row.w)
    df <- sweep(df, 2, sqrt(col.w), "*")
    if (!transpose) {
      df <- crossprod(df, df)
    }
    else {
      df <- tcrossprod(df, df)
    }
    eig1 <- eigen(df, symmetric = TRUE)
    eig <- eig1$values
    rank <- sum((eig / eig[1]) > tol)
    # if (scannf) {
    #   # if (exists("ade4TkGUIFlag")) {
    #   #   nf <- ade4TkGUI::chooseaxes(eig, rank)
    #   # }
    #   # else {
    #     # barplot(eig[1:rank])
    #     # cat("Select the number of axes: ")
    #     # nf <- as.integer(readLines(n = 1))
    #     # # messageScannf(call, nf)
    #   }
    # }
    if (nf <= 0)
      nf <- 2
    if (nf > rank)
      nf <- rank
    if (full)
      nf <- rank
    res$eig <- eig[1:rank]
    res$rank <- rank
    res$nf <- nf
    col.w[which(col.w == 0)] <- 1
    row.w[which(row.w == 0)] <- 1
    dval <- sqrt(res$eig)[1:nf]
    if (!transpose) {
      col.w <- 1 / sqrt(col.w)
      auxi <- eig1$vectors[, 1:nf] * col.w
      auxi2 <- sweep(df.ori, 2, res$cw, "*")
      auxi2 <- data.frame(auxi2 %*% auxi)
      auxi <- data.frame(auxi)
      names(auxi) <- paste("CS", (1:nf), sep = "")
      row.names(auxi) <- make.unique(names(res$tab))
      res$c1 <- auxi
      names(auxi2) <- paste("Axis", (1:nf), sep = "")
      row.names(auxi2) <- row.names(res$tab)
      res$li <- auxi2
      res$co <- sweep(res$c1, 2, dval, "*")
      names(res$co) <- paste("Comp", (1:nf), sep = "")
      res$l1 <- sweep(res$li, 2, dval, "/")
      names(res$l1) <- paste("RS", (1:nf), sep = "")
    } else {
      row.w <- 1 / sqrt(row.w)
      auxi <- eig1$vectors[, 1:nf] * row.w
      auxi2 <- t(sweep(df.ori, 1, res$lw, "*"))
      auxi2 <- data.frame(auxi2 %*% auxi)
      auxi <- data.frame(auxi)
      names(auxi) <- paste("RS", (1:nf), sep = "")
      row.names(auxi) <- row.names(res$tab)
      res$l1 <- auxi
      names(auxi2) <- paste("Comp", (1:nf), sep = "")
      row.names(auxi2) <- make.unique(names(res$tab))
      res$co <- auxi2
      res$li <- sweep(res$l1, 2, dval, "*")
      names(res$li) <- paste("Axis", (1:nf), sep = "")
      res$c1 <- sweep(res$co, 2, dval, "/")
      names(res$c1) <- paste("CS", (1:nf), sep = "")
    }
    res$call <- call
    class(res) <- c(type, "dudi")
    return(res)
  }







#' @noRd
multispati <-
  function(dudi,
           listw,
           scannf = TRUE,
           nfposi = 2,
           nfnega = 0) {
    if (!inherits(dudi, "dudi"))
      stop("object of class 'dudi' expected")
    if (!inherits(listw, "listw"))
      stop("object of class 'listw' expected")
    if (listw$style != "W")
      stop("object of class 'listw' with style 'W' expected")
    NEARZERO <- 1e-14

    dudi$cw <- dudi$cw
    fun <- function(x) {
      spdep::lag.listw(listw, x, TRUE)}
    tablag <- apply(dudi$tab, 2, fun)
    covar <- t(tablag) %*% as.matrix((dudi$tab * dudi$lw))
    covar <- (covar + t(covar)) / 2
    covar <- covar * sqrt(dudi$cw)
    covar <- t(t(covar) * sqrt(dudi$cw))
    covar <- eigen(covar, symmetric = TRUE)
    res <- list()
    res$eig <- covar$values[abs(covar$values) > NEARZERO]
    ndim <- length(res$eig)
    covar$vectors <- covar$vectors[, abs(covar$values) > NEARZERO]

    # if (scannf) {
    #   barplot(res$eig)
    #   cat("Select the number of axes with positive spatial autocorrelation: ")
    #   nfposi <- as.integer(readLines(n = 1))
    #
    #   cat("Select the number of axes with negative spatial autocorrelation: ")
    #   nfnega <- as.integer(readLines(n = 1))
    # }

    if (nfposi <= 0)
      nfposi <- 1
    if (nfnega <= 0)
      nfnega <- 0

    if (nfposi > sum(res$eig > 0)) {
      nfposi <- sum(res$eig > 0)
      warning(paste("There are only", sum(res$eig > 0), "positive factors."))
    }
    if (nfnega > sum(res$eig < 0)) {
      nfnega <- sum(res$eig < 0)
      warning(paste("There are only", sum(res$eig < 0), "negative factors."))
    }

    res$nfposi <- nfposi
    res$nfnega <- nfnega
    agarder <-
      c(1:nfposi, if (nfnega > 0)
        (ndim - nfnega + 1):ndim
        else
          NULL)
    dudi$cw[which(dudi$cw == 0)] <- 1
    auxi <- data.frame(covar$vectors[, agarder] / sqrt(dudi$cw))
    names(auxi) <- paste("CS", agarder, sep = "")
    row.names(auxi) <- names(dudi$tab)
    res$c1 <- auxi
    auxi <- as.matrix(auxi) * dudi$cw
    auxi1 <- as.matrix(dudi$tab) %*% auxi
    auxi1 <- data.frame(auxi1)
    names(auxi1) <- names(res$c1)
    row.names(auxi1) <- row.names(dudi$tab)
    res$li <- auxi1
    auxi1 <- as.matrix(tablag) %*% auxi
    auxi1 <- data.frame(auxi1)
    names(auxi1) <- names(res$c1)
    row.names(auxi1) <-  row.names(dudi$tab)
    res$ls <- auxi1
    auxi <- as.matrix(res$c1) * unlist(dudi$cw)
    auxi <- data.frame(t(as.matrix(dudi$c1)) %*% auxi)
    row.names(auxi) <- names(dudi$li)
    names(auxi) <- names(res$li)
    res$as <- auxi
    res$call <- match.call()
    class(res) <- "multispati"
    return(res)
  }


#' @noRd
summary.multispati <- function(object, ...) {
  norm.w <- function(X, w) {
    f2 <- function(v)
      sum(v * v * w) / sum(w)
    norm <- apply(X, 2, f2)
    return(norm)
  }

  if (!inherits(object, "multispati"))
    stop("to be used with 'multispati' object")

  cat("\nMultivariate Spatial Analysis\n")
  cat("Call: ")
  print(object$call)

  appel <- as.list(object$call)
  dudi <- eval.parent(appel$dudi)
  listw <- eval.parent(appel$listw)

  ## les scores de l'analyse de base
  nf <- dudi$nf
  eig <- dudi$eig[1:nf]
  cum <- cumsum(dudi$eig) [1:nf]
  ratio <- cum / sum(dudi$eig)
  w <-
    apply(dudi$l1,
          2,
          spdep::lag.listw,
          x = listw,
          zero.policy = TRUE)
  moran <- apply(w * as.matrix(dudi$l1) * dudi$lw, 2, sum)
  res <- data.frame(
    var = eig,
    cum = cum,
    ratio = ratio,
    moran = moran
  )
  cat("\nScores from the initial duality diagram:\n")
  print(res)

  ## les scores de l'analyse spatiale
  ## on recalcule l'objet en gardant tous les axes
  eig <- object$eig
  nfposi <- object$nfposi
  nfnega <- object$nfnega
  nfposimax <- sum(eig > 0)
  nfnegamax <- sum(eig < 0)

  ms <- multispati(
    dudi = dudi,
    listw = listw,
    scannf = FALSE,
    nfposi = nfposimax,
    nfnega = nfnegamax
  )

  ndim <- dudi$rank
  nf <- nfposi + nfnega
  agarder <-
    c(1:nfposi, if (nfnega > 0)
      (ndim - nfnega + 1):ndim
      else
        NULL)
  varspa <- norm.w(ms$li, dudi$lw)
  moran <- apply(as.matrix(ms$li) * as.matrix(ms$ls) * dudi$lw, 2, sum)
  res <- data.frame(eig = eig,
                    var = varspa,
                    moran = moran / varspa)

  cat("\nMultispati eigenvalues decomposition:\n")
  print(res[agarder, ])
  return(invisible(res))
}


#' @noRd
print.multispati <- function(x, ...)
{
  cat("Multispati object \n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$nfposi:", x$nfposi, "axis-components saved")
  cat("\n$nfnega:", x$nfnega, "axis-components saved")
  #cat("\n$rank: ")
  #cat(x$rank)
  cat("\nPositive eigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else
    cat("\n")
  cat("Negative eigenvalues: ")
  l0 <- sum(x$eig <= 0)
  cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else
    cat("\n")
  cat('\n')
  sumry <- array("", c(1, 4), list(1, c("vector", "length",
                                        "mode", "content")))
  sumry[1,] <- c('$eig', length(x$eig), mode(x$eig), 'eigen values')

  print(sumry, quote = FALSE)
  cat("\n")
  sumry <-
    array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
  sumry[1,] <-
    c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
  sumry[2,] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
  sumry[3,] <-
    c("$ls", nrow(x$ls), ncol(x$ls), 'lag vector coordinates')
  sumry[4,] <-
    c("$as",
      nrow(x$as),
      ncol(x$as),
      'inertia axes onto multispati axes')


  print(sumry, quote = FALSE)
  cat("other elements: ")
  if (length(names(x)) > 8)
    cat(names(x)[9:(length(names(x)))], "\n")
  else
    cat("NULL\n")
}


## This function is from e1071 package fixed so can be run in R-4.2.0

#' @noRd
fclustIndex <- function(y, x, index = "all")
{
  clres <- y
  gath.geva <- function(clres, x) {
    xrows <- dim(clres$me)[1]
    xcols <- dim(clres$ce)[2]
    ncenters <- dim(clres$centers)[1]
    scatter <- array(0, c(xcols, xcols, ncenters))
    scatternew <- array(0, c(xcols, xcols, ncenters))
    fhv <- as.double(0)
    apd <- as.double(0)
    pd <- as.double(0)
    control <- as.double(0)
    for (i in 1:ncenters) {
      paronomastis <- as.double(0)
      paronomastis2 <- as.double(0)
      for (j in 1:xrows) {
        paronomastis <- paronomastis + clres$me[j, i]
        diff <- x[j,] - clres$ce[i,]
        scatternew[, , i] <- clres$me[j, i] * (t(t(diff)) %*%
                                                 t(diff))
        scatter[, , i] <- scatter[, , i] + scatternew[,
                                                      , i]
      }
      scatter[, , i] <- scatter[, , i] / paronomastis
      for (j in 1:xrows) {
        diff <- x[j,] - clres$ce[i,]
        control <- (t(diff) %*% solve(scatter[, , i])) %*%
          t(t(diff))
        if (control < 1)
          paronomastis2 <- paronomastis2 + clres$me[j,
                                                    i]
      }
      fhv <- fhv + sqrt(det(scatter[, , i]))
      apd <- apd + paronomastis2 / sqrt(det(scatter[, , i]))
      pd <- pd + paronomastis2
    }
    pd <- pd / fhv
    apd <- apd / ncenters
    retval <-
      list(
        fuzzy.hypervolume = fhv,
        average.partition.density = apd,
        partition.density = pd
      )
    return(retval)
  }
  xie.beni <- function(clres) {
    xrows <- dim(clres$me)[1]
    minimum <- -1
    error <- clres$within
    ncenters <- dim(clres$centers)[1]
    for (i in 1:(ncenters - 1)) {
      for (j in (i + 1):ncenters) {
        diff <- clres$ce[i,] - clres$ce[j,]
        diffdist <- t(diff) %*% t(t(diff))
        if (minimum == -1)
          minimum <- diffdist
        if (diffdist < minimum)
          minimum <- diffdist
      }
    }
    xiebeni <- error / (xrows * minimum)
    return(xiebeni)
  }
  fukuyama.sugeno <- function(clres) {
    xrows <- dim(clres$me)[1]
    ncenters <- dim(clres$centers)[1]
    error <- clres$within
    k2 <- as.double(0)
    meancenters <- apply(clres$ce, 2, mean)
    for (i in 1:ncenters) {
      paronomastis3 <- as.double(0)
      for (j in 1:xrows) {
        paronomastis3 <- paronomastis3 + (clres$me[j,
                                                   i] ^ 2)
      }
      diff <- clres$ce[i,] - meancenters
      diffdist <- t(diff) %*% t(t(diff))
      k2 <- k2 + paronomastis3 * diffdist
    }
    fukuyamasugeno <- error - k2
    return(fukuyamasugeno)
  }
  partition.coefficient <- function(clres) {
    xrows <- dim(clres$me)[1]
    partitioncoefficient <- sum(apply(clres$me ^ 2, 1, sum)) / xrows
    return(partitioncoefficient)
  }
  partition.entropy <- function(clres) {
    xrows <- dim(clres$me)[1]
    ncenters <- dim(clres$centers)[1]
    partitionentropy <- 0
    for (i in 1:xrows) {
      for (k in 1:ncenters) {
        if (clres$me[i, k] != 0)
          partitionentropy <- partitionentropy + (clres$me[i,
                                                           k] * log(clres$me[i, k]))
      }
    }
    partitionentropy <- partitionentropy / ((-1) * xrows)
    return(partitionentropy)
  }
  separation.index <- function(clres, x) {
    xrows <- dim(clres$me)[1]
    xcols <- dim(x)[2]
    ncenters <- dim(clres$centers)[1]
    maxcluster <- double(ncenters)
    minimum <- -1
    for (i in 1:ncenters) {
      maxcluster[i] <- max(dist(matrix(x[clres$cl == i],
                                       ncol = xcols)))
    }
    maxdia <- maxcluster[rev(order(maxcluster))[1]]
    for (i in 1:(ncenters - 1)) {
      for (j in (i + 1):(ncenters)) {
        for (m in 1:xrows) {
          if (clres$cl[m] == i) {
            for (l in 1:xrows) {
              if (clres$cl[l] == j) {
                diff <- x[m,] - x[l,]
                diffdist <- sqrt(t(diff) %*% t(t(diff)))
                fraction <- diffdist / maxdia
                if (minimum == -1)
                  minimum <- fraction
                if (fraction < minimum)
                  minimum <- fraction
              }
            }
          }
        }
      }
    }
    return(minimum)
  }
  proportion.exponent <- function(clres) {
    k <- dim(clres$centers)[2]
    xrows <- dim(clres$me)[1]
    bexp <- as.integer(1)
    for (j in 1:xrows) {
      greatint <- as.integer(1 / max(clres$me[j,]))
      aexp <- as.integer(0)
      for (l in 1:greatint) {
        aexp <- aexp + (-1) ^ (l + 1) * (gamma(k + 1) / (gamma(l +
                                                                 1) * gamma(k - l + 1))) * (1 - l * max(clres$me[j,])) ^
          (k - 1)
      }
      bexp <- bexp * aexp
    }
    proportionexponent <- -log(bexp)
    return(proportionexponent)
  }
  index <-
    pmatch(
      index,
      c(
        "gath.geva",
        "xie.beni",
        "fukuyama.sugeno",
        "partition.coefficient",
        "partition.entropy",
        "proportion.exponent",
        "separation.index",
        "all"
      )
    )
  if (any(is.na(index)))
    stop("invalid clustering index")
  if (any(index == -1))
    stop("ambiguous index")
  vecallindex <- numeric(9)
  if (any(index == 1) || any(index == 8)) {
    gd <- gath.geva(clres, x)
    vecallindex[1] <- gd$fuzzy
    vecallindex[2] <- gd$average
    vecallindex[3] <- gd$partition
  }
  if (any(index == 2) || any(index == 8))
    vecallindex[4] <- xie.beni(clres)
  if (any(index == 3) || any(index == 8))
    vecallindex[5] <- fukuyama.sugeno(clres)
  if (any(index == 4) || any(index == 8))
    vecallindex[6] <- partition.coefficient(clres)
  if (any(index == 5) || any(index == 8))
    vecallindex[7] <- partition.entropy(clres)
  if (any(index == 6) || any(index == 8))
    vecallindex[8] <- proportion.exponent(clres)
  if (any(index == 7) || any(index == 8))
    vecallindex[9] <- separation.index(clres, x)
  names(vecallindex) <- c("fhv", "apd", "pd", "xb", "fs", "pc",
                          "pe", "pre", "si")
  if (any(index < 8)) {
    if (any(index == 1))
      vecallindex <- vecallindex[1:3]
    else
      vecallindex <- vecallindex[index + 2]
  }
  return(vecallindex)
}
