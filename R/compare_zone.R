
#' Compare spatial zone means
#'
#' @param data \code{sf} object with zones
#' @param variable \code{character} or \code{sf} object to use for mean comparison
#' @param zonesCol \code{character} colname from data were zone are specified
#' @param alpha \code{numeric} Significance level to use for comparison
#' @param returnLSD \code{logical} when LSD calculates with spatial variance should be returned
#' @param join function to use for st_join if variable is \code{sf} object
#' @param grid_dim \code{numeric} grid dimentins to estimate spatial variance
#'
#' @return \code{list} with differences and descriptive_stat
#' @references
#' Paccioretti, P., Córdoba, M., & Balzarini, M. (2020).
#' FastMapping: Software to create field maps and identify management zones
#' in precision agriculture. Computers and Electronics in Agriculture, 175,
#' 105556 https://doi.org/10.1016/j.compag.2020.105556.
#' @export
#'
#' @example inst/examples/compare_zone.R
#'

compare_zone <- function(data,
                         variable,
                         zonesCol,
                         alpha = 0.05,
                         join = sf::st_nearest_feature,
                         returnLSD = FALSE,
                         grid_dim) {

  descriptive <- varKrigDescr(data,
                              variable,
                              grid_dim = grid_dim)

  comparisson_table <-
    validVarKrig(data,
                 variable,
                 zonesCol,
                 descriptive,
                 returnLSD = returnLSD,
                 join = join)

  list('differences' = comparisson_table$Diferencias,
       'descriptive_stat' = comparisson_table$Descriptivo)

}

#' Auxiliary function for zone compare
#'
#' @noRd
validVarKrig <-
  function(clusterData,
           dataToCompare,
           numCluster,
           descStat,
           join = sf::st_nearest_feature,
           returnLSD) {

    MedCV <-
      paste0(round(descStat$MediaVar, 2), " (", round(descStat$CVVar, 1), ")")
    MedCV <-
      data.table::data.table(
        "Variable" = descStat$Variable,
                 MedCV,
                 "n" = descStat$n)

    descStatTable <- data.table::dcast(data = MedCV,
                          n ~ Variable ,
                          value.var = "MedCV")

    ##### ComparacionMedias ----

    if (inherits(dataToCompare, "sf")) {
      MeansComparissons <-
        lapply(colnames(sf::st_drop_geometry(dataToCompare)), function(x) {
          makeMeanComparisson(dataToCompare[,x],
                              Clasif = clusterData,
                              numCluster = numCluster,
                              descStat = descStat,
                              join = join,
                              retDMS = returnLSD)

        })
      names(MeansComparissons) <- colnames(sf::st_drop_geometry(dataToCompare))

    } else {
      MeansComparissons <-
        lapply(
          dataToCompare,
          makeMeanComparisson,
          Clasif = clusterData,
          numCluster = numCluster,
          descStat = descStat,
          retDMS = returnLSD
        )

      names(MeansComparissons) <- dataToCompare
    }



    return(list(Diferencias = MeansComparissons,
                Descriptivo = descStatTable))
  }

#' Auxiliary function for zone compare
#'
#' @noRd
varKrigDescr <-
  function(datos_predsf,
           dataToCompare,
           grid_dim) {
    is_newData <- FALSE
    # If variable is not from the same dataset, an interpolation is made assuming that
    # new data is not interpolated.
    if (inherits(dataToCompare, "sf")) {
      if (sf::st_is_longlat(dataToCompare)) {
        stop("Please provide a projected sf object", call. = FALSE)
      }
      is_newData <- TRUE
      datos_predsf <- dataToCompare
      dataToCompare <- colnames(sf::st_drop_geometry(dataToCompare))
    }
    if (missing(grid_dim))
      grid_dim <- NULL
    Krig <-
      lapply(dataToCompare, function(columna) {
        if (columna == "geometry") {
          return(NULL)
        }

        formulaKr <- stats::as.formula(paste(columna, "~ 1"))

        # Ajuste de semivariograma empírico
        semiv_emp <- gstat::variogram(formulaKr, datos_predsf)

        # Ajuste de semivariograma teóricos sin valores iniciales de los parámetros

        suppressWarnings({
           mdls <-
          gstat::fit.variogram(semiv_emp, gstat::vgm(c("Exp", "Sph", "Gau")))

           if (is_newData) {
             if (is.null(grid_dim)) {
               grid_dim <-  mean(diff(sf::st_bbox(datos_predsf),lag = 2)/100)
             }
             if (requireNamespace("stars", quietly = TRUE)) {
               grd <- sf::st_bbox(datos_predsf)
               grd <- stars::st_as_stars(grd, dx = grid_dim)
               grd <- sf::st_crop(grd, datos_predsf)

               kriging_cv <-
                 gstat::krige(formulaKr,
                              datos_predsf,
                              grd,
                              mdls,
                              nmin = 7,
                              nmax = 25)
             } else {
               message("stars package is suggested to do this comparisson")
               is_newData <- !is_newData
             }


           }
           if (!is_newData) {
             kriging_cv <-
               gstat::krige.cv(
                 formulaKr,
                 datos_predsf,
                 nfold = 10,
                 nmin = 7,
                 nmax = 25,
                 model = mdls
               )
           }

        medianKriging <- stats::median(sqrt(kriging_cv$var1.var), na.rm = T)
        })


        media <- mean(datos_predsf[[columna]], na.rm = T)
        Cv <-
          stats::sd(datos_predsf[[columna]], na.rm = T) / media * 100
        tamPunt <- sum(!is.na(datos_predsf[[columna]]))
        data.frame(
          "Variable" = as.character(columna),
          "MedianaK" = medianKriging,
          "MediaVar" = media,
          "CVVar" = Cv,
          "n" = tamPunt
        )

      })

    descStat <- do.call(rbind, Krig)
    descStat
  }

#' Auxiliary function for zone compare
#'
#' @noRd
makeMeanComparisson <-
  function(VariableEstudiada,
           Clasif,
           numCluster,
           descStat = descStat,
           alpha = 0.05,
           join = sf::st_nearest_feature,
           retDMS = FALSE) {

    if (inherits(VariableEstudiada, "sf")) {

      myData <- sf::st_join(VariableEstudiada,
                            Clasif[, numCluster],
                            join = join)
      medias <- stats::aggregate(
        myData[, colnames(myData) != numCluster],
        list(myData[[numCluster]]),
        mean, na.rm = TRUE)
      VariableEstudiada <- colnames(sf::st_drop_geometry(VariableEstudiada))
    } else {
      medias <-
        stats::aggregate(Clasif[,VariableEstudiada],
                         list(Clasif[[numCluster]]),
                         mean, na.rm = TRUE)
    }


    names(medias)[1] <- c(numCluster)
    medias <- medias[order(medias[[VariableEstudiada]], medias[[numCluster]]), ]
    ## Test ----

    comb <- utils::combn(nrow(medias), 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    pvalue <- dif
    sdtdif <- dif
    # sig <- rep(" ", nn)

    DMS <-
      descStat[which(descStat[["Variable"]] == VariableEstudiada), "MedianaK"] *
      stats::qnorm(1 - alpha / 2)# * 2 It is too much?? Why is 2 here??

    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      dif[k] <- medias[i, ][[2]] - medias[j, ][[2]]
      sdtdif[k] <- descStat[which(descStat[["Variable"]] == VariableEstudiada), "MedianaK"]
      pvalue[k] <- abs(dif[k][[1]]) <= DMS
    }

    MisLetras <-
      CalcLetras(myMeans = medias,
                 pvalue = pvalue,
                 alpha = alpha)
    if (retDMS) {
      return(list("TablaMedias" = MisLetras, "DMS" = DMS))
    }
    MisLetras
  }
#' Auxiliary function for zone compare
#'
#' @noRd
CalcLetras <- function(myMeans, pvalue, alpha = 0.05) {
  #Adapted from agricolae package
  lastC <-
    function(x) {
      y <- sub(" +$", "", x)
      p1 <- nchar(y)
      cc <- substr(y, p1, p1)
      return(cc)
    }

  Qm <- matrix(1, ncol = nrow(myMeans), nrow = nrow(myMeans))
  p <- pvalue
  k <- 0
  for (i in 1:(nrow(myMeans) - 1)) {
    for (j in (i + 1):nrow(myMeans)) {
      k <- k + 1
      Qm[i, j] <- p[k]
      Qm[j, i] <- p[k]
    }
  }



  n <- nrow(myMeans)
  z <- myMeans
  letras <- c(letters[1:26],LETTERS[1:26],1:9,c(".","+","-","*","/","#","$",
                                                "%","&","^","[","]",":","@",
                                                ";","_","?","!","=","#",rep(" ",2000)))
  w <- z[order(z[[2]], decreasing = FALSE),]
  M <- rep("", n)
  k <- 1
  # k1 <- 0
  j <- 1
  i <- 1
  cambio <- n
  cambio1 <- 0
  chequeo = 0
  M[1] <- letras[k]
  q <- as.numeric(rownames(w)) #Check
  while (j < n) {
    chequeo <- chequeo + 1
    if (chequeo > n)
      break
    for (i in j:n) {
      stest <- Qm[q[i], q[j]] > alpha
      if (stest) {
        if (lastC(M[i]) != letras[k]) {
          M[i] <- paste(M[i], letras[k], sep = "")
        }
      }
      else {
        k <- k + 1
        cambio <- i
        cambio1 <- 0
        ja <- j
        for (jj in cambio:n)
          M[jj] <- paste(M[jj], "", sep = "")
        M[cambio] <- paste(M[cambio], letras[k], sep = "")
        for (v in ja:cambio) {
          if (Qm[q[v], q[cambio]] <= alpha) {
            j <- j + 1
            cambio1 <- 1
          }
          else
            break
        }
        break
      }
    }
    if (cambio1 == 0)
      j <- j + 1
  }

  w <- data.frame(w, stat = M)
  if (k > 81
  )
    cat(
      "\n",
      k,
      "groups are estimated.The number of groups exceeded the maximum of 81 labels. Change to group=FALSE.\n"
    )
  invisible(w)
}

