respDiv <- function(model, newdata, gradient = c("sst"), unit = "bird1", resp = c("weight"), ndraws = 1000, type = c("link","response"), raw = FALSE, width = 0.95){

  epreds <- tidybayes::add_epred_draws(model, newdata = newdata,
                                       re_formula = as.formula(paste0("~ (1|",unit,")")),
                                       resp = resp,
                                       ndraws = ndraws) |>
    dplyr::ungroup() |>
    dplyr::select(-c(.chain,.iteration))
  
  epreds <- if(raw){
   epreds
    }else{
   epreds  |>
        #dplyr::reframe(.epred = mean(.epred), .by = c(!!as.symbol(gradient), !!as.symbol(unit), .draw)) |>
        dplyr::group_by(!!as.symbol(unit), !!as.symbol(gradient)) |>
        ggdist::median_qi(.exclude = c(".row",".chain", ".iteration", ".draw"),
                          .width = width) |> #estimate 80% quartiles for the derivatives
        dplyr::ungroup()
  }
  
  m <- marginaleffects::slopes(model,
                               newdata = newdata,
                               re_formula = as.formula(paste0("~ (1|",unit,")")),
                               resp = resp,
                               ndraws = ndraws,
                               variables = gradient,
                               by = c(gradient,unit),
                               type = type)

  dr <- marginaleffects::get_draws(m) |>
    dplyr::select(drawid, draw,dplyr::matches(gradient), dplyr::matches(unit)) |>
    tidyr::pivot_wider(id_cols = c(dplyr::matches(gradient),drawid),names_from = !!as.symbol(unit), values_from = "draw")

  dr[[".dissimilarity"]] <- apply(
    X = dr[, grepl(paste(unique(newdata[[unit]]),collapse = "|"),colnames(dr))],
    MARGIN = 1,
    FUN = function(x) {
      #https://doi.org/10.1111/2041-210X.14087
      if (length(unique(x)) == 1) {
        dis = 0
      } else{
        d <- dist(x, diag = T) # euclidean distance matrix
        Z <- exp(as.matrix(-d)) # similarity matrix
        nspecies <- length(x)
        p = matrix(1 / nspecies, nspecies) # relative abundance matrix. Needs to be changed if evenness is of interest
        lenq = 1 # initialises hill number. Need to fix if we want any evenness indices
        qq <- seq(length = lenq,
                  from = 0,
                  by = .11)

        # Initialise the Zp matrix to zero
        Zp = matrix(0, nspecies, 1)

        # Compute Zp
        for (i in 1:nspecies) {
          for (j in 1:nspecies) {
            Zp[i, 1] <- Zp[i, 1] + Z[i, j] * p[j, 1]
          }
        }

        # Initialise the Diversity matrix to zero
        Dqz = matrix(0, lenq , 1)

        for (iq in 1:lenq)  {
          q <- qq[iq]

          for (zpi in 1:length(Zp[, 1])) {
            if (Zp[zpi, 1] > 0)
              (
                Dqz[iq, 1] <- Dqz[iq, 1] + p[zpi, 1] * (Zp[zpi, 1])^(q - 1))
          }

          Dqz[iq, 1] <- Dqz[iq, 1]^(1 / (1 - q))

        }
        dis <- Dqz[iq, 1]
      }
      dis

    }
  )

  dr[[".divergence"]] <- apply(
    X = dr[, grepl(paste(unique(newdata[[unit]]),collapse = "|"),colnames(dr))],
    MARGIN = 1,
    FUN = function(x) {
      #https://doi.org/10.1111/2041-210X.14087
      if (length(unique(x)) == 1) {
        div = 0
      } else{
        div1 <- max(x) - min(x)
        div2 <- abs(abs(max(x)) - abs(min(x)))
        div <- (div1 - div2) / div1
      }

      div

    }
  )
  
  dr[[".imbalance"]] <- apply(
    X = dr[, grepl(paste(unique(newdata[[unit]]),collapse = "|"),colnames(dr))],
    MARGIN = 1,
    FUN = function(x) {
      #https://doi.org/10.1111/ele.70224
      S = length(x)
      abs((1/S)*sum(x))
      #abs(mean(x))
    }
  )
  
  dissim <- if(raw){
    dr[,c(gradient,"drawid",".dissimilarity")]
  }else{
    dr[,c(gradient,"drawid",".dissimilarity")] |>
      dplyr::group_by(!!as.symbol(gradient)) |>
      ggdist::median_qi(
        .exclude = c(".chain", ".iteration", "drawid", ".row"),
        .width = width
      ) |> #estimate 80% quartiles for the derivatives
      dplyr::ungroup()
  }
  
  diverg <- if(raw){
    dr[,c(gradient,"drawid",".divergence")]
  }else{
    dr[,c(gradient,"drawid",".divergence")] |>
      dplyr::group_by(!!as.symbol(gradient)) |>
      ggdist::median_qi(
        .exclude = c(".chain", ".iteration", "drawid", ".row"),
        .width = width
      ) |> #estimate 80% quartiles for the derivatives
      dplyr::ungroup()
  }
  
  imbalance <- if(raw){
    dr[,c(gradient,"drawid",".imbalance")]
  }else{
    dr[,c(gradient,"drawid",".imbalance")] |>
      dplyr::group_by(!!as.symbol(gradient)) |>
      ggdist::median_qi(
        .exclude = c(".chain", ".iteration", "drawid", ".row"),
        .width = width
      ) |> #estimate 80% quartiles for the derivatives
      dplyr::ungroup()
  }

  m <- if(raw){
    marginaleffects::get_draws(m) |>
      dplyr::select(drawid, draw,dplyr::matches(gradient), dplyr::matches(unit)) |>
      dplyr::rename(estimate = draw)
  }else{
    m |>
      as.data.frame() |>
      dplyr::select(!!as.symbol(gradient), !!as.symbol(unit), estimate, conf.low, conf.high) |>
      dplyr::rename(.lower = conf.low, .upper = conf.high)
  }
  

  list("derivs" = m,
       "epreds" = epreds,
       "dissim" = dissim,
       "diverg" = diverg,
       "imbalance" = imbalance
  )

}


