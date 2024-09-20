#' Helper to adjust styling of a plot
#'
#' @param sero - color for seroprevalence line
#' @param ci - color for confidence interval
#' @param foi - color for force of infection line
#' @param sero_line - linetype for seroprevalence line
#' @param foi_line - linetype for force of infection line
#'
#' @export
set_plot_style <- function(sero = "blueviolet", ci = "royalblue1", foi = "#fc0328", sero_line = "solid", foi_line = "dashed"){
    list(
      scale_colour_manual(
        values = c("sero" = sero, "foi" = foi)
      ),
      scale_linetype_manual(
        values = c("sero" = sero_line, "foi" = foi_line)
      ),
      scale_fill_manual(
        values = c("ci" =ci)
      ),
      labs(linetype = "Line", colour = "Line", fill="Fill color")
    )
}

# ==== Helper function to compute plot data ====
plot_data <- function(x){
  if(x$datatype == "linelisting"){
    # transform data before plotting
    df_ <- transform_data(x$df$age, x$df$pos)
    age <- df_$t
    pos <- df_$pos
    tot <- df_$tot
    # use pre-aggregated age for FOI
  }else if (x$datatype == "aggregated"){
    age <- x$df$age
    pos <- x$df$pos
    tot <- x$df$tot
  }

  return(data.frame(
    "age" = age, "pos" = pos, "tot" = tot
  ))
}

#=== Helper function for plotting =====
plot_util <- function(age, pos, tot, sero, foi, cex = 20){
  plot <- ggplot(data = data.frame(age, pos, tot), aes(x = age, y = pos/tot)) +
    geom_point(size = cex*(pos)/max(tot), shape = 1)+
    labs(y="Seroprevalence", x="Age")+
    coord_cartesian(xlim=c(0,max(age)), ylim=c(0, 1)) +
    scale_y_continuous(
      name = "Seroprevalence",
      sec.axis = sec_axis(~.*1, name = " Force of infection")
    ) + set_plot_style()

  # === Add seroprevalence layer
  if (is(sero, "data.frame")){
    if("ymax" %in% colnames(sero)){
      plot <- plot + geom_smooth(aes(x = x, y = y, ymin = ymin, ymax = ymax, col = "sero", linetype="sero",
                                     fill = "ci"), data=sero,
                                 stat="identity",lwd=0.5)
    }else{
      # --- Handle cases where CI for seroprevalence is not computable & length of age for foi differs from provided age vector
      plot <- plot + geom_line(aes(x = x, y = y, col = "sero", linetype="sero"), data = sero,
                               stat="identity",lwd=0.5)
    }
  }else{
    # --- Simply plot seroprevalence line if CI cannot be computed
    plot <- plot + geom_line(aes(x = age, y = sero, col = "sero", linetype="sero"),
                             stat="identity",lwd=0.5)
  }

  # === Add foi layers
  if (is(foi, "data.frame")){
    if ("ymax" %in% colnames(foi)){
      # --- Handle cases where FOI is a data.frame (with CI)
      # plot <- plot + geom_smooth(aes_auto(foi, col = "foi", linetype="foi", fill="ci"), data=foi,
      #                          stat="identity",lwd=0.5)
      plot <- plot + geom_smooth(aes(x = x, y=y, ymin =ymin, ymax = ymax, col = "foi", linetype="foi", fill="ci"), data=foi,
                               stat="identity",lwd=0.5)
    }else{
      # --- Handle cases where CI for FOI is not computable & length of age for foi differs from provided age vector
      plot <- plot + geom_line(aes(x = x, y=y, col = "foi", linetype="foi"), data=foi,
                               stat="identity",lwd=0.5)
    }
  }else if (length(age) != length(foi)){
    # --- handle some cases when length of age differs from length of foi
    age <- age[c(-1,-length(age))]
    foi <- data.frame(x = age, y = foi)
    plot <- plot + geom_line(aes(x = x, y=y, col = "foi", linetype="foi"), data = foi,
                             lwd = 0.5)
  }else{
    # --- Simply plot foi
    plot <- plot + geom_line(aes(x = age, y = foi, col = "foi", linetype="foi"),
                lwd = 0.5)
  }
  plot
}

# Force using the generic plot function to override
plot <- function(x, ...) {
  UseMethod("plot")
}

#### Polynomial model ####
#' plot() overloading for polynomial model
#'
#' @param x the polynomial model object
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.polynomial_model <- function(x, ...) {
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex

  out.DF <- compute_ci(x)

  if(x$datatype == "linelisting"){
    # use pre-aggregated age for FOI
    foi <- data.frame(x = x$df$age, y = as.numeric(x$foi))
  }else if (x$datatype == "aggregated"){
    foi <- as.numeric(x$foi)
  }

  to_plot <- plot_data(x)

  with(x$df, {
    plot_util(age = to_plot$age, pos = to_plot$pos, tot = to_plot$tot, sero = out.DF, foi = foi, cex = cex)
  })

}

#### Non-linear ####

#### Farrington model ####
#' plot() overloading for Farrington model
#'
#' @param x the Farrington model object.
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.farrington_model <- function(x,...) {
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex
  # out.DF <- compute_ci(x)

  to_plot <- plot_data(x)


  if(x$datatype == "linelisting"){
    # use pre-aggregated age for FOI & sero
    foi <- data.frame(x = x$df$age, y = x$foi)
    sero <- data.frame(x = x$df$age, y = x$sp)
  }else if (x$datatype == "aggregated"){
    foi <- x$foi
    sero <- x$sp
  }

  with(x$df, {
    plot_util(age = to_plot$age, pos = to_plot$pos, tot = to_plot$tot, sero = sero, foi = foi, cex = cex)
  })
}

#### Weibull model ####
#' plot() overloading for Weibull model
#'
#' @param x the Weibull model object.
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.weibull_model <- function(x, ...) {
  # df_ <- transform_data(x$df$t, x$df$spos)
  # names(df_)[names(df_) == "t"] <- "exposure"
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex

  out.DF <- compute_ci.weibull_model(x)

  to_plot <- plot_data(x)

  if(x$datatype == "linelisting"){
    # use pre-aggregated age for FOI & sero
    foi <- data.frame(x = x$df$age, y = x$foi)
  }else if (x$datatype == "aggregated"){
    foi <- x$foi
  }

  plot_util(age = to_plot$age, pos = to_plot$pos, tot = to_plot$tot, sero = out.DF, foi = data.frame(x = x$df$age, y = x$foi), cex = cex)
}

#### Fractional polynomial model ####
#' plot() overloading for fractional polynomial model
#'
#' @param x the fractional polynomial model object.
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.fp_model <- function(x,...) {
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex

  out.DF <- compute_ci.fp_model(x)
  to_plot <- plot_data(x)

  with(x$df, {
    plot_util(age = to_plot$age, pos = to_plot$pos, tot = to_plot$tot, sero = out.DF, foi = x$foi, cex = cex)
  })
}

#### Non-parametric ####

#### Local polynomial model ####
#' plot() overloading for local polynomial model
#'
#' @param x the local polynomial model object.
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.lp_model <- function(x, ...) {
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex

  out.DF <- compute_ci.lp_model(x)
  to_plot <- plot_data(x)

  if(x$datatype == "linelisting"){
    # use pre-aggregated age for FOI
    foi <- data.frame(x = x$df$age, y = as.numeric(x$foi))
  }else if (x$datatype == "aggregated"){
    foi <- x$foi
  }

  with(x$df, {
    plot_util(age = to_plot$age, pos = to_plot$pos, tot = to_plot$tot, sero = out.DF, foi = foi, cex=cex)
  })
}

#### Hierarchical Bayesian model ####
#' plot() overloading for hierarchical_bayesian_model
#'
#' @param x hierarchical_bayesian_model object created by serosv.
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.hierarchical_bayesian_model <- function(x,  ...){
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex

  with(x$df, {
    plot_util(age = age, pos = pos, tot = tot, sero = x$sp, foi = x$foi, cex=cex)
  })
}


#### Penalized splines ####
#' plot() overloading for penalized spline
#'
#' @param x the penalized_spline_model object
#' @param ... arbitrary params.
#' @import ggplot2
#' @importFrom methods is
#'
#' @export
plot.penalized_spline_model <- function(x, ...){
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex
  ci <- compute_ci.penalized_spline_model(x)

  out.DF <- ci[[1]]
  out.FOI <- ci[[2]]

  to_plot <- plot_data(x)

  with(x$df, {
    plot_util(age = to_plot$age, pos = to_plot$pos, tot = to_plot$tot, sero = out.DF, foi = out.FOI, cex=cex)
  })

}


#### Mixture model ####
#' plot() overloading for mixture model
#'
#' @param x the mixture_model
#' @param ... arbitrary params.
#' @import ggplot2
#'
#' @export
plot.mixture_model <- function(x, ...){
  ci_layer <- function(x, y, xmin, xmax, fill="royalblue1"){
    idx_lim <- which(x>xmin & x<xmax)
    x <- x[idx_lim]
    y <- y[idx_lim]

    # add values to make sure shape is filled at y=0
    x <- c(x[1], x, x[length(x)])
    y <- c(0, y, 0)
    geom_polygon(aes(x=x, y=y), fill = fill, alpha=0.3)
  }

  ci <- compute_ci.mixture_model(x)

  with(x$df, {
    # ---- code to compute probability density from mixgroup output
    ntot <- sum(count)
    m <- length(count)
    iwid <- antibody_level[2:(m - 1)] - antibody_level[1:(m - 2)]
    iwid <- c(2 * iwid[1], iwid, 2 * iwid[m - 2])
    idens <- (count/iwid)/ntot
    # ---- end code from mixdist


    # --- Get plotting data from parameter to have prettier plot =))
    x_coord <- seq(min(antibody_level), max(antibody_level[is.finite(antibody_level)]), length.out=100)
    fitted_susceptible <- dnorm(x_coord, x$info$parameters[1, ]$mu, x$info$parameters[1, ]$sigma)*x$info$parameters[1, ]$pi
    fitted_infected <- dnorm(x_coord, x$info$parameters[2, ]$mu, x$info$parameters[2, ]$sigma)*x$info$parameters[2, ]$pi


    ggplot() +
      ci_layer(x = x_coord, y = fitted_susceptible, xmin = ci$susceptible$lower_bound,
               xmax = ci$susceptible$upper_bound, fill="forestgreen") +
      ci_layer(x = x_coord, y = fitted_infected, xmin = ci$infected$lower_bound,
               xmax = ci$infected$upper_bound, fill="blueviolet") +
      geom_step(aes(x=antibody_level-iwid, y = idens)) +
      geom_line(aes(x=x_coord, y=fitted_susceptible, col = "susceptible")) +
      geom_line(aes(x=x_coord, y=fitted_infected, col = "infected")) +
      annotate("segment", x=x$info$parameters[1, ]$mu, y=0.01, yend=-0.01, colour="#fc0328",)+
      annotate("segment", x=x$info$parameters[2, ]$mu, y=0.01, yend=-0.01, colour="#fc0328",)+
      labs(x="Log(Antibody level+1)", y="Probability Density") +
      scale_color_manual(
        values = c("susceptible"= "forestgreen", "infected"="blueviolet")
      )+
      coord_cartesian(xlim = c(0, max(antibody_level[is.finite(antibody_level)])) ,ylim = c(0, max(idens)))
  })

}

#### Estimate sero from mixture model #####
#' plot() overloading for result of estimate_from_mixture
#'
#' @param x the mixture_model
#' @param ... arbitrary params.
#' @import ggplot2
#'
#' @export
plot.estimate_from_mixture <- function(x, ... ){
  cex <- if (is.null(list(...)[["cex"]])) 20 else list(...)$cex
  age <- x$df$age

  returned_plot <- ggplot()

  if(!is.null(x$df$threshold_status)){
    aggregated <- transform_data(round(x$df$age), x$df$threshold_status)

    returned_plot <-  returned_plot +
      geom_point(aes( x = t, y = pos/tot, size = cex*(pos)/max(tot) ), data = aggregated,
                 shape = 1, show.legend = FALSE)
  }

  returned_plot <- returned_plot +
    geom_line(aes(x = age, y = x$sp, col = "sero", linetype = "sero")) +
    geom_line(aes(x = foi_x, y = foi, col = "foi", linetype = "foi"), data=x$foi)

  returned_plot + set_plot_style() + labs(x = "Age", y="Seroprevalence")
}


#### GCV values ####
#' Plotting GCV values with respect to different nn-s and h-s parameters.
#'
#' Refers to section 7.2.
#'
#' @param age the age vector.
#' @param pos the pos vector.
#' @param tot the tot vector.#'
#' @param nn_seq Nearest neighbor sequence.
#' @param h_seq Smoothing parameter sequence.
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' df <- mumps_uk_1986_1987
#' plot_gcv(
#'   df$age, df$pos, df$tot,
#'   nn_seq = seq(0.2, 0.8, by=0.1),
#'   h_seq = seq(5, 25, by=1)
#' )
#'
#' @import locfit patchwork ggplot2
#' @import graphics
#'
#' @export
plot_gcv <- function(age, pos, tot, nn_seq, h_seq, kern="tcub", deg=2) {
  y <- pos/tot

  # --- Plot for nn seq
  res <-  cbind(nn_seq, summary(gcvplot(y~age, family="binomial", alpha=nn_seq)))
  nn_plot <- ggplot() +
    geom_line(aes(x = res[,1], y = res[,3]), col = "royalblue") +
    labs(x="nn (% Neighbors)", y="GCV")


  # --- Plot for h seq
  h_seq_ <- cbind(rep(0, length(h_seq)), h_seq)
  h_res <- cbind(h_seq_[,2],summary(gcvplot(y~age,family="binomial",alpha=h_seq_)))
  h_plot <- ggplot() +
    geom_line(aes(x = h_res[,1], y = h_res[,3]), col = "royalblue") +
    labs(x="h (Bandwidth)", y="GCV")

  # --- Combine 2 plots
  nn_plot + h_plot + plot_layout(ncol=2)
}
