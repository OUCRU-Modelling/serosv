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

#--- Helper function for plotting
plot_util <- function(age, pos, tot, sero, foi){
  plot <- ggplot(data = data.frame(age, pos, tot), aes(x = age, y = pos/tot)) +
    geom_point(size = 20*(pos)/max(tot), shape = 1)+
    labs(y="Seroprevalence", x="Age")+
    coord_cartesian(xlim=c(0,max(age)), ylim=c(0, 1)) +
    scale_y_continuous(
      name = "Seroprevalence",
      sec.axis = sec_axis(~.*1, name = " Force of infection")
    ) + set_plot_style()

  # === Add seroprevalence layer
  if (class(sero) == "data.frame"){
    # --- Sero is dataframe when confidence interval is computed
    plot <- plot + geom_smooth(aes_auto(sero, col = "sero", linetype="sero",
                                        fill = "ci"), data=sero,
                               stat="identity",lwd=0.5)
  }else{
    # --- Simply plot seroprevalence line if CI cannot be computed
    plot <- plot + geom_line(aes(x = age, y = sero, col = "sero", linetype="sero"),
                             stat="identity",lwd=0.5)
  }

  # === Add foi layers
  if (class(foi) == "data.frame"){
    if ("ymax" %in% colnames(foi)){
      # --- Handle cases like that of weibull model
      plot <- plot + geom_smooth(aes_auto(foi, col = "foi", linetype="foi", fill="ci"), data=foi,
                               stat="identity",lwd=0.5)
    }else{
      # --- Handle cases like that of weibull model
      plot <- plot + geom_line(aes_auto(foi, col = "foi", linetype="foi"), data=foi,
                               stat="identity",lwd=0.5)
    }
  }else if (length(age) != length(foi)){
    # --- handle some cases when length of age differs from length of foi
    age <- age[c(-1,-length(age))]
    foi <- data.frame(x = age, y = foi)
    plot <- plot + geom_line(aes_auto(foi, col = "foi", linetype="foi"), data = foi,
                             lwd = 0.5)
  }else{
    # --- Simply plot foi
    plot <- plot + geom_line(aes(x = age, y = foi, col = "foi", linetype="foi"),
                lwd = 0.5)
  }
  plot
}

#### Polynomial model ####
#' plot() overloading for polynomial model
#'
#' @param x the polynomial model object
#' @param ... arbitrary params.
#' @import ggplot2
#'
#' @export
plot.polynomial_model <- function(x, ...) {
  out.DF <- compute_ci(x)

  with(x$df, {
    plot_util(age = x$df$Age, pos = x$df$Pos, tot = x$df$Tot, sero = out.DF, foi = as.numeric(x$foi))
  })

}

#### Non-linear ####

#### Farrington model ####
#' plot() overloading for Farrington model
#'
#' @param x the Farrington model object.
#' @param ... arbitrary params.
#'
#' @export
plot.farrington_model <- function(x, ...) {
  # out.DF <- compute_ci(x)

  with(x$df, {
    plot_util(age = age, pos = pos, tot = tot, sero = x$sp, foi = x$foi)
  })
}

#### Weibull model ####
#' plot() overloading for Weibull model
#'
#' @param x the Weibull model object.
#' @param ... arbitrary params.
#'
#' @export
plot.weibull_model <- function(x, ...) {
  df_ <- transform(x$df$t, x$df$spos)
  names(df_)[names(df_) == "t"] <- "exposure"

  out.DF <- compute_ci.weibull_model(x)
  plot_util(age = df_$exposure, pos = df_$pos, tot = df_$tot, sero = out.DF, foi = data.frame(x = x$df$t, y = x$foi))
}

#### Fractional polynomial model ####
#' plot() overloading for fractional polynomial model
#'
#' @param x the fractional polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.fp_model <- function(x, ...) {
  out.DF <- compute_ci.fp_model(x)

  with(x$df, {
    plot_util(age = age, pos = pos, tot = tot, sero = out.DF, foi = x$foi)
  })
}

#### Non-parametric ####

#### Local polynomial model ####
#' plot() overloading for local polynomial model
#'
#' @param x the local polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.lp_model <- function(x, ...) {
  out.DF <- compute_ci.lp_model(x)

  with(x$df, {
    plot_util(age = age, pos = pos, tot = tot, sero = out.DF, foi = x$foi)
  })
}

#### Penalized splines ####

#' plot() overloading for penalized spline fitted with GLMM
#'
#' @param x the glmm_ps_model object
#' @param ... arbitrary params.
#'
#' @export
plot.glmm_ps_model <- function(x, ...){
  ci <- compute_ci.glmm_ps_model(x)

  out.DF <- ci[[1]]
  out.FOI <- ci[[2]]

  with(x$df, {
    plot_util(age = age, pos = pos, tot = tot, sero = out.DF, foi = out.FOI)
  })

}

#### Mixture model ####
#' plot() overloading for mixture model
#'
#' @param x the mixture_model
#' @param ... arbitrary params.
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


#### Bivariate Dale model ####
#' plot() overloading for bivariate_dale_model
#'
#' @param x the bivariate_dale_model object
#' @param ... arbitrary params including y1, y2 (label for predictors), plot_type ("ci" or "sp")
#'
#' @import patchwork ggplot2
#' @importFrom stringr str_interp
#'
#' @return patchwork object
#' @export
plot.bivariate_dale_model <- function(x, ...) {
  y1 <- if (is.null(list(...)[["y1"]])) "Y1" else list(...)$y1
  y2 <- if (is.null(list(...)[["y2"]])) "Y2" else list(...)$y2
  plot_type <- if (is.null(list(...)[["plot_type"]])) "ci" else list(...)$plot_type

  ci_plot <- function(age, ci, title="Marginal seroprevalence", ylabel="Seroprevalence"){
    ci_df <- data.frame(val = ci[1,], lower_bound = ci[2,], upper_bound = ci[3,])
    ggplot()+
      ggtitle(title) +
      # geom_point(size = 7*(pos)/max(tot), shape = 1)+
      labs(y="Seroprevalence", x="Age")+
      coord_cartesian(xlim=c(0, max(age)), ylim=c(0, max(ci_df$val)))+
      geom_smooth(aes(x = age, y = ci_df$val, ymin = ci_df$lower_bound, ymax = ci_df$upper_bound), stat="identity",
                  col = "blueviolet",fill = "royalblue1",
                  lwd = .5) +
      scale_y_continuous(
        name = ylabel
      )
  }

  cond_plot <- function(age, sp, foi, title = "Conditional and Marginal Seroprevalence and FOI"){
    sp <- data.frame(sp)
    foi <- data.frame(foi)
    ggplot(mapping = aes(x = age)) +
      ggtitle(title) +
      labs(y="Seroprevalence", x="Age")+
      coord_cartesian(xlim=c(0, max(age)), ylim=c(0, 1))+
      scale_colour_manual(
        values = c("Marginal"="royalblue1",
                   "Conditioning on occurence"="#1ce8bc",
                   "Conditioning on non-occurence"="#fc0328")
      )+
      scale_linetype_manual(
        values = c("Seroprevalence"="solid", "FOI"="dashed")
      )+
      geom_line(aes(y = sp[,1], col="Marginal", linetype = "Seroprevalence")) +
      geom_line(aes(y = sp[,2], col ="Conditioning on occurence", linetype = "Seroprevalence")) +
      geom_line(aes(y = sp[,3], col ="Conditioning on non-occurence", linetype = "Seroprevalence")) +
      geom_line(aes(x = age[c(-1, -length(age))], y = foi[,1], col ="Marginal", linetype = "FOI")) +
      geom_line(aes(x = age[c(-1, -length(age))], y = foi[,2], col ="Conditioning on occurence", linetype = "FOI")) +
      geom_line(aes(x = age[c(-1, -length(age))], y = foi[,3], col ="Conditioning on non-occurence", linetype = "FOI")) +
      scale_y_continuous(
        name = "Seroprevalence",
        sec.axis = sec_axis(~.*1, name = " Force of infection")
      )

  }

  data_points_layer <- function(age, pos, tot){
    geom_point(aes(x = age, y = pos/ tot), size = 7*(pos)/max(tot), shape = 1)
  }

  plot_list <- list()

  # --- Plot marginal seroprevalence and CI
  plot_list[["y1m"]] <- ci_plot(age = x$df$age, ci = x$ci$y1m, title = str_interp("Marginal seroprevalence of ${y1}")) +
    data_points_layer(age = x$df$age, pos=rowSums(x$df$y[,c("PN", "PP")]),tot=rowSums(x$df$y))

  plot_list[["y2m"]] <- ci_plot(age = x$df$age, ci = x$ci$y2m, title = str_interp("Marginal seroprevalence of ${y2}")) +
    data_points_layer(age = x$df$age, pos=rowSums(x$df$y[,c("NP", "PP")]),tot=rowSums(x$df$y))

  plot_list[["or"]] <- ci_plot(age = x$df$age, ci = x$ci$or, title = "Odd ratio", ylabel = "Odd ratio") +
    geom_hline(aes(yintercept = 1), col = "#fc0328", lwd = 0.5, linetype="dashed")

  # --- Plot conditional vs marginal sp and FOI
  plot_list[["y1condy2"]] <- cond_plot(age = x$df$age,
                                       sp = x$sp[c("y1m", "y1condy2pos", "y1condy2neg")],
                                       foi = x$foi[c("y1m", "y1condy2pos", "y1condy2neg")],
                                       title = str_interp("SP and FOI for ${y1} given ${y2}"))

  plot_list[["y2condy1"]] <- cond_plot(age = x$df$age,
                                       sp = x$sp[c("y2m", "y2condy1pos", "y2condy1neg")],
                                       foi = x$foi[c("y2m", "y2condy1pos", "y2condy1neg")],
                                       title = str_interp("SP and FOI for ${y2} given ${y1}"))


  # --- Plot joint sp and FOI
  plot_list[["joint"]] <- ci_plot(age = x$df$age, ci = x$ci$pi11, title = "Joint Seroprevalence and FOI") +
    geom_line(aes(x = x$df$age[c(-1, -length(x$df$age))], y = x$foi$joint),
              col ="blueviolet", linetype="dashed") +
    data_points_layer(age = x$df$age, pos = x$df$y$PP, tot = rowSums(x$df$y))

  # --- Return plots depending on options
  if (plot_type == "ci"){
    return(
      wrap_plots(plot_list[c("y1m", "y2m", "or")]) + guide_area() + plot_layout(guides = "collect")
      & theme(plot.title = element_text(size = "11"))
    )
  }else if (plot_type == "sp"){
    return(
      wrap_plots(plot_list[c("y1condy2", "y2condy1", "joint")]) + guide_area() + plot_layout(guides = "collect")
      & theme(plot.title = element_text(size = "11"))
    )
  }else{
    return(plot_list)
  }

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
