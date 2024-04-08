#### Polynomial model ####
#' plot() overloading for polynomial model
#'
#' @param x the polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.polynomial_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
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
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
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
  CEX_SCALER <- 4 # arbitrary number for better visual

  df_ <- transform(x$df$t, x$df$spos)
  names(df_)[names(df_) == "t"] <- "exposure"

  with(c(x$df, df_), {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      exposure,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="exposure", ylab="seroprevalence",
      xlim=c(0, max(exposure)), ylim=c(0,1)
    )
    lines(t, x$sp, lwd=2)
    lines(t, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=10), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#### Fractional polynomial model ####
#' plot() overloading for fractional polynomial model
#'
#' @param x the fractional polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.fp_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age[c(-1,-length(age))], x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
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
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
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
#'
#' @examples
#' plot(model, y1 = "Parvo B19", y2 = "VZV", plot_type = "sp")
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
#' @import locfit
#' @import graphics
#'
#' @export
plot_gcv <- function(age, pos, tot, nn_seq, h_seq, kern="tcub", deg=2) {
  par(mfrow=c(1,2),lwd=2,las=1,cex.axis=1,cex.lab=1.1,mgp=c(2,0.5, 0),mar=c(3.1,3.1,3.1,3))

  y <- pos/tot
  res = cbind(nn_seq, summary(gcvplot(y~age, family="binomial", alpha=nn_seq)))
  plot(res[,1],res[,3],type="n",xlab="nn (% Neighbors)",ylab=" ")
  lines(res[,1],res[,3])
  mtext(side=2, "GCV", las=3, line=2.4, cex=0.9)

  h_seq_ <- cbind(rep(0, length(h_seq)), h_seq)
  res=cbind(h_seq_[,2],summary(gcvplot(y~age,family="binomial",alpha=h_seq_)))
  plot(res[,1],res[,3],type="n",xlab="h (Bandwidth)",ylab=" ")
  lines(res[,1],res[,3])
  mtext(side=2, "GCV", las=3, line=3, cex=0.9)
  # reset the plotting parameter back to its default state
  par(mfrow = c(1, 1))
}
