est_foi <- function(t, sp)
{
  dsp <- diff(sp)/diff(t)
  foi <- approx(
    (t[-1]+t[-length(t)])/2,
    dsp,
    t[c(-1,-length(t))]
  )$y/(1-sp[c(-1,-length(t))])

  foi
}

est.foi <- function(t, sp)
{
  t_ <- sort(unique(t))
  sp_ <- (sp[order(t)])[duplicated(sort(t))==F]
  dsp <- diff(sp_)/diff(t_)
  foi <- approx(
    (t_[-1]+t_[-length(t_)])/2,
    dsp,
    t_[c(-1,-length(t_))]
  )$y/(1-sp_[c(-1,-length(t_))])

  list(x=t_, y=foi)
}

transform <- function(t, spos, round.by=NULL) {
  if (!missing(round.by) & !is.null(round.by)) {
    t <- round(t, digits=round.by)
  }
  df <- data.frame(t, spos)
  df_agg <- df %>%
    group_by(t) %>%
    summarize(
      pos = sum(spos),
      tot = n()
    )
  df_agg
}
