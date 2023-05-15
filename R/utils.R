
rename_col <- function(df, old, new)
{
  names(df)[names(df) == old] <- new
  df
}
