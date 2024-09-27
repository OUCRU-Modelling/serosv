# Serological survey data of infectious diseases
#
# These samples are cross-sectional surveys providing information about
# whether or not the individual has been infected before that time point.
# The data is taken from "Modeling Infectious Disease Parameters Based on
# Serological and Social Contact Data" - N. Hens et al. 2013
# data are normalized and cleaned for the purpose of this package


#' Hepatitis A serological data from Belgium in 1993 and 1994 (aggregated)
#'
#' A study of the prevalence of HAV antibodies conducted in the Flemish
#' Community of Belgium in 1993 and early 1994
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#' }
#'
#' @examples
#' # Reproduce Fig 4.1 (upper left panel), p. 63
#' age <- hav_be_1993_1994$age
#' pos <- hav_be_1993_1994$pos
#' tot <- hav_be_1993_1994$tot
#' plot(
#'   age, pos / tot,
#'   pty = "s", cex = 0.06 * tot, pch = 16, xlab = "age",
#'   ylab = "seroprevalence", xlim = c(0, 86), ylim = c(0, 1)
#' )
#'
#' @source Beutels, M., Van Damme, P., Aelvoet, W. et al. Prevalence of
#' hepatitis A, B and C in the Flemish population. Eur J Epidemiol 13,
#' 275-280 (1997). \doi{doi:10.1023/A:1007393405966}
"hav_be_1993_1994"

#' Hepatitis A serological data from Belgium in 2002 (line listing)
#'
#' A subset of the serological dataset of Varicella-Zoster Virus (VZV) and
#' Parvovirus B19 in Belgium where only individuals living in Flanders were
#' selected
#'
#' @format A data frame with 2 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#' }
#'
#' @examples
#' # Reproduce Fig 4.1 (upper right panel), p. 63
#' library(dplyr)
#' df <- hav_be_2002 %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(
#'   df$age, df$pos / df$tot,
#'   pty = "s", cex = 0.06 * df$tot, pch = 16, xlab = "age",
#'   ylab = "seroprevalence", xlim = c(0, 86), ylim = c(0, 1)
#' )
#'
#' @source Thiry, N., Beutels, P., Shkedy, Z. et al. The seroepidemiology
#' of primary varicella-zoster virus infection in Flanders (Belgium).
#' Eur J Pediatr 161, 588-593 (2002).
#' \doi{doi:10.1007/s00431-002-1053-2}
"hav_be_2002"

#' Hepatitis A serological data from Bulgaria in 1964 (aggregated)
#'
#' A cross-sectional survey conducted in 1964 in Bulgaria. Samples were
#' collected from schoolchildren and blood donors.
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#' }
#'
#' @examples
#' # Reproduce Fig 4.1 (lower panel), p. 63
#' age <- hav_bg_1964$age
#' pos <- hav_bg_1964$pos
#' tot <- hav_bg_1964$tot
#' plot(
#'     age, pos / tot,
#'     pty = "s", cex = 0.08 * tot, pch = 16, xlab = "age",
#'     ylab = "seroprevalence", xlim = c(0, 86), ylim = c(0, 1)
#' )
#'
#' @source Keiding, Niels. "Age-Specific Incidence and Prevalence: A
#' Statistical Perspective." Journal of the Royal Statistical Society.
#' Series A (Statistics in Society) 154, no. 3 (1991): 371-412.
#' \doi{doi:10.2307/2983150}
"hav_bg_1964"

#' Hepatitis B serological data from Russia in 1999 (aggregated)
#'
#' A seroprevalence study conducted in St. Petersburg
#' (more information in the book)
#'
#' @format A data frame with 4 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#'  \item{gender}{Gender of cohort (unsure what 1 and 2 means)}
#' }
#'
#' @examples
#' # Reproduce Fig 4.2, p. 65
#' library(dplyr)
#' hbv_ru_1999$age <- trunc(hbv_ru_1999$age / 1) * 1
#' hbv_ru_1999$age[hbv_ru_1999$age > 40] <- trunc(
#'   hbv_ru_1999$age[hbv_ru_1999$age > 40] / 5
#' ) * 5
#' df <- hbv_ru_1999 %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(pos), tot = sum(tot))
#' plot(
#'   df$age, df$pos / df$tot,
#'   cex = 0.05 * df$tot, pch = 16, xlab = "age",
#'   ylab = "seroprevalence", xlim = c(0, 72)
#' )
#'
#' @source Mukomolov, S., L. Shliakhtenko, I. Levakova, and E. Shargorodskaya.
#' Viral hepatitis in Russian federation. An analytical overview. Technical
#' Report 213 (3), 3rd edn. St Petersburg Pasteur Institute,
#' St Petersburg, 2000.
"hbv_ru_1999"

#' Hepatitis C serological data from Belgium in 2006 (line listing)
#'
#' A study of HCV infection among injecting drug users. All injecting drug users
#' were interviewed by means of a standardized face-to-face interview and
#' information on their socio-demographic status, drug use history, drug use,
#' and related risk behavior was recorded
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{dur}{Duration of injection/Exposure time (years)}
#'  \item{seropositive}{If the individual is seropositive or not}
#' }
#'
#' @examples
#' # Reproduce Fig 4.3, p. 66
#' library(dplyr)
#' # snapping age to aggregated age group
#' # (credit: https://stackoverflow.com/a/12861810)
#' groups <- c(0.5:24.5)
#' range <- 0.5
#' low <- findInterval(hcv_be_2006$dur, groups)
#' high <- low + 1
#' low_diff <- hcv_be_2006$dur - groups[ifelse(low == 0, NA, low)]
#' high_diff <- groups[ifelse(high == 0, NA, high)] - hcv_be_2006$dur
#' mins <- pmin(low_diff, high_diff, na.rm = TRUE)
#' pick <- ifelse(!is.na(low_diff) & mins == low_diff, low, high)
#' hcv_be_2006$dur <- ifelse(
#'   mins <= range + .Machine$double.eps, groups[pick], hcv_be_2006$dur
#' )
#' hcv_be_2006 <- hcv_be_2006 %>%
#'   group_by(dur) %>%
#'   summarise(tot = n(), pos = sum(seropositive))
#'
#' plot(
#'   hcv_be_2006$dur, hcv_be_2006$pos / hcv_be_2006$tot,
#'   cex = 0.1 * hcv_be_2006$tot, pch = 16,
#'   xlab = "duration of injection (years)",
#'   ylab = "seroprevalence", xlim = c(0, 25), ylim = c(0, 1)
#' )
#'
#' @source Mathei, C., Shkedy, Z., Denis, B., Kabali, C., Aerts,
#' M., Molenberghs, G., Van Damme, P. and Buntinx, F. (2006), Evidence for
#' a substantial role of sharing of injecting paraphernalia other than
#' syringes/needles to the spread of hepatitis C among injecting drug users.
#' Journal of Viral Hepatitis, 13: 560-570.
#' \doi{doi:10.1111/j.1365-2893.2006.00725.x}
"hcv_be_2006"

#' Mumps serological data from the UK in 1986 and 1987 (aggregated)
#'
#' a large survey of prevalence of antibodies to mumps and rubella viruses in
#' the UK. The survey, covering subjects from 1 to over 65 years of age,
#' provides information on the prevalence of antibody by age
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#' }
#'
#' @examples
#' # Reproduce Fig 4.4 (left panel), p. 67
#' age <- mumps_uk_1986_1987$age
#' pos <- mumps_uk_1986_1987$pos
#' tot <- mumps_uk_1986_1987$tot
#' plot(age, pos / tot,
#'   cex = 0.008 * tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 45), ylim = c(0, 1)
#' )
#'
#' @source  Morgan-Capner P, Wright J, Miller C L, Miller E. Surveillance of
#' antibody to measles, mumps, and rubella by age. British Medical Journal
#' 1988; 297 :770 \doi{doi:10.1136/bmj.297.6651.770}
"mumps_uk_1986_1987"

#' Parvo B19 serological data from Belgium from 2001-2003 (line listing)
#'
#' A seroprevalence survey testing for parvovirus B19 IgG antibody, performed on
#' large representative national serum banks in Belgium, England and Wales,
#' Finland, Italy, and Poland. The sera were collected between 1995 and 2004 and
#' were obtained from residual sera submitted for routine laboratory testing.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#'  \item{year}{Year surveyed}
#'  \item{gender}{Gender of individual}
#'  \item{parvouml}{Parvo B19 antibody units per ml}
#' }
#'
#' @examples
#' # Reproduce Fig 4.5 (left upper panel), p. 68
#' library(dplyr)
#' df <- parvob19_be_2001_2003 %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(df$age, df$pos / df$tot,
#'   cex = 0.02 * df$tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 82), ylim = c(0, 1)
#' )
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"parvob19_be_2001_2003"

#' Parvo B19 serological data from England and Wales in 1996 (line listing)
#'
#' A seroprevalence survey testing for parvovirus B19 IgG antibody, performed on
#' large representative national serum banks in Belgium, England and Wales,
#' Finland, Italy, and Poland. The sera were collected between 1995 and 2004 and
#' were obtained from residual sera submitted for routine laboratory testing.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#'  \item{year}{Year surveyed}
#'  \item{gender}{Gender of individual}
#'  \item{parvouml}{Parvo B19 antibody units per ml}
#' }
#'
#' @examples
#' # Reproduce Fig 4.5 (right upper panel), p. 68
#' # NB: This figure will look different to that of in the book, since we
#' # believe that the original authors has made some errors in specifying
#' # the sample size of the dots.
#' library(dplyr)
#' df <- parvob19_ew_1996 %>%
#'   mutate(age = round(age)) %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(df$age, df$pos / df$tot,
#'   cex = 0.02 * df$tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 82), ylim = c(0, 1)
#' )
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"parvob19_ew_1996"

#' Parvo B19 serological data from Finland from 1997-1998 (line listing)
#'
#' A seroprevalence survey testing for parvovirus B19 IgG antibody, performed on
#' large representative national serum banks in Belgium, England and Wales,
#' Finland, Italy, and Poland. The sera were collected between 1995 and 2004 and
#' were obtained from residual sera submitted for routine laboratory testing.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#'  \item{year}{Year surveyed}
#'  \item{gender}{Gender of individual}
#'  \item{parvouml}{Parvo B19 antibody units per ml}
#' }
#'
#' @examples
#' # Reproduce Fig 4.5 (left bottom panel), p. 68
#' # NB: This figure will look different to that of in the book, since we
#' # believe that the original authors has made some errors in specifying
#' # the sample size of the dots.
#' library(dplyr)
#' df <- parvob19_fi_1997_1998 %>%
#'   mutate(age = round(age)) %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(df$age, df$pos / df$tot,
#'   cex = 0.07 * df$tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 82), ylim = c(0, 1)
#' )
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"parvob19_fi_1997_1998"

#' Parvo B19 serological data from Italy from 2003-2004 (line listing)
#'
#' A seroprevalence survey testing for parvovirus B19 IgG antibody, performed on
#' large representative national serum banks in Belgium, England and Wales,
#' Finland, Italy, and Poland. The sera were collected between 1995 and 2004 and
#' were obtained from residual sera submitted for routine laboratory testing.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#'  \item{year}{Year surveyed}
#'  \item{gender}{Gender of individual}
#'  \item{parvouml}{Parvo B19 antibody units per ml}
#' }
#'
#' @examples
#' # Reproduce Fig 4.5 (middle bottom panel), p. 68
#' # NB: This figure will look different to that of in the book, since we
#' # believe that the original authors has made some errors in specifying
#' # the sample size of the dots.
#' library(dplyr)
#' df <- parvob19_it_2003_2004 %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(df$age, df$pos / df$tot,
#'   cex = 0.07 * df$tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 82), ylim = c(0, 1)
#' )
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"parvob19_it_2003_2004"

#' Parvo B19 serological data from Poland from 1995-2004 (line listing)
#'
#' A seroprevalence survey testing for parvovirus B19 IgG antibody, performed on
#' large representative national serum banks in Belgium, England and Wales,
#' Finland, Italy, and Poland. The sera were collected between 1995 and 2004 and
#' were obtained from residual sera submitted for routine laboratory testing.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#'  \item{year}{Year surveyed}
#'  \item{gender}{Gender of individual}
#'  \item{parvouml}{Parvo B19 antibody units per ml}
#' }
#'
#' @examples
#' # Reproduce Fig 4.5 (right bottom panel), p. 68
#' # NB: This figure will look different to that of in the book, since we
#' # believe that the original authors has made some errors in specifying
#' # the sample size of the dots.
#' library(dplyr)
#' df <- parvob19_pl_1995_2004 %>%
#'   mutate(age = round(age)) %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(df$age, df$pos / df$tot,
#'   cex = 0.07 * df$tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 82), ylim = c(0, 1)
#' )
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"parvob19_pl_1995_2004"

#' Rubella serological data from the UK in 1986 and 1987 (aggregated)
#'
#' Prevalence of rubella in the UK, obtained from a large survey of prevalence
#' of antibodies to both mumps and rubella viruses.
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#' }
#'
#' @examples
#' # Reproduce Fig 4.4 (middle panel), p. 67
#' age <- rubella_uk_1986_1987$age
#' pos <- rubella_uk_1986_1987$pos
#' tot <- rubella_uk_1986_1987$tot
#' plot(age, pos / tot,
#'   cex = 0.008 * tot, pch = 16, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 45), ylim = c(0, 1)
#' )
#'
#' @source  Morgan-Capner P, Wright J, Miller C L, Miller E. Surveillance of
#' antibody to measles, mumps, and rubella by age. British Medical Journal
#' 1988; 297 :770 \doi{doi:10.1136/bmj.297.6651.770}
"rubella_uk_1986_1987"

#' Tuberculosis serological data from the Netherlands 1966-1973 (aggregated)
#'
#' A study of tuberculosis conducted in the Netherlands. Schoolchildren, aged
#' between 6 and 18 years, were tested using the tuberculin skin test.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#'  \item{gender}{Gender of cohort (unsure what 0 and 1 means)}
#'  \item{birthyr}{Birth year of cohort}
#' }
#'
#' @examples
#' # Reproduce Fig 4.6, p.70
#' age <- tb_nl_1966_1973$age
#' birthyr <- tb_nl_1966_1973$birthyr
#' pos <- tb_nl_1966_1973$pos
#' tot <- tb_nl_1966_1973$tot
#' # left panel
#' plot(age, pos / tot,
#'   pch = 16, cex = 0.00005 * tot, xlab = "age",
#'   ylab = "prevalence", xlim = c(6, 18)
#' )
#' # right panel
#' plot(birthyr, pos / tot,
#'   pch = 16, cex = 0.00005 * tot, xlab = "year", ylab = "prevalence"
#' )
#'
#' @source Nagelkerke, N., Heisterkamp, S., Borgdorff, M., Broekmans, J. and
#' Van Houwelingen, H. (1999), Semi-parametric estimation of age-time specific
#' infection incidence from serial prevalence data. Statist. Med., 18: 307-320.
#' \doi{doi:10.1002/(SICI)1097-0258(19990215)18:3<307::AID-SIM15>3.0.CO;2-Z}
"tb_nl_1966_1973"

#' VZV serological data from Belgium (Flanders) from 1999-2000 (aggregated)
#'
#' Age-specific seroprevalence of VZV antibodies, assessed in Flanders (Belgium)
#' between October 1999 and April 2000. This population was stratified by age
#' in order to obtain approximately 100 observations per age group.
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#' }
#'
#' @examples
#' # Reproduce Fig 4.7 (left panel), p.71
#' age <- vzv_be_1999_2000$age
#' pos <- vzv_be_1999_2000$pos
#' tot <- vzv_be_1999_2000$tot
#' plot(age, pos / tot,
#'   cex = 0.036 * tot, pch = 19, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 45), ylim = c(0, 1)
#' )
#'
#' @source Thiry, N., Beutels, P., Shkedy, Z. et al. The seroepidemiology
#' of primary varicella-zoster virus infection in Flanders (Belgium).
#' Eur J Pediatr 161, 588-593 (2002).
#' \doi{doi:10.1007/s00431-002-1053-2}
"vzv_be_1999_2000"

#' VZV serological data from Belgium from 2001-2003 (line listing)
#'
#' The survey is the same as the one used to study the seroprevalence of
#' parvovirus B19 in Belgium, as described above.
#'
#' @format A data frame with 4 variables:
#' \describe{
#'  \item{age}{Age of individual}
#'  \item{seropositive}{If the individual is seropositive or not}
#'  \item{gender}{Gender of individual}
#'  \item{VZVmIUml}{VZV milli international units per ml}
#' }
#'
#' @examples
#' # Reproduce Fig 4.7 (right panel), p.71
#' library(dplyr)
#' df <- vzv_be_2001_2003 %>%
#'   mutate(age = round(age)) %>%
#'   group_by(age) %>%
#'   summarise(pos = sum(seropositive), tot = n())
#' plot(df$age, df$pos / df$tot,
#'   cex = 0.036 * df$tot, pch = 19, xlab = "age", ylab = "seroprevalence",
#'   xlim = c(0, 45), ylim = c(0, 1)
#' )
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"vzv_be_2001_2003"

#' Rubella - Mumps data from the UK (aggregated)
#'
#'
#' @format A data frame with 5 variables:
#' \describe{
#'  \item{age}{Age group}
#'  \item{NN}{Number of individuals negative to rubella and mumps}
#'  \item{NP}{Number of individuals negative to rubella and positive to mumps}
#'  \item{PN}{Number of individuals positive to rubella and negative to mumps}
#'  \item{PP}{Number of individuals positive to rubella and mumps}
#' }
#'
#'
#' @source  Morgan-Capner P, Wright J, Miller C L, Miller E. Surveillance of
#' antibody to measles, mumps, and rubella by age. British Medical Journal
#' 1988; 297 :770 \doi{doi:10.1136/bmj.297.6651.770}
"rubella_mumps_uk"


#' VZV and Parvovirus B19 serological data in Belgium (line listing)
#' @format A data frame with 7 variables:
#' \describe{
#'  \item{id}{ID of individual}
#'  \item{age}{Age of individual}
#'  \item{gender}{Gender of individual}
#'  \item{parvouml}{Parvo B19 antibody units per ml}
#'  \item{parvo_res}{If an individual is positive for parvovirus B19}
#'  \item{VZVmUIml}{VZV milli international units per ml}
#'  \item{vzv_res}{If an individual is positive for VZV}
#' }
#'
#' @source MOSSONG, J., N. HENS, V. FRIEDERICHS, I. DAVIDKIN, M. BROMAN,
#' B. LITWINSKA, J. SIENNICKA, et al. "Parvovirus B19 Infection in Five
#' European Countries: Seroepidemiology, Force of Infection and Maternal Risk
#' of Infection." Epidemiology and Infection 136, no. 8 (2008): 1059-68.
#' \doi{doi:10.1017/S0950268807009661}
"vzv_parvo_be"
