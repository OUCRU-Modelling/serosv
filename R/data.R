#' Serological survey data of infectious diseases
#'
#' These samples are cross-sectional surveys providing information about
#' whether or not the individual has been infected before that time point.
#' The data is taken from "Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data" - N. Hens et al. 2013
#' data are normalized and cleaned for the purpose of this package


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
#' @source <https://doi.org/10.1023/A:1007393405966>
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
#' @source <https://doi.org/10.1007/s00431-002-1053-2>
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
#' @source <https://doi.org/10.2307/2983150>
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
#' @source <Mukomolov, S., et al. Viral hepatitis in Russian federation. An
#' analytical overview. Technical Report 213 (3), 3rd edn. St Petersburg Pasteur
#' Institute, St Petersburg, 2000.>
"hbv_ru_1999"

#' Hepatitis C serological data from Belgium in 2006 (aggregated)
#'
#' A study of HCV infection among injecting drug users. All injecting drug users
#' were interviewed by means of a standardized face-to-face interview and
#' information on their socio-demographic status, drug use history, drug use,
#' and related risk behavior was recorded
#'
#' @format A data frame with 3 variables:
#' \describe{
#'  \item{dur}{Duration of injection/Exposure time (years)}
#'  \item{pos}{Number of seropositive individuals}
#'  \item{tot}{Total number of individuals surveyed}
#' }
#' @source <https://doi.org/10.1111/j.1365-2893.2006.00725.x>
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
#' @source <https://doi.org/10.1136/bmj.297.6651.770>
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
#' @source <http://doi.org/10.1017/S0950268807009661>
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
#' @source <http://doi.org/10.1017/S0950268807009661>
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
#' @source <http://doi.org/10.1017/S0950268807009661>
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
#' @source <http://doi.org/10.1017/S0950268807009661>
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
#' @source <http://doi.org/10.1017/S0950268807009661>
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
#' @source <https://doi.org/10.1136/bmj.297.6651.770>
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
#' @source <https://doi.org/10.1002/(SICI)1097-0258(19990215)18:3%3C307::AID-SIM15%3E3.0.CO;2-Z>
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
#' @source <https://doi.org/10.1007/s00431-002-1053-2>
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
#' @source <http://doi.org/10.1017/S0950268807009661>
"vzv_be_2001_2003"
