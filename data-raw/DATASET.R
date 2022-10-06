library(dplyr)
library(magrittr)
library(usethis)


################ HEPATITIS A ################
## HAV in Belgium (Flanders) in 1993-1994
hav_be_1993_1994 <- read.table(
    "data-raw\\serobook_datasets\\HAV-BE.dat",
    header = TRUE
)
# Neg means seropositive! (possible mislabelling by original authors)
hav_be_1993_1994 <- hav_be_1993_1994 %>%
    mutate(age = Age, pos = Neg, tot = Tot, .keep = "none")
use_data(hav_be_1993_1994, overwrite = TRUE)

## HAV in Belgium (Flanders) in 2002
hav_be_2002 <- read.table(
    "data-raw\\serobook_datasets\\hepatitis1993-2002.dat",
    header = TRUE
)
hav_be_2002 <- hav_be_2002 %>%
    filter(birth.cohort + age == 2002) %>%
    mutate(age = age, seropositive = status, .keep = "none")
use_data(hav_be_2002, overwrite = TRUE)

## HAV in Bulgaria in 1964
hav_bg_1964 <- read.table(
    "data-raw\\serobook_datasets\\HAV-BUL.dat",
    header = TRUE
)
hav_bg_1964 <- hav_bg_1964 %>%
    mutate(age = Age, pos = Pos, tot = Tot, .keep = "none")
use_data(hav_bg_1964, overwrite = TRUE)


################ HEPATITIS B ################
## HBV in Russia (St. Petersburg) in 1999
hbv_ru_1999 <- read.table(
    "data-raw\\serobook_datasets\\seroprevalencedataHepB.txt",
    header = TRUE
)
hbv_ru_1999 <- hbv_ru_1999 %>%
    relocate(pos, .after = age)
use_data(hbv_ru_1999, overwrite = TRUE)


################ HEPATITIS C ################
## HCV in Belgium in 2006
hcv_be_2006 <- read.table(
    "data-raw\\serobook_datasets\\hcvdat.txt",
    header = FALSE
)
hcv_be_2006 <- hcv_be_2006 %>%
    mutate(dur = V1, pos = V2 * V4, tot = V2, .keep = "none")
use_data(hcv_be_2006, overwrite = TRUE)


################ MUMPS ################
## Mumps in UK in 1986-1987
mumps_uk_1986_1987 <- read.table(
    "data-raw\\serobook_datasets\\MUMPSUK.dat",
    header = TRUE
)
mumps_uk_1986_1987 <- mumps_uk_1986_1987 %>%
    mutate(age = age, pos = pos, tot = ntot, .keep = "none")
use_data(mumps_uk_1986_1987, overwrite = TRUE)


################ PARVOVIRUS B19 ################
parvob19_eu_1995_2004 <- read.table(
    "data-raw\\serobook_datasets\\B19-countries.dat",
    header = TRUE
)
parvob19_eu_1995_2004 <- parvob19_eu_1995_2004 %>%
    mutate(
        age = age,
        seropositive = parvores,
        year = stringr::str_split(year, "_", simplify = TRUE)[1, 1],
        gender = sex,
        parvouml = parvouml,
        country = country,
        .keep = "none"
    ) %>%
    relocate(seropositive, .after = age) %>%
    relocate(parvouml, .after = gender) %>%
    filter(!is.na(seropositive), !is.na(age))

## Parvo B19 in Belgium in 2001-2003
parvob19_be_2001_2003 <- parvob19_eu_1995_2004 %>%
    filter(country == "be") %>%
    mutate(country = NULL)
use_data(parvob19_be_2001_2003, overwrite = TRUE)

## Parvo B19 in England and Wales in 1996
parvob19_ew_1996 <- parvob19_eu_1995_2004 %>%
    filter(country == "ew") %>%
    mutate(country = NULL)
use_data(parvob19_ew_1996, overwrite = TRUE)

## Parvo B19 in Finland in 1997-1998
parvob19_fi_1997_1998 <- parvob19_eu_1995_2004 %>%
    filter(country == "fi") %>%
    mutate(country = NULL)
use_data(parvob19_fi_1997_1998, overwrite = TRUE)

## Parvo B19 in Italy in 2003-2004
parvob19_it_2003_2004 <- parvob19_eu_1995_2004 %>%
    filter(country == "it") %>%
    mutate(country = NULL)
use_data(parvob19_it_2003_2004, overwrite = TRUE)

## Parvo B19 in Poland in 1995-2004
parvob19_pl_1995_2004 <- parvob19_eu_1995_2004 %>%
    filter(country == "pl") %>%
    mutate(country = NULL)
use_data(parvob19_pl_1995_2004, overwrite = TRUE)


################ RUBELLA ################
## Rubella in the UK in 1986-1987
rubella_uk_1986_1987 <- read.table(
    "data-raw\\serobook_datasets\\Rubella-UK.dat",
    header = TRUE
)
rubella_uk_1986_1987 <- rubella_uk_1986_1987 %>%
    mutate(age = Age, pos = Pos, tot = Pos + Neg, .keep = "none")
use_data(rubella_uk_1986_1987, overwrite = TRUE)


################ TUBERCULOSIS ################
## Tuberculosis in the Netherlands in 1966-1973
tb_nl_1966_1973 <- read.table(
    "data-raw\\serobook_datasets\\tb.dat",
    header = TRUE
)
tb_nl_1966_1973 <- tb_nl_1966_1973 %>%
    mutate(
        age = AGE, pos = PPD, tot = N,
        gender = SEX, birthyr = 1900 + BRTHYR, .keep = "none"
    )
use_data(tb_nl_1966_1973, overwrite = TRUE)


################ VARICELLA ZOSTER VIRUS ################
## VZV in Belgium (Flanders) in 1999-2000
vzv_be_1999_2000 <- read.table(
    "data-raw\\serobook_datasets\\VZV-Flanders.dat",
    header = TRUE
)
vzv_be_1999_2000 <- vzv_be_1999_2000 %>%
    mutate(age = age, pos = pos, tot = ntot, .keep = "none")
use_data(vzv_be_1999_2000, overwrite = TRUE)

## VZV in Belgium in 2001-2003
vzv_be_2001_2003 <- read.table(
    "data-raw\\serobook_datasets\\VZV-B19-BE.dat",
    header = TRUE
)
vzv_be_2001_2003 <- vzv_be_2001_2003 %>%
    mutate(
        age = age, seropositive = VZVres,
        gender = sex, VZVmIUml = VZVmUIml, .keep = "none"
    ) %>%
    filter(!is.na(seropositive)) %>%
    relocate(gender, .after = seropositive)
use_data(vzv_be_2001_2003, overwrite = TRUE)
