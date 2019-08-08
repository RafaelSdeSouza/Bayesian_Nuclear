######################################################################
# Author: Christian Iliadis (07/13/2019)
######################################################################
#
# ANALYSIS OF 7Be(n,p)7Li DATA
#
# REACTION MODEL:
# - incoherent sum of single-level, two-channel R-matrix expressions
#   plus a constant background term
# - boundary condition parameter is not varied independently; we
#   chose E_lambda = E_r
# - we allow for a constant reduced cross section background term
# - each resonance is sampled by multiplying R-matrix expression by
#   an "indicator" variable that is sampled according a Bernoulli
#   density [resonance is either needed or not needed, with equal
#   prior probability]
#
# DATA TREATMENT:
# - we split each data sets into two parts: a "relative cross
#   section" part, and an "absolute cross section" part; the relative
#   part is modeled with a normalization factor [e.g., y.norm1] for
#   which the prior is sampled from a very broad density; the absolute
#   part consists of just one data point. It does have a statistical
#   uncertainty. It is modeled with a normalization factor [e.g.,
#   y.norm11] for which the prior is sampled from a highly informative
#   lognormal density, where the lognormal parameters mu and sigma
#   a chosen according to the systematic uncertainty reported for a
#   given experiment
#
# OUTPUT:
# - in the plot of the reduced cross section, the red data points show
#   the "absolute" cross section of each data set
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library("rjags")
## for block updating [we do not need to center predictor variables]
load.module("glm")

## load external R-matrix functions
load.module("nuclear")
# random number seed
# set.seed(3000)
require(MASS)
require(RcppGSL)

######################################################################
## FUNCTIONS
######################################################################
#

######################################################################
## DATA INPUT: ENERGY, sqrt(E) * SIGMA
## [statistical uncertainties only!]
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b;

### RELATIVE CROSS SECTIONS:

## DATA SET 1: damone 2018;
##             includes 7Li ground and first excited state
tab1     <- read.table("exp_damone_rel.dat", header=FALSE)
obsx1    <- c(tab1[,1])
obsy1    <- c(tab1[,2])
errobsy1 <- c(tab1[,3])

## DATA SET 2: gibbons & macklin 1959; data with Eplab>2371 keV excluded;
##             also excluded data below Encm=0.01 MeV; use relativistic
##             results;
##             includes only 7Li ground state
tab2     <- read.table("exp_gibbons_rel.dat", header=FALSE)
obsx2    <- c(tab2[,1])
obsy2    <- c(tab2[,2])
errobsy2 <- c(tab2[,3])

## DATA SET 3: martin-hernandez 2018; absolute pn cross directly obtained
##             from author; converted to np using relativistic kinematics;
##             data below Encm=0.002 MeV are excluded because uncertainty
##             in 7Be mass impacts results;
##             includes only 7Li ground state
tab3     <- read.table("exp_martin_rel.dat", header=FALSE)
obsx3    <- c(tab3[,1])
obsy3    <- c(tab3[,2])
errobsy3 <- c(tab3[,3])

## DATA SET 4: koehler 1988; for energies below 8e-5 MeV [constant reduced
##             cross section] data have been rebinned by adding 5 data points
##             and calculating weighted mean;
##             includes 7Li ground and first excited state
tab4     <- read.table("exp_koehler_rel.dat", header=FALSE)
obsx4    <- c(tab4[,1])
obsy4    <- c(tab4[,2])
errobsy4 <- c(tab4[,3])

### ABSOLUTE NORMALIZATIONS:

## DATA SET 10: koehler 1988; thermal cross section
tab10     <- read.table("exp_koehler_abs.dat", header=FALSE)
obsx10    <- c(tab10[,1])
obsy10    <- c(tab10[,2])
errobsy10 <- c(tab10[,3])

## DATA SET 11: damone 2018; thermal cross section
tab11     <- read.table("exp_damone_abs.dat", header=FALSE)
obsx11    <- c(tab11[,1])
obsy11    <- c(tab11[,2])
errobsy11 <- c(tab11[,3])

## DATA SET 12: gibbons 1959; lowest-energy data point
tab12     <- read.table("exp_gibbons_abs.dat", header=FALSE)
obsx12    <- c(tab12[,1])
obsy12    <- c(tab12[,2])
errobsy12 <- c(tab12[,3])

## DATA SET 13: hernandez 2019; lowest-energy data point
tab13     <- read.table("exp_martin_abs.dat", header=FALSE)
obsx13    <- c(tab13[,1])
obsy13    <- c(tab13[,2])
errobsy13 <- c(tab13[,3])

## DATA SET 14: cervena 1989; thermal cross section
tab14     <- read.table("exp_cervena_abs.dat", header=FALSE)
obsx14    <- c(tab14[,1])
obsy14    <- c(tab14[,2])
errobsy14 <- c(tab14[,3])

## DATA SET 15: tomandl 2019; thermal cross section
tab15     <- read.table("exp_tomandl_abs.dat", header=FALSE)
obsx15    <- c(tab15[,1])
obsy15    <- c(tab15[,2])
errobsy15 <- c(tab15[,3])

######################################################################
# JAGS MODEL
######################################################################
# rjags ----->
cat('model {

######################################
#
# LIKELIHOODS
#
######################################
# - careful: dnorm is differently defined in R and JAGS!
# - precision=sigma^(-2)
# - in a for loop, make sure **all** variables on the LEFT of an
#   expression has the index [i]
# - systematic error as normalization factor y.norm...

## RELATIVE DATA

for (i in 1:length(obsx1)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy1[i] ~ dnorm(ya1[i], pow(errobsy1[i], -2))
  # ...subject to extrinsic scatter, if any:
  ya1[i] ~ dnorm(ym1[i], pow(yscat1, -2))
  # ...subject to syst uncertainties:
  ym1[i] <- yt1[i] * 10^y.norm1
  # true sigma [calculated from theory]:
  yt1[i] <- hbg + sqrt(obsx1[i]) *
            (
            sigma7Benpx(obsx1[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx1[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx1[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx1[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx1[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx1[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx1[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx2)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy2[i] ~ dnorm(ya2[i], pow(errobsy2[i], -2))
  # ...subject to extrinsic scatter, if any:
  ya2[i] ~ dnorm(ym2[i], pow(yscat2, -2))
  # ...subject to syst uncertainties:
  ym2[i] <- yt2[i] * 10^y.norm2
  # true sigma [calculated from theory]:
  yt2[i] <- hbg + sqrt(obsx2[i]) *
            (
            sigma7Benpx(obsx2[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx2[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx2[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx2[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx2[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx2[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx2[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx3)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy3[i] ~ dnorm(ya3[i], pow(errobsy3[i], -2))
  # ...subject to extrinsic scatter, if any:
  ya3[i] ~ dnorm(ym3[i], pow(yscat3, -2))
  # ...subject to syst uncertainties:
  ym3[i] <- yt3[i] * 10^y.norm3
  # true sigma [calculated from theory]:
  yt3[i] <- hbg + sqrt(obsx3[i]) *
            (
            sigma7Benpx(obsx3[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx3[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx3[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx3[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx3[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx3[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx3[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx4)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy4[i] ~ dnorm(ya4[i], pow(errobsy4[i], -2))
  # ...subject to extrinsic scatter, if any:
  ya4[i] ~ dnorm(ym4[i], pow(yscat4, -2))
  # ...subject to syst uncertainties:
  ym4[i] <- yt4[i] * 10^y.norm4
  # true sigma [calculated from theory]:
  yt4[i] <- hbg + sqrt(obsx4[i]) *
            (
            sigma7Benpx(obsx4[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx4[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx4[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx4[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx4[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx4[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx4[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

## ABSOLUTE DATA [no "extrinsic scatter" needed since there is only one data point]

for (i in 1:length(obsx10)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy10[i] ~ dnorm(ym10[i], pow(errobsy10[i], -2))
  # ...subject to syst uncertainties:
  ym10[i] <- yt10[i] * y.norm10
  # true sigma [calculated from theory]:
  yt10[i] <- hbg + sqrt(obsx10[i]) *
            (
            sigma7Benpx(obsx10[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx10[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx10[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx10[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx10[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx10[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx10[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx11)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy11[i] ~ dnorm(ym11[i], pow(errobsy11[i], -2))
  # ...subject to syst uncertainties:
  ym11[i] <- yt11[i] * y.norm11
  # true sigma [calculated from theory]:
  yt11[i] <- hbg + sqrt(obsx11[i]) *
            (
            sigma7Benpx(obsx11[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx11[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx11[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx11[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx11[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx11[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx11[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx12)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy12[i] ~ dnorm(ym12[i], pow(errobsy12[i], -2))
  # ...subject to syst uncertainties:
  ym12[i] <- yt12[i] * y.norm12
  # true sigma [calculated from theory]:
  yt12[i] <- hbg + sqrt(obsx12[i]) *
            (
            sigma7Benpx(obsx12[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx12[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx12[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx12[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx12[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx12[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx12[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx13)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy13[i] ~ dnorm(ym13[i], pow(errobsy13[i], -2))
  # ...subject to syst uncertainties:
  ym13[i] <- yt13[i] * y.norm13
  # true sigma [calculated from theory]:
  yt13[i] <- hbg + sqrt(obsx13[i]) *
            (
            sigma7Benpx(obsx13[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx13[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx13[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx13[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx13[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx13[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx13[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx14)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy14[i] ~ dnorm(ym14[i], pow(errobsy14[i], -2))
  # ...subject to syst uncertainties:
  ym14[i] <- yt14[i] * y.norm14
  # true sigma [calculated from theory]:
  yt14[i] <- hbg + sqrt(obsx14[i]) *
            (
            sigma7Benpx(obsx14[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx14[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx14[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx14[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx14[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx14[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx14[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx15)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy15[i] ~ dnorm(ym15[i], pow(errobsy15[i], -2))
  # ...subject to syst uncertainties:
  ym15[i] <- yt15[i] * y.norm15
  # true sigma [calculated from theory]:
  yt15[i] <- hbg + sqrt(obsx15[i]) *
            (
            sigma7Benpx(obsx15[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx15[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx15[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx15[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx15[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx15[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx15[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

######################################
#
# PRIORS
#
######################################
#
# NUCLEAR PRIORS
# parameters: ecm, e0, ga, gb, ra, rb, xj, xla, xlb

# channel radii
  ra ~ dnorm(4.0, pow(0.5, -2))T(0,)
  rb ~ dnorm(4.0, pow(0.5, -2))T(0,)
##  ra ~ dunif(3, 5)
##  rb ~ dunif(3, 5)
#  ra <- 5
#  rb <- 5

# Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
#
# for the reduced width priors, we will assume sigma = 0.5 * wl, so
# that Wigner limit corresponds to 2 * sigma
#
# neutron channel: wl_n = 41.80159/(0.88186 a_c^2) = 47.40160/a_c^2
# proton channel:  wl_p = 41.80159/(0.88079 a_c^2) = 47.45920/a_c^2
#
  wl_n <- 47.40160*pow(ra, -2)
  wl_p <- 47.45920*pow(rb, -2)

###########
# RESONANCE 1: 18.910 MeV, Jp=2-
# resonance spin, orbital angular momenta
  xj_1  <- 2
  xla_1 <- 0
  xlb_1 <- 0
# energy eigenvalue
  e0_1 ~ dnorm(0.0, pow(0.1, -2))T(0,)
# reduced widths
  ga_1 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_1 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# RESONANCE 2: 19.069 MeV, Jp=3+
# resonance spin, orbital angular momenta
  xj_2  <- 3
  xla_2 <- 1
  xlb_2 <- 1
# energy eigenvalue
  e0_2 ~ dnorm(0.150, pow(0.025, -2))T(0,)
# reduced widths
  ga_2 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_2 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# RESONANCE 3: 19.235 MeV, Jp=3+
# resonance spin, orbital angular momenta
  xj_3  <- 3
  xla_3 <- 1
  xlb_3 <- 1
# energy eigenvalue
  e0_3 ~ dnorm(0.336, pow(0.01, -2))T(0,)
# reduced widths
  ga_3 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_3 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# RESONANCE 4: 19.400 MeV, Jp=1-
# resonance spin, orbital angular momenta
  xj_4  <- 1
  xla_4 <- 0
  xlb_4 <- 0
# energy eigenvalue
  e0_4 ~ dnorm(0.51, pow(0.1, -2))T(0,)
# reduced widths
  ga_4 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_4 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# RESONANCE 5: 19.850 MeV, Jp=4+
# resonance spin, orbital angular momenta
  xj_5  <- 4
  xla_5 <- 3
  xlb_5 <- 3
# energy eigenvalue
  e0_5 ~ dnorm(0.96, pow(0.1, -2))T(0,)
# reduced widths
  ga_5 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_5 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# RESONANCE 6: 20.100 MeV, Jp=2+
# resonance spin, orbital angular momenta
  xj_6  <- 2
  xla_6 <- 1
  xlb_6 <- 1
# energy eigenvalue
  e0_6 ~ dnorm(1.23, pow(0.1, -2))T(0,)
# reduced widths
  ga_6 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_6 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# RESONANCE 7: 20.200 MeV, Jp=0+
# resonance spin, orbital angular momenta
  xj_7  <- 0
  xla_7 <- 1
  xlb_7 <- 1
# energy eigenvalue
  e0_7 ~ dnorm(1.32, pow(0.1, -2))T(0,)
# reduced widths
  ga_7 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)
  gb_7 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)

###########
# Constant background
##  hbg ~ dnorm(0.0, pow(0.2, -2))T(0,)
  hbg <- 0

######################################
#
# DATA SET PRIORS

###########
# RELATIVE DATA SETS
# extrinsic scatter
  yscat1 ~ dnorm(0.0, pow(2, -2))T(0,)     # 2 [sqrt(MeV)]b
  yscat2 ~ dnorm(0.0, pow(2, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(2, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(2, -2))T(0,)

# systematic normalization factor for S-factor:
# log(): natural logarithm

# chose very broad [weakly informative] for the data sets...
# sampled values are logs of multiplication factors
  y.norm1 ~ dunif(-1, 1)
  y.norm2 ~ dunif(-1, 1)
  y.norm3 ~ dunif(-1, 1)
  y.norm4 ~ dunif(-1, 1)

###########
## ABSOLUTE DATA SETS
# ...and informative for the absolute cross section points
  y.norm10 ~ dlnorm(logmu10, pow(logsigma10, -2))
  logmu10 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma10 <- log(1.020)  # factor uncertainty is 1.02, i.e. 2.0% for Koe88

  y.norm11 ~ dlnorm(logmu11, pow(logsigma11, -2))
  logmu11 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma11 <- log(1.10)   # factor uncertainty is 1.10, i.e. 10% for Dam18

  y.norm12 ~ dlnorm(logmu12, pow(logsigma12, -2))
  logmu12 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma12 <- log(1.050)  # factor uncertainty is 1.050, i.e. 5.0% for Gib59

  y.norm13 ~ dlnorm(logmu13, pow(logsigma13, -2))
  logmu13 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma13 <- log(1.051)  # factor uncertainty is 1.051, i.e. 5.1% for Mar19

  y.norm14 ~ dlnorm(logmu14, pow(logsigma14, -2))
  logmu14 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma14 <- log(1.085)  # factor uncertainty is 1.085, i.e. 8.5% for Cer89

  y.norm15 ~ dlnorm(logmu15, pow(logsigma15, -2))
  logmu15 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma15 <- log(1.032)  # factor uncertainty is 1.032, i.e. 3.2% for Tom19

}', file={f <- tempfile()})

# dunif(min=0, max=1): gives uniform density
# pow(a,b) = a^b
######################################################################
# n.adapt:  number of iterations in the chain for adaptation
#           [JAGS will use to choose the sampler and to assure optimum
#           mixing of the MCMC chain; will be discarded]
# update(): performs the burn-in on each chain by running the MCMC for
#           n.burn iterations without saving any of the posterior samples
# coda.samples(): runs each MCMC chain for the number of iterations
#           specified by n.iter, but it does not save every iteration;
#           instead, it saves only ever nth iteration, where n is given
#           by thin
# n.chains: number of mcmc chains

n.chains <- 3
n.adapt  <- 500
n.burn   <- 500
n.iter   <- 500
thin     <- 1

# "f": is the model specification from above;
# data = list(...): define all data elements that are referenced in the
# JAGS model;

ourmodel <- jags.model(f, data = list(
                 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1
                ,'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2
                ,'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3
                ,'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4
                ,'obsx10' = obsx10, 'obsy10' = obsy10, 'errobsy10' = errobsy10
                ,'obsx11' = obsx11, 'obsy11' = obsy11, 'errobsy11' = errobsy11
                ,'obsx12' = obsx12, 'obsy12' = obsy12, 'errobsy12' = errobsy12
                ,'obsx13' = obsx13, 'obsy13' = obsy13, 'errobsy13' = errobsy13
                ,'obsx14' = obsx14, 'obsy14' = obsy14, 'errobsy14' = errobsy14
                ,'obsx15' = obsx15, 'obsy15' = obsy15, 'errobsy15' = errobsy15
                                     ),
#               inits = list(e0 = 0.08, ga = 3.2, gb = 0.13, ra = 5.5, rb = 3.6),
                n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn)

# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel,
                 variable.names=c(
                  'e0_1', 'ga_1', 'gb_1'
                 ,'e0_2', 'ga_2', 'gb_2'
                 ,'e0_3', 'ga_3', 'gb_3'
                 ,'e0_4', 'ga_4', 'gb_4'
                 ,'e0_5', 'ga_5', 'gb_5'
                 ,'e0_6', 'ga_6', 'gb_6'
                 ,'e0_7', 'ga_7', 'gb_7'
                 ,'hbg', 'ra', 'rb'
                 ,'yscat1', 'y.norm1'
                 ,'yscat2', 'y.norm2'
                 ,'yscat3', 'y.norm3'
                 ,'yscat4', 'y.norm4'
                 ,'y.norm10', 'y.norm11', 'y.norm12'
                 ,'y.norm13', 'y.norm14', 'y.norm15'
#                 , 'z_1', 'z_2', 'z_3', 'z_4', 'z_5', 'z_6', 'z_7'
                                  ),
                 n.iter=n.iter, thin=thin)

# <---- rjags
######################################################################
######################################################################
# OUTPUT NUMERICAL RESULTS TO SCREEN
######################################################################
cat("", "\n")    # output empty line

# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain)
show(effectiveChainLength)

cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n")

######################################################################
# TRACES AND MARGINALIZED POSTERIORS
######################################################################
pdf("MCMC_7Benp_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# OUTPUT RESULTS TO FILES
######################################################################
# !!! make sure to check the order of the parameters in the MCMC output !!!

# matrix samplesmat contains the samples from all n.chains
# as.matrix() strips MCMC attributes from mcmc object
samplesmat = as.matrix(mcmcChain)
write.matrix(samplesmat,"7Benp_SAMP")

######################################################################



