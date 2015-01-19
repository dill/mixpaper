Amakihi point transect data
===========================

Data used in paper "Mixture models for distance sampling detection functions" by David L. Miller and Len Thomas. Taken from Marques et al (2007) and were provided by Steve Fancy.

These are distance sampling point transect observations of 1485 Hawaii Amakihi (Hemignathus virens). Data include distances and covariates (see below).

Survey details
==============


Adapted from the data description in Marques et al (2007):

Data were collected as part of a larger survey on the island of Hawaii (Fancy et al. 1997). Multispecies point transect surveys were performed at 7 survey periods between July 1992 and April 1995. There were 41 point count stations, although they were not all surveyed in some survey periods.


Data format
===========

Comma separated value file with one row for each of the 1485 observations, each with the following 7 columns:

survey   : Factor with 7 levels giving the date of the survey (e.g. "Apr 93").
object   : Unique identifier per observation.
distance : Exact distances to observed birds in metres.
obs      : Factor with 3 levels giving observer ("SGF","TJS","TKP")
mas      : Minutes after sunrise that the observation took place.
has      : Hours (integer) after sunrise that the observation took place.
observed : Column of 1s indicated that the individual was observed (needed for analysis software).



References
==========

Fancy, SG, TJ Snetsinger, and JD Jacoby. Translocation of the Palila, an endangered Hawaiian honeycreeper. Pacific Conservation Biology 3 (1997) 39-46.

Marques, TA, L Thomas, SG Fancy, and ST Buckland. Improving Estimates of Bird Density Using Multiple-Covariate Distance Sampling. The Auk 124, no. 4 (2007): 1229–43.
