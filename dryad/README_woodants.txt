Wood ant line transect data
===========================

Data used in paper "Mixture models for distance sampling detection functions" by David L. Miller and Len Thomas. Taken from Borkin et al (2012) and provided by Kerry Borkin.

Data are line transect observations of 144 wood ant (Formica lugubris and Formica aquilonia) nests in Abernethy Forest, Scotland. Both distance to ant nests and other covariates were collected.


Survey details
==============

Adapted from Borkin et al (2012):

Surveys took place on 13 days between between 22 August and 9 September 2003 using a line transect protocol on 17 transects in the Abernethy Forest, Scotland (57°15 ́N, 3°40 ́W).


Data format
===========

Comma separated value file with one row for each of the 144 observations, each with the following 7 columns:

habitat   : 4 level stand type: "stand initiation", "stem exclusion", "understorey re-initiation" and "old-growth".
transect  : ID of the transect where the observation was made.
distance  : Exact perpendicular distances to observed ant nests in metres.
nest.size : Size of nest, calculated as half-width multiplied by height.
species   : Factor with 2 levels ("aqu" and "lug") giving the ant species.
object    : Unique observation identifier.
detected  : Column of 1s indicated that the individual was observed (needed for analysis software).


References
==========

Borkin, KM, RW Summers, and L Thomas. Surveying Abundance and Stand Type Associations of Formica Aquilonia And F. Lugubris (Hymenoptera: Formicidae) Nest Mounds Over an Extensive Area: Trialing a Novel Method. European Journal of Entomology 109, no. 1 (2012): 47–53.
