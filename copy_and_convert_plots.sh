#!/bin/sh

# copy and convert pdf plots to eps

# simulation plots
cp simulations/pa-plot-mmds.pdf Figure3.pdf
cp simulations/pa-plot-cds.pdf Figure-S1.pdf
cp simulations/pa-plot-combined.pdf Figure-S2.pdf
pdftops -eps -pagecrop Figure-S1.pdf
pdftops -eps -pagecrop Figure-S2.pdf
pdftops -eps -pagecrop Figure3.pdf
mv Figure-S1.eps PLOS\ submission/
mv Figure-S2.eps PLOS\ submission/
mv Figure3.eps PLOS\ submission/

# copy pdf files
cp dsmixtures-appendix-s1.pdf PLOS\ submission/
cp dsmixtures-appendix-s2.pdf PLOS\ submission/
cp dsmixtures-text-s1.pdf PLOS\ submission/
cp dsmixtures-text-s2.pdf PLOS\ submission/
cp dsmixtures.pdf PLOS\ submission/
