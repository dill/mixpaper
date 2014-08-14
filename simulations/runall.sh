#!/bin/bash
export R_LIBS="~/R/library:${R_LIBS}"
export R_LIBS_USER="~/R/library:${R_LIBS}"


cd nocov
R -f runsim.R
cd ..

cd pt
R -f runsim.R
cd ..

cd 3point
R -f runsim-MS.R
cd ..

cd hazard
R -f runsim.R
cd ..

cd eps
R -f runsim.R
cd ..

cd covar
R -f covsim1-MS.R
R -f covsim1-MS-noadj.R
R -f covsim2-MS.R
R -f covsim2-MS-noadj.R
cd ..

cd hazard
R -f runsim.R
cd ..

