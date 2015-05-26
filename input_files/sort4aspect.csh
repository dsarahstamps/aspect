#!/bin/csh

if ($#argv < 2) then
  echo ""
  echo "This script will resampled a file and format it for an aspect input file."
  echo ""
  echo "usage: [in_file] [out_file]"
  echo ""
  echo "[1] [in_file] = Input file assumes we start with lat, lon, value format."
  echo "[2] [out_file] = Output file. Include the file extension."
  echo ""
  exit
endif

# Assuming that we are starting with a file in lat, lon, value 
# grid to 3 arc minutes using surface (set .gmtdefaults to 6 decimal points)
surface $1 -G$2.grd -: -I3m -R25/35.95/-11.983300/6.98333
grd2xyz $2.grd > $2.tmp 

# grd2xyz will put in lon lat value format 
awk '{print $2,$1,$3}' $2.tmp | sort -k1,1nr -k2,2 > $2 
