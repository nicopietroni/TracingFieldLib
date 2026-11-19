#!/bin/bash

# optional angle for the features runs: leave empty to use the default of
# quad_commandline, or set a value (e.g. ANGLE=30), or pass it as first
# argument: ./run_all_quadrangulation.sh 30
ANGLE="${1:-}"

#abc
./run_quadrangulation_organic_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/abc/symm/no_features/
./run_quadrangulation_features_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/abc/symm/features/ $ANGLE
#./run_quadrangulation_organic.sh  /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/abc/no_symm/no_features/
./run_quadrangulation_features.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/abc/no_symm/features/ $ANGLE

#moldex
./run_quadrangulation_features.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/Moldex/no_symm/ $ANGLE
./run_quadrangulation_features_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/Moldex/symm/ $ANGLE

#tingi
./run_quadrangulation_features_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/tingi/symm/features/ $ANGLE

#tractive
./run_quadrangulation_organic_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/Tractive/symm/

#ansys
./run_quadrangulation_features.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/Ansys/ $ANGLE

#300
./run_quadrangulation_organic_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/300/symm/no_features/
./run_quadrangulation_features_symm.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/300/symm/features/ $ANGLE
./run_quadrangulation_organic.sh  /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/300/no_symm/no_features/
./run_quadrangulation_features.sh /Users/134750/Desktop/00_GeoProjects/Quadrangulation_Test/300/no_symm/features/ $ANGLE


