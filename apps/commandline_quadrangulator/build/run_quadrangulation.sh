#!/bin/bash

# Directory containing the .obj files
INPUT_DIR="./test_data/with_features"

./quad_commandline ../test_data/no_features/aquadom.obj -nofeatures
./quad_commandline ../test_data/no_features/beetle.obj -nofeatures
./quad_commandline ../test_data/no_features/botanic.obj -nofeatures
./quad_commandline ../test_data/no_features/bunny.obj -nofeatures
./quad_commandline ../test_data/no_features/calla.obj -nofeatures
./quad_commandline ../test_data/no_features/chair.obj -nofeatures
./quad_commandline ../test_data/no_features/chinese_lion24K.obj -nofeatures
./quad_commandline ../test_data/no_features/dome.obj -nofeatures
./quad_commandline ../test_data/no_features/ellipse.obj -nofeatures
./quad_commandline ../test_data/no_features/ellipsetorus.obj -nofeatures
./quad_commandline ../test_data/no_features/fertility.obj -nofeatures
./quad_commandline ../test_data/no_features/funnel.obj -nofeatures
./quad_commandline ../test_data/no_features/indorelax.obj -nofeatures
./quad_commandline ../test_data/no_features/lilium_2.obj -nofeatures
./quad_commandline ../test_data/no_features/manta.obj -nofeatures
./quad_commandline ../test_data/no_features/six_s0.obj -nofeatures
./quad_commandline ../test_data/no_features/spring.obj -nofeatures
./quad_commandline ../test_data/no_features/trebol.obj -nofeatures
./quad_commandline ../test_data/no_features/twist.obj -nofeatures


./quad_commandline ../test_data/with_features/12.obj -angle 70
./quad_commandline ../test_data/with_features/Bracket.obj -angle 70
./quad_commandline ../test_data/with_features/anti_backlash_nut.obj
./quad_commandline ../test_data/with_features/bathtub.obj -angle 70
./quad_commandline ../test_data/with_features/blockLP.obj
./quad_commandline ../test_data/with_features/bottle.obj
./quad_commandline ../test_data/with_features/bracket_CAD.obj
./quad_commandline ../test_data/with_features/fandisk.obj
./quad_commandline ../test_data/with_features/fusee_Lp.obj -angle 40
./quad_commandline ../test_data/with_features/joint.obj
./quad_commandline ../test_data/with_features/mount.obj
./quad_commandline ../test_data/with_features/mouse.obj
./quad_commandline ../test_data/with_features/part20k.obj
./quad_commandline ../test_data/with_features/wineglass.obj

echo "✅ All files processed."
