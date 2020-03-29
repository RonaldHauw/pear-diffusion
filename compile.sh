# !/bin/bash
# ./compile.sh

echo "starting command line compilation with flags -Wall and -O3"
#echo " - component.hpp"
#g++ -c ../component.hpp
#echo " - grid.hpp"
#g++ -c ../grid.hpp
#echo " - diffusion.hpp"
#g++ -c ../diffusion.hpp
#echo " - reaction.hpp"
#g++ -c ../reaction.hpp
#echo " - rdc.hpp"
#g++ -c ../rdc.hpp
#echo " - nlsolver.hpp"
#g++ -c ../nlsolver.hpp
echo " --> making executable"
g++ -Wall -std=c++17 -O3    software/component.hpp \
                            software/diffusion.hpp \
                            software/grid.hpp \
                            software/rdc.hpp \
                            software/reaction.hpp \
                            software/nlsolver.hpp \
                            software/main.cpp
echo "compilatioin succeeded"
cp a.out pear_diffusion
echo "executable named pear_diffusion in main folder"
rm a.out    software/component.hpp.gch \
            software/diffusion.hpp.gch \
            software/grid.hpp.gch \
            software/nlsolver.hpp.gch \
            software/rdc.hpp.gch \
            software/reaction.hpp.gch
echo "cleanup done"