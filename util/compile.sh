# !/bin/bash
# ./compile.sh

echo "Starting command line compilation with flags -Wall and -O3"
g++ -Wall -std=c++17 -O3    software/component.hpp \
                            software/diffusion.hpp \
                            software/grid.hpp \
                            software/rdc.hpp \
                            software/respiration.hpp \
                            software/reaction.hpp \
                            software/nlsolver.hpp \
                            software/main.cpp
echo "Compilation succeeded"
cp a.out pear_diffusion
echo "Executable named pear_diffusion in main folder"
rm a.out    software/component.hpp.gch \
            software/diffusion.hpp.gch \
            software/grid.hpp.gch \
            software/nlsolver.hpp.gch \
            software/rdc.hpp.gch \
            software/respiration.hpp.gch \
            software/reaction.hpp.gch
echo "Cleanup done"

