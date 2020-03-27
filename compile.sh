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
g++ -Wall -std=c++17 -O3  c++/component.hpp c++/diffusion.hpp c++/grid.hpp c++/rdc.hpp c++/reaction.hpp c++/nlsolver.hpp c++/main.cpp
echo "compilatioin succeeded"
cp a.out pear_diffusion
echo "executable named pear_diffusion in main folder"
rm a.out c++/component.hpp.gch c++/diffusion.hpp.gch c++/grid.hpp.gch c++/nlsolver.hpp.gch c++/rdc.hpp.gch c++/reaction.hpp.gch
echo "cleanup done"