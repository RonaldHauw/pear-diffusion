/* Main file for project mathematical engineering: reaction-diffusion model of a pear.
 *
 * Authors:
 *  - Octave Oliviers
 *  - Agathe van Lamsweerde
 *  - Ronald Hauwaerts
 *
 * How to compile:
 *  - Via the make file
 *      >> make pear
 *  - Via shell script
 *      >> ./compile.sh
 */

#include <iostream>
#include "component.hpp"
#include "diffusion.hpp"
#include "grid.hpp"
#include "reaction.hpp"
#include "rdc.hpp"
#include "nlsolver.hpp"

int main(int argc, char* argv[]){

    typedef double d_type; // change data type here

    std::cout<<"Ello fellas, com estas?"<<std::endl;
    const char* grid_name = "../conference_pear.txt";
    pear::grid<d_type> grid(grid_name);
    pear::component<d_type> co2("CO_2", grid);
    pear::component<d_type> o2("O_2", grid);
    pear::diffusion<d_type> diff_co2(co2, grid, 0.4);
    pear::diffusion<d_type> diff_o2(o2, grid, 1.0);
    pear::respiration_o2<d_type> react_o2(o2, co2, grid, 0.0, 0.0, 0.0);
    pear::respiration_co2<d_type> react_co2(co2, o2, grid, 0.0, 0.0, 0.0);



    pear::rdc<d_type> equation(o2, co2, diff_o2, diff_co2, react_o2, react_co2);

    pear::nlsolver<d_type, pear::rdc<d_type>> nlsolve(equation);


    // nlsolve.solve(lasolver="cg", linesearch="wolfe");

    // export_solution(co2.concentration, o2.concentration)




}