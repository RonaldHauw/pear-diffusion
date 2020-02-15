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

int main(int argc, char* argv[]){

    typedef double d_type; // change data type here

    std::cout<<"Ello fellas, com estas?"<<std::endl;
    const char* grid_name = "../conference_pear.txt";
    pear::grid<d_type> grid(grid_name);
    pear::component co2("CO_2");
    pear::component o2("O_2");
    pear::diffusion<d_type> diff_co2(co2, grid, 0.4);
    pear::diffusion<d_type> diff_o2(o2, grid, 1.0);
    pear::respiration_o2<d_type> react_o2(o2, co2, grid, 0.0, 0.0, 0.0);
    pear::respiration_co2<d_type> react_co2(co2, o2, grid, 0.0, 0.0, 0.0);



    pear::rdc<d_type> equation(o2, co2, diff_o2, diff_co2, react_o2, react_co2);

    // equation.add_component(o2)
    // equation.add_component(co2)
    // equation.add_diffusion(co2, diff_co2)
    // equation.add_diffusion(o2, diff_o2)
    // equation.add_reaction(o2, co2, reaction_o2)
    // equation.add_reaction(co2, o2, reaction_co2)
    // equation.solve(method="BFGS")
    // export_solution(equation.sol)




}