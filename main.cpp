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
#include <random>
#include "nlsolver.hpp"
#include "eigen-3.3.7/Eigen/Dense"


int main(int argc, char* argv[]){

    typedef double d_type; // change data type here
    typedef Eigen::Matrix<d_type, Eigen::Dynamic, 1> vec_type; // define vec_type here
    typedef Eigen::Matrix<d_type, Eigen::Dynamic, Eigen::Dynamic> mat_type; // define mat_type here

    std::cout<<"Ello ello ello, com estas?"<<std::endl;

    std::string grid_name = "../grid/HCTmesh";
    pear::grid<d_type> grid(grid_name);

    // allocate memory for the solution
    vec_type conc;
    conc.resize(grid.length()+grid.length(), 1);




    pear::component<d_type, vec_type> co2("CO_2", grid, conc, 0, grid.length(), 1);
    pear::component<d_type, vec_type> o2("O_2", grid, conc, grid.length(), grid.length()*2, 1);
    pear::diffusion<d_type, vec_type> diff_co2(co2, grid, 0.4);
    pear::diffusion<d_type, vec_type> diff_o2(o2, grid, 1.0);
    pear::respiration<d_type, vec_type> resp_co2_o2(co2, o2, grid, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    pear::rdc<d_type, vec_type, mat_type> equation(diff_o2, diff_co2, resp_co2_o2);


    //pear::nlsolver<d_type, pear::rdc<d_type, vec_type>> nlsolve(equation);
    //nlsolve.solve();


    // nlsolve.solve(lasolver="cg", linesearch="wolfe");

    // export_solution(co2.concentration, o2.concentration)




}