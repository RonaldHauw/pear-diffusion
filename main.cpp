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

int main(){

    std::cout<<"Ello fellas, com estas?"<<std::endl;

    // grid = fem.grid(conference_pear.txt)
    // co2 = fem.component("CO_2")
    // o2 = fem.component("O_2")
    // diff_co2 = fem.diffusion(co2, grid, sigma = 0.5)
    // diff_o2 = fem.diffusion(o2, grid, sigma = 1.0)
    // mm_o2 = fem.michealismenten(alpha = 1.0, beta = 0.5, gamma = 2.5)
    // resp_co2 = fem.respiration(mm_o2, alpha = 1.0, beta = 0.5)
    // reaction_o2 = fem.reaction(o2, co2, grid, mm_o2)
    // reaction_co2 = fem.reaction(co2, o2, grid, resp_co2)
    // equation = fem.rdc()
    // equation.add_component(o2)
    // equation.add_component(co2)
    // equation.add_diffusion(co2, diff_co2)
    // equation.add_diffusion(o2, diff_o2)
    // equation.add_reaction(o2, co2, reaction_o2)
    // equation.add_reaction(co2, o2, reaction_co2)
    // equation.solve(method="BFGS")
    // export_solution(equation.sol)

    // illustration of some methods

    // grid.coord(number = 1)
    // >> tws::vector, length = 2: [0.0, 1.0]
    // diff_co2.evaluate(concentration = xco2)
    // >> tws::vector, dim = 1, length = N
    // diff_co2.jacobian(concentration = xco2)
    // >> tws::matrix, dim = 2, size = (N, N)




}