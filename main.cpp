#include <iostream>
#include "component.hpp"
#include "diffusion.hpp"
#include "grid.hpp"
#include "reaction.hpp"
#include "rdc.hpp"
#include "nlsolver.hpp"

int main(int argc, char* argv[]){

    typedef double d_type; // change data type here
    typedef std::vector<d_type> v_type; // change data type here

    // Diffusion parameters
    d_type sigma_u_r = 2.8e-10;
    d_type sigma_u_z = 1.10e-9;

    d_type sigma_v_r = 2.32e-9;
    d_type sigma_v_z = 6.97e-9;

    // Respiration parameters
    d_type v_mu_ref = 2.39e-4;
    d_type e_a_vmu_ref = 80200;
    d_type v_mfv_ref = 1.61e-4;
    d_type e_a_vmfv_ref = 56700;

    d_type k_mu = 0.4103;
    d_type k_mv = 27.2438;
    d_type k_mfu = 0.1149;
    d_type r_q = 0.97;

    // Boundary parameters
    d_type varrho_u = 7e-7;
    d_type varrho_v = 7.5e-7;

    d_type p_atm = 101300;
    d_type T_ref = 293.15;
    d_type R_g = 8.314;

    // Ambient conditions: Default is Orchard
    d_type T = 298.15; // 25 degrees Celsius
    d_type eta_u = 20.8e-2;
    d_type eta_v = 0.04e-2;

    if (std::string(argv[1]) == "-ShelfLife") {
        // Ambient conditions
        d_type T = 293.15; // 25 degrees Celsius
        d_type eta_u = 20.8e-2;
        d_type eta_v = 0;
    };

    if (std::string(argv[1]) == "-Refrigerator") {
        // Ambient conditions
        d_type T = 280.15; // 25 degrees Celsius
        d_type eta_u = 20.8e-2;
        d_type eta_v = 0;
    };

    if (std::string(argv[1]) == "-Precooling") {
        // Ambient conditions
        d_type T = 272.15; // 25 degrees Celsius
        d_type eta_u = 20.8e-2;
        d_type eta_v = 0;
    };

    if (std::string(argv[1]) == "-DisorderInducing") {
        // Ambient conditions
        d_type T = 272.15; // 25 degrees Celsius
        d_type eta_u = 2e-2;
        d_type eta_v = 5e-2;
    };

    if (std::string(argv[1]) == "-OptimalCA") {
        // Ambient conditions
        d_type T = 272.15; // 25 degrees Celsius
        d_type eta_u = 2e-2;
        d_type eta_v = 0.7e-2;
    };

    for (int i = 2; i < argc; ++i) {
        if (i+1<argc){
            if (std::string(argv[i]) == "-sigma_u_r") {sigma_u_r = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-sigma_u_z") {sigma_u_z = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-sigma_v_r"){sigma_v_r = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-sigma_v_z") {sigma_v_z = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-v_mu_ref") {v_mu_ref = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-e_a_vmu_ref") {e_a_vmu_ref = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-v_mfv_ref") {v_mfv_ref = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-e_a_vmfv_ref") {e_a_vmfv_ref =  std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-k_mu") {k_mu =  std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-k_mv") {k_mv = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-k_mfu") {k_mfu = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-r_q") {r_q = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-varrho_u") {varrho_u = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-varrho_v") {varrho_v = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-p_atm") {p_atm = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-T") {T = 273.15 + std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-eta_u") {varrho_v = std::stod(std::string(argv[i+1]));i++;};
            if (std::string(argv[i]) == "-eta_v") {varrho_v = std::stod(std::string(argv[i+1]));i++;};
        } // end if
    } // end for

    // Respiration parameters
    d_type v_mu = v_mu_ref * exp(e_a_vmu_ref/R_g * (1/T_ref - 1/T) );
    d_type v_mfv = v_mfv_ref * exp(e_a_vmfv_ref/R_g * (1/T_ref - 1/T) );
    v_type respiration_param{v_mu, v_mfv, k_mu, k_mv, k_mfu, r_q};

    // Boundary parameters
    d_type c_u_amb = p_atm * eta_u / (R_g*T);
    d_type c_v_amb = p_atm * eta_v / (R_g*T);

    std::cout<<"Ola fellas, com estas?"<<std::endl;
    const char* grid_name = "../conference_pear.txt";
    pear::grid<d_type> grid(grid_name);
    pear::component<d_type> co2("CO_2", grid);
    pear::component<d_type> o2("O_2", grid);
    pear::diffusion<d_type> diff_co2(co2, grid, sigma_u_r, sigma_u_z);
    pear::diffusion<d_type> diff_o2(o2, grid, sigma_v_r, sigma_v_z);
    pear::respiration_o2<d_type> react_o2(o2, co2, grid, 0.0, 0.0, 0.0);
    pear::respiration_co2<d_type> react_co2(co2, o2, grid, 0.0, 0.0, 0.0);



    pear::rdc<d_type> equation(o2, co2, diff_o2, diff_co2, react_o2, react_co2);

    pear::nlsolver<d_type, pear::rdc<d_type>> nlsolve(equation);


    // nlsolve.solve(lasolver="cg", linesearch="wolfe");

    // export_solution(co2.concentration, o2.concentration)


int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
