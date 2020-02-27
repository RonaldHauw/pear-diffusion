#include <iostream>
#include "component.hpp"
#include "diffusion.hpp"
#include "grid.hpp"
#include "reaction.hpp"
#include "rdc.hpp"
#include <random>
#include "nlsolver.hpp"
#include "eigen-3.3.7/Eigen/Dense"







int main(int argc, char* argv[]) {

    typedef double d_type; // change data type here
    typedef Eigen::Matrix<d_type, Eigen::Dynamic, 1> vec_type; // define vec_type here
    typedef Eigen::Matrix<d_type, Eigen::Dynamic, Eigen::Dynamic> mat_type; // define mat_type here

    std::cout << "Ello ello ello, com estas?" << std::endl;

    std::string grid_name = "grid/HCTmesh";
    pear::grid<d_type> grid(grid_name);

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
    d_type r_u = 7e-7;
    d_type r_v = 7.5e-7;

    d_type p_atm = 101300;
    d_type T_ref = 293.15;
    d_type R_g = 8.314;

    // Ambient conditions: Default is Orchard
    d_type T = 298.15; // 25 degrees Celsius
    d_type eta_u = 20.8e-2;
    d_type eta_v = 0.04e-2;


    if (argc>1) {
        if (std::string(argv[1]) == "-ShelfLife") {
            T = 293.15; // 20 degrees Celsius
            eta_u = 20.8e-2;
            eta_v = 0;
        };

        if (std::string(argv[1]) == "-Refrigerator") {
            T = 280.15; // 7 degrees Celsius
            eta_u = 20.8e-2;
            eta_v = 0;
        };

        if (std::string(argv[1]) == "-Precooling") {
            T = 272.15; // -1 degrees Celsius
            eta_u = 20.8e-2;
            eta_v = 0;
        };

        if (std::string(argv[1]) == "-DisorderInducing") {
            T = 272.15; // -1 degrees Celsius
            eta_u = 2e-2;
            eta_v = 5e-2;
        };

        if (std::string(argv[1]) == "-OptimalCA") {
            T = 272.15; // -1 degrees Celsius
            eta_u = 2e-2;
            eta_v = 0.7e-2;
        };
    }; // if argc>1


    // Respiration parameters
    d_type v_mu = v_mu_ref * exp(e_a_vmu_ref / R_g * (1 / T_ref - 1 / T));
    d_type v_mfv = v_mfv_ref * exp(e_a_vmfv_ref / R_g * (1 / T_ref - 1 / T));
    std::vector<d_type> respiration_param = {v_mu, v_mfv, k_mu, k_mv, k_mfu, r_q};

    // Boundary parameters
    d_type c_u_amb = p_atm * eta_u / (R_g * T);
    d_type c_v_amb = p_atm * eta_v / (R_g * T);
    std::vector<d_type> diffusion_o2_param = {sigma_u_r, sigma_u_z, r_u, c_u_amb};
    std::vector<d_type> diffusion_co2_param = {sigma_v_r, sigma_v_z, r_v, c_v_amb};

    // allocate memory for the solution
    vec_type conc;
    conc.resize(grid.nb_nodes(), 1);


    std::cout << "Ola fellas, com estas?" << std::endl;

    pear::component<d_type, vec_type> o2("CO_2", grid, conc, 0, grid.nb_nodes(), 1);
    pear::component<d_type, vec_type> co2("O_2", grid, conc, grid.nb_nodes(), grid.nb_nodes() * 2, 1);
    pear::diffusion<d_type, vec_type, mat_type> diff_o2(o2, grid, diffusion_o2_param);
    pear::diffusion<d_type, vec_type, mat_type> diff_co2(co2, grid, diffusion_co2_param);


    pear::respiration<d_type, vec_type, mat_type> resp_co2_o2(co2, o2, grid, respiration_param);

    pear::rdc<d_type, vec_type, mat_type> equation(diff_o2, diff_co2, resp_co2_o2);

    mat_type K;
    K.resize(grid.nb_nodes() + grid.nb_nodes(), grid.nb_nodes() + grid.nb_nodes());




    pear::nlsolver<d_type, pear::rdc<d_type, vec_type, mat_type>, vec_type, mat_type> nlsolve(equation);
    nlsolve.solve();

    //std::cout<<conc<<std::endl;


    // nlsolve.solve(lasolver="cg", linesearch="wolfe");

    // export_solution(co2.concentration, o2.concentration)


};
