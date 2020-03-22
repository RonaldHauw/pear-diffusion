#include <iostream>
#include "component.hpp"
#include "diffusion.hpp"
#include "grid.hpp"
#include "reaction.hpp"
#include "rdc.hpp"
#include <random>
#include "nlsolver.hpp"
#include "eigen/Eigen/Dense"

template<typename d_type, typename vec_type>
int export_solution(std::string const file_name, pear::grid<d_type> grid, std::vector<pear::component<d_type, vec_type>> components){
    for(int i = 0; i < components.size(); i++){
        pear::component<d_type, vec_type> comp = components[i];
        std::ofstream file;
        file.open(file_name+"_"+comp.name()+".txt");
        for (int n = 0; n < grid.nb_nodes(); n++){
            file
            <<n
            <<" "<<grid.node(n)[0]
            <<" "<<grid.node(n)[1]
            <<" "<<comp.cons()(n)
            <<std::endl;
        }
        file.close();
    }
    return 1;
}


int main(int argc, char* argv[]) {

    typedef double d_type; // change data type here
    typedef Eigen::Matrix<d_type, Eigen::Dynamic, 1> vec_type; // define vec_type here
    typedef Eigen::Matrix<d_type, Eigen::Dynamic, Eigen::Dynamic> mat_type; // define mat_type here

    std::cout << "Ello ello ello, com estas?" << std::endl;

    std::string grid_name = "matlab-prototype/mesh/HCTmesh3";
    pear::grid<d_type> grid(grid_name);

    // Diffusion parameters
    d_type sigma_u_r = 2.8e-10;
    d_type sigma_u_z = 1.10e-9;

    // to remove!
    //sigma_u_r = sigma_u_r*1e10;
    //sigma_u_z = sigma_u_z*1e10;

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

    // to remove!
    //r_u = r_u*1e7;

    d_type p_atm = 101300;
    d_type T_ref = 293.15;
    d_type R_g = 8.314;

    // Ambient conditions: Default is Orchard
    d_type T = 298.15; // 25 degrees Celsius
    d_type eta_u = 20.8e-2;
    d_type eta_v = 0.04e-2;


    int nl_maxit = 10;
    int set_environment = 0;
    d_type steplength = 1.;
    d_type alpha = 1.;
    for (int i = 0; i <argc; i++) {
        if (std::string(argv[i]) == "-ShelfLife" && set_environment == 0) {
            T = 293.15; // 20 degrees Celsius
            eta_u = 20.8e-2;
            eta_v = 0;
            std::cout<<"Shelflife environment chosen"<<std::endl;
            set_environment = 1;
        };

        if (std::string(argv[i]) == "-Refrigerator" && set_environment == 0) {
            T = 280.15; // 7 degrees Celsius
            eta_u = 20.8e-2;
            eta_v = 0;
            std::cout<<"Refrigerator environment chosen"<<std::endl;
            set_environment = 1;
        };

        if (std::string(argv[i]) == "-Precooling" && set_environment == 0) {
            T = 272.15; // -1 degrees Celsius
            eta_u = 20.8e-2;
            eta_v = 0;
            std::cout<<"Precooling environment chosen"<<std::endl;
            set_environment = 1;
        };

        if (std::string(argv[i]) == "-DisorderInducing" && set_environment == 0) {
            T = 272.15; // -1 degrees Celsius
            eta_u = 2e-2;
            eta_v = 5e-2;
            std::cout<<"DisorderInducing environment chosen"<<std::endl;
            set_environment = 1;
        };

        if (std::string(argv[i]) == "-OptimalCA" && set_environment == 0) {
            T = 272.15; // -1 degrees Celsius
            eta_u = 2e-2;
            eta_v = 0.7e-2;
            std::cout<<"OptimalCA environment chosen"<<std::endl;
            set_environment = 1;
        };
        if (argc - i > 1) { // at least two arguments remain
            if (std::string(argv[i]) == "-maxit") {
                nl_maxit = std::stoi(argv[i + 1]);
                std::cout<<"Setting maximum nonlinear iterations to: "<<nl_maxit<<std::endl;
            }
            if (std::string(argv[i]) == "-vmuref") {
                v_mu_ref = std::stod(argv[i + 1]);
                std::cout<<"Settig v_mu_ref to: "<<v_mu_ref<<std::endl;
            }
            if (std::string(argv[i]) == "-vmfvref") {
                v_mfv_ref = std::stod(argv[i + 1]);
                std::cout<<"Settig v_mfv_ref to: "<<v_mfv_ref<<std::endl;
            }
            if (std::string(argv[i]) == "-steplength") {
                steplength = std::stod(argv[i + 1]);
                std::cout<<"Settig steplength to: "<<steplength<<std::endl;
            }
            if (std::string(argv[i]) == "-anl") {
                alpha = std::stod(argv[i + 1]);
                std::cout<<"Adaptive nonlinearity: "<<alpha<<std::endl;
            }
        }

    };

    if (set_environment == 0){std::cout<<"Orchard environment chosen"<<std::endl;};


    // Respiration parameters
    d_type v_mu = v_mu_ref * exp(e_a_vmu_ref / R_g * (1. / T_ref - 1. / T));
    d_type v_mfv = v_mfv_ref * exp(e_a_vmfv_ref / R_g * (1. / T_ref - 1. / T));
    std::vector<d_type> respiration_param = {v_mu, v_mfv, k_mu, k_mv, k_mfu, r_q};

    // Boundary parameters // remove last hard coded number changed to zero!
    d_type c_u_amb = p_atm * eta_u / (R_g * T);
    d_type c_v_amb = p_atm * eta_v / (R_g * T);
    std::vector<d_type> diffusion_o2_param = {sigma_u_r, sigma_u_z, r_u, c_u_amb};
    std::vector<d_type> diffusion_co2_param = {sigma_v_r, sigma_v_z, r_v, c_v_amb};

    // allocate memory for the solution
    std::cout<<"pear::main(): allocating memory to store the solution"<<std::endl;
    std::cout<<"       - vec_type of size "<<grid.nb_nodes()*2<<std::endl;
    vec_type conc;
    conc.resize(grid.nb_nodes()*2, 1);
    //conc.setZero();

    // linearisatie
    // sparseheid patroon
    // niet homogene oplossing
    // test juiste implementatie Jacobiaan met eindige differenties

    // automatically deallocated pointers

    pear::component<d_type, vec_type> o2("O_2",   grid, conc, 0,        grid.nb_nodes(),     1);
    pear::component<d_type, vec_type> co2("CO_2", grid, conc, grid.nb_nodes(), grid.nb_nodes() * 2, 1);
    o2.set_initial(c_u_amb);
    co2.set_initial(c_v_amb);


    pear::diffusion<d_type, vec_type, mat_type> diff_o2(o2, grid, diffusion_o2_param);
    pear::diffusion<d_type, vec_type, mat_type> diff_co2(co2, grid, diffusion_co2_param);


    pear::respiration<d_type, vec_type, mat_type> resp_co2_o2(co2, o2, grid, respiration_param);

    pear::rdc<d_type, vec_type, mat_type> equation(diff_o2, diff_co2, resp_co2_o2);

    pear::nlsolver<d_type, pear::rdc<d_type, vec_type, mat_type>, vec_type, mat_type> nlsolve(equation);

    nlsolve.solve(nl_maxit, steplength, alpha);


    int exp_flag = export_solution("prototype/mesh/solution", grid, std::vector<pear::component<d_type, vec_type>>{o2, co2});

    if (exp_flag == 1){std::cout<<"solutions exported"<<std::endl; }

};
