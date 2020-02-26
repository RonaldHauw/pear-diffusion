#ifndef REACTION_HPP
#define REACTION_HPP

#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "grid.hpp"


namespace pear {


    template <typename d_type, typename vec_type, typename mat_type>
    class respiration{
    public:

        respiration(pear::component<d_type, vec_type> & co2,
                pear::component<d_type, vec_type> & o2,
                pear::grid<d_type> & grid,
                std::vector<d_type> & param)
                : o2_(o2)
                , co2_(co2)
                , grid_(grid)
                , v_mu_(param[0])
                , v_mfv_(param[1])
                , k_mu_(param[2])
                , k_mv_(param[3])
                , k_mfu_(param[4])
                , r_q_(param[5])
        {
            std::cout<<"Respiration between O2 and CO2"<<std::endl;
        }

        d_type R_u(d_type C_u, d_type C_v){
            return v_mu_*C_u/(k_mu_+C_u)/(1+C_v/k_mv_);
        };

        d_type R_v(d_type C_u, d_type C_v){
            return r_q_*R_u(C_u, C_v) + v_mfv_/(1+C_u/k_mfu_);
        };

        void f_o2(vec_type & H) {
            for (int t = 1; t < grid_.nb_elements(); t++) {
                std::vector<d_type> elem_nodes = grid_.element(t);

                d_type r1 = grid_.node(elem_nodes[0])[0];   d_type z1 = grid_.node(elem_nodes[0])[1];
                d_type r2 = grid_.node(elem_nodes[1])[0];   d_type z2 = grid_.node(elem_nodes[1])[1];
                d_type r3 = grid_.node(elem_nodes[2])[0];   d_type z3 = grid_.node(elem_nodes[2])[1];

                d_type area = (r2*z3 - r3*z2) - (r1*z3-r3*z1) + (r1*z2-z1*r2) ;
                H[0] = H[0] + r1 * R_u(o2_.concentration(elem_nodes[0]), co2_.concentration(elem_nodes[2])) * area / 6. ;
                H[1] = H[1] + r2 * R_u(o2_.concentration(elem_nodes[1]), co2_.concentration(elem_nodes[2])) * area / 6. ;
                H[2] = H[2] + r3 * R_u(o2_.concentration(elem_nodes[2]), co2_.concentration(elem_nodes[2])) * area / 6. ;

            };
        };

        void f_co2(mat_type & H, d_type const C_u, d_type const C_v){
            for (int t = 1; t < grid_.nb_elements(); t++) {
                std::vector<d_type> elem_nodes = grid_.element(t);

                d_type r1 = grid_.node(elem_nodes[0])[0];   d_type z1 = grid_.node(elem_nodes[0])[1];
                d_type r2 = grid_.node(elem_nodes[1])[0];   d_type z2 = grid_.node(elem_nodes[1])[1];
                d_type r3 = grid_.node(elem_nodes[2])[0];   d_type z3 = grid_.node(elem_nodes[2])[1];

                d_type area = (r2*z3 - r3*z2) - (r1*z3-r3*z1) + (r1*z2-z1*r2) ;
                H[0] = H[0] - r1 * R_v(o2_.concentration(elem_nodes[0]), co2_.concentration(elem_nodes[2])) * area / 6. ;
                H[1] = H[1] - r2 * R_v(o2_.concentration(elem_nodes[1]), co2_.concentration(elem_nodes[2])) * area / 6. ;
                H[2] = H[2] - r3 * R_v(o2_.concentration(elem_nodes[2]), co2_.concentration(elem_nodes[2])) * area / 6. ;

            };
        };

        void J_o2(mat_type J){
            std::cout<<"empty function"<<std::endl;
        };

        void J_co2(mat_type J){
            std::cout<<"empty function"<<std::endl;
        };


    private:
        pear::component<d_type, vec_type> & o2_;
        pear::component<d_type, vec_type> & co2_;
        pear::grid<d_type> & grid_;
        d_type v_mu_, v_mfv_, k_mu_, k_mv_, k_mfu_, r_q_;

    };




}



#endif