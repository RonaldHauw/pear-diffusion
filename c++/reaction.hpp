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
#include "eigen/Eigen/Dense"
#include "component.hpp"
#include "diffusion.hpp"


namespace pear {


    template <typename d_type, typename vec_type, typename mat_type>
    class respiration{
    public:

        respiration(pear::component<d_type, vec_type, mat_type> & co2,
                pear::component<d_type, vec_type, mat_type> & o2,
                pear::grid<d_type, mat_type> & grid,
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
                , alpha_(1.) // To suppress non-linearity
        {
            std::cout<<"Respiration between O2 and CO2"<<std::endl;
        }

        // Respiration dynamics: functions

        d_type R_u(d_type C_u, d_type C_v){
            return alpha_*v_mu_*C_u/(k_mu_+C_u)/(1+C_v/k_mv_);
        };

        d_type R_v(d_type C_u, d_type C_v){
            return r_q_*R_u(C_u, C_v) + alpha_*v_mfv_/(1+C_u/k_mfu_);
        };

        // Respiration dynamics: derivatives

        d_type dR_u_u(d_type C_u, d_type C_v){
            return alpha_*v_mu_ / (k_mu_ + C_u) / (1. + C_v/k_mv_) * (1. - C_u/(k_mu_+C_u));
        };

        d_type dR_u_v(d_type C_u, d_type C_v){
            return -1 / k_mv_ * alpha_*v_mu_ * C_u / (k_mu_ + C_u) / (1 + C_v/k_mv_) / (1 + C_v/k_mv_);
        };

        d_type dR_v_u(d_type C_u, d_type C_v){
            return r_q_*dR_u_u(C_u, C_v) - 1/k_mfu_* alpha_*v_mfv_ /(1+C_u/k_mfu_) /(1+C_u/k_mfu_);
        };

        d_type dR_v_v(d_type C_u, d_type C_v){
            return r_q_*dR_u_v(C_u, C_v);
        };

        void f(Eigen::Ref<vec_type>  H_o2, Eigen::Ref<vec_type>  H_co2) {
            for (int t = 1; t < grid_.nb_elements()+1; t++) {

                std::vector<int> elem_nodes = grid_.element(t);

                d_type r1 = grid_.node(elem_nodes[0])[0];   d_type z1 = grid_.node(elem_nodes[0])[1];
                d_type r2 = grid_.node(elem_nodes[1])[0];   d_type z2 = grid_.node(elem_nodes[1])[1];
                d_type r3 = grid_.node(elem_nodes[2])[0];   d_type z3 = grid_.node(elem_nodes[2])[1];

                d_type area = (r2*z3 - r3*z2) - (r1*z3-r3*z1) + (r1*z2-z1*r2) ;

                H_o2(elem_nodes[0]-1) += r1 * R_u(o2_.concentration(elem_nodes[0]), co2_.concentration(elem_nodes[0])) * area / 6. ;
                H_o2(elem_nodes[1]-1) += r2 * R_u(o2_.concentration(elem_nodes[1]), co2_.concentration(elem_nodes[1])) * area / 6. ;
                H_o2(elem_nodes[2]-1) += r3 * R_u(o2_.concentration(elem_nodes[2]), co2_.concentration(elem_nodes[2])) * area / 6. ;


                H_co2(elem_nodes[0]-1) -= r1 * R_v(o2_.concentration(elem_nodes[0]), co2_.concentration(elem_nodes[0])) * area / 6. ;
                H_co2(elem_nodes[1]-1) -= r2 * R_v(o2_.concentration(elem_nodes[1]), co2_.concentration(elem_nodes[1])) * area / 6. ;
                H_co2(elem_nodes[2]-1) -= r3 * R_v(o2_.concentration(elem_nodes[2]), co2_.concentration(elem_nodes[2])) * area / 6. ;

            };
        };

        void J(mat_type & dH) { // checked with Matlab !
            for (int t = 1; t < grid_.nb_nodes()+1; t++) {
                std::vector<int> elements = grid_.elements_for_node(t);
                d_type area = 0.;
                for (int e = 1; e < static_cast<signed int>(elements.size())+1; e++) {

                        std::vector<int> nodes = grid_.element(elements[e - 1]);

                        d_type r1 = grid_.node(nodes[0])[0]; d_type z1 = grid_.node(nodes[0])[1];
                        d_type r2 = grid_.node(nodes[1])[0]; d_type z2 = grid_.node(nodes[1])[1];
                        d_type r3 = grid_.node(nodes[2])[0]; d_type z3 = grid_.node(nodes[2])[1];
                        area += (r2 * z3 - r3 * z2) - (r1 * z3 - r3 * z1) + (r1 * z2 - z1 * r2);

                };

                dH.coeffRef(t-1, t-1)                                            +=  area * grid_.node(t)[0] * dR_u_u(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH.coeffRef(t-1, t + grid_.nb_nodes() -1)                        +=  area * grid_.node(t)[0] * dR_u_v(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH.coeffRef(t + grid_.nb_nodes() -1, t-1)                        += -area * grid_.node(t)[0] * dR_v_u(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH.coeffRef(t + grid_.nb_nodes() -1, t + grid_.nb_nodes() -1)    += -area * grid_.node(t)[0] * dR_v_v(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
            };

        };

        void suppress_nonlinearity(d_type alpha){
            alpha_ = alpha;
            std::cout<<"pear::respiration.suppress_nonlinearity(): alpha_ = "<<alpha_<<std::endl;
        };

    private:
        pear::component<d_type, vec_type, mat_type> & o2_;
        pear::component<d_type, vec_type, mat_type> & co2_;
        pear::grid<d_type, mat_type> & grid_;
        d_type v_mu_, v_mfv_, k_mu_, k_mv_, k_mfu_, r_q_, alpha_;
    };




}



#endif