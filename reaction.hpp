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

        // Respiration dynamics: functions

        d_type R_u(d_type C_u, d_type C_v){
            return v_mu_*C_u/(k_mu_+C_u)/(1+C_v/k_mv_);
        };

        d_type R_v(d_type C_u, d_type C_v){
            return r_q_*R_u(C_u, C_v) + v_mfv_/(1+C_u/k_mfu_);
        };

        // Respiration dynamics: derivatives

        d_type dR_u_u(d_type C_u, d_type C_v){
            //std::cout<<" v_mu = "<<v_mu_<<" k_mu_ = "<<k_mu_<<" k_mv_ = "<<k_mv_;
            return v_mu_ / (k_mu_ + C_u) / (1. + C_v/k_mv_) * (1. - C_u/(k_mu_+C_u));
        };

        d_type dR_u_v(d_type C_u, d_type C_v){
            return -1 / k_mv_ * v_mu_ * C_u / (k_mu_ + C_u) / (1 + C_v/k_mv_) / (1 + C_v/k_mv_);
        };

        d_type dR_v_u(d_type C_u, d_type C_v){
            return r_q_*dR_u_u(C_u, C_v) - 1/k_mfu_* v_mfv_ /(1+C_u/k_mfu_) /(1+C_u/k_mfu_);
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

        void J(Eigen::Ref<mat_type> dH) {
            for (int t = 1; t < grid_.nb_nodes()+1; t++) {
                std::vector<int> elements = grid_.elements_for_node(t);
                std::cout<<" elements = ["<<elements.size()<<"]"<<std::endl;
                std::cout<<" node = "<<t<<std::endl;
                d_type area = 0.;
                for (int e = 1; e < elements.size()+1; e++) {

                        std::vector<int> nodes = grid_.element(elements[e - 1]);

                        d_type r1 = grid_.node(nodes[0])[0]; d_type z1 = grid_.node(nodes[0])[1];
                        d_type r2 = grid_.node(nodes[1])[0]; d_type z2 = grid_.node(nodes[1])[1];
                        d_type r3 = grid_.node(nodes[2])[0]; d_type z3 = grid_.node(nodes[2])[1];
                        area += (r2 * z3 - r3 * z2) - (r1 * z3 - r3 * z1) + (r1 * z2 - z1 * r2);
                        //std::cout<<"r1 = "<<r1<<"  r2 = "<<r2<<" r3 = "<<r3;
                        //std::cout<<" z1  = "<<z1<<" z2 = "<<z2<<" z3 = "<<z3<<std::endl;
                        //std::cout<<"area it "<<area<<std::endl;
                };

                dH(t-1, t-1)                                            += area * grid_.node(t)[0] * dR_u_u(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                std::cout<<"area "<<area<<"  druu "<<dR_u_u(o2_.concentration(t), co2_.concentration(t) )<<"  r(i) "<<grid_.node(t)[0]<<std::endl;
                //std::cout<<"dvar_1 "<<area * grid_.node(t)[0] * dR_u_u(o2_.concentration(t), co2_.concentration(t) ) / 6.<<" "<<t<<" "<<t<<std::endl;
                std::cout<<"Cco2 = "<<o2_.concentration(t)<<" Co2 = "<<co2_.concentration(t)<<std::endl;

                dH(t-1, t + grid_.nb_nodes() -1)                        += area * grid_.node(t)[0] * dR_u_v(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH(t + grid_.nb_nodes() -1, t-1)                        -= area * grid_.node(t)[0] * dR_v_u(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH(t + grid_.nb_nodes() -1, t + grid_.nb_nodes() -1)    -= area * grid_.node(t)[0] * dR_v_v(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
            };
        };

    private:
        pear::component<d_type, vec_type> & o2_;
        pear::component<d_type, vec_type> & co2_;
        pear::grid<d_type> & grid_;
        d_type v_mu_, v_mfv_, k_mu_, k_mv_, k_mfu_, r_q_;
    };




}



#endif