#ifndef PEAR_DIFFUSION_RESPIRATION_HPP
#define PEAR_DIFFUSION_RESPIRATION_HPP


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

        /* Constructor for respiration
         *
         * IN :    -Two #component#'s describing the components to react
         *         -A #grid# describing the mesh of the problem
         *         -A #vector# (of length 6X1) describing the different physical parameters necessary in the reaction
         *                 model; v_mu, v_mfv, k_mu, k_mv, k_mfu, r_q
         * OUT :   the private variables  which correspond to the above mentioned quantities
         *          - a start value alpha of 1, to suppress non-linearity
         */
        respiration(std::vector<d_type> & param)
                : o2_(o2)
                , co2_(co2)
                , grid_(grid)
                , v_mu_(param[0])
                , v_mfv_(param[1])
                , k_mu_(param[2])
                , k_mv_(param[3])
                , k_mfu_(param[4])
                , r_q_(param[5])
                , alpha_(1.)
        {}

        // Respiration dynamics: functions
        d_type R_u(d_type C_u, d_type C_v){
            return alpha_*v_mu_*C_u/(k_mu_+C_u)/(1+C_v/k_mv_);};
        d_type R_v(d_type C_u, d_type C_v){
            return r_q_*R_u(C_u, C_v) + alpha_*v_mfv_/(1+C_u/k_mfu_);};

        // Respiration dynamics: derivatives
        d_type dR_u_u(d_type C_u, d_type C_v){
            return alpha_*v_mu_ / (k_mu_ + C_u) / (1. + C_v/k_mv_) * (1. - C_u/(k_mu_+C_u));};
        d_type dR_u_v(d_type C_u, d_type C_v){
            return -1 / k_mv_ * alpha_*v_mu_ * C_u / (k_mu_ + C_u) / (1 + C_v/k_mv_) / (1 + C_v/k_mv_);};
        d_type dR_v_u(d_type C_u, d_type C_v){
            return r_q_*dR_u_u(C_u, C_v) - 1/k_mfu_* alpha_*v_mfv_ /(1+C_u/k_mfu_) /(1+C_u/k_mfu_);};
        d_type dR_v_v(d_type C_u, d_type C_v){
            return r_q_*dR_u_v(C_u, C_v);};


        void suppress_nonlinearity(d_type alpha){
            alpha_ = alpha;
        };

    private:
        d_type v_mu_, v_mfv_, k_mu_, k_mv_, k_mfu_, r_q_, alpha_;
    };




}












#endif //PEAR_DIFFUSION_RESPIRATION_HPP
