#ifndef RDC_HPP
#define RDC_HPP

#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "grid.hpp"
#include "diffusion.hpp"
#include "reaction.hpp"
#include "component.hpp"
#include "eigen/Eigen/Dense"


namespace pear {

    template <typename d_type, typename vec_type, typename mat_type>
    class rdc{
    public:

        /* Constructor for the equation
         *
         * ! Warning: extremely specific to the equation !
         *
         * IN:  - Two #diffusion# describing the diffusion model of each of the components O2 and C02
         *      - A #reaction# describing the reaction model between the two components 02 and CO2
         * OUT : the private variables corresponding to the above quantities
         *
         */
        rdc(pear::diffusion<d_type, vec_type, mat_type> diff_o2,
                pear::diffusion<d_type, vec_type, mat_type> diff_co2,
                pear::respiration<d_type, vec_type, mat_type> resp)
                : diff_o2_(diff_o2)
                , diff_co2_(diff_co2)
                , resp_(resp)
                {
                }

        /* Construction of the function
        *
        * IN: - A #vector# describing the evaluation of the (non-linear) equation at each gridpoint
        *     - A #sparse_matrix# to be used as intermediate workspace in the algebraic computations
        * OUT: void. Upon completion, the #vector# has been re-initialised (as the memory is recycled) and
        *              contains the diffusion of each of the components as well as the existing reactions between them
        *
        */
        void f(vec_type & x, mat_type & workmat){


            x.setZero();
            diff_o2_.setSparsityPattern(workmat);


            diff_o2_.f(x.segment(diff_o2_.cons_start(), diff_o2_.nb_nodes()));
            diff_co2_.f(x.segment(diff_co2_.cons_start(), diff_co2_.nb_nodes()));
            diff_o2_.J(workmat);
            diff_co2_.J(workmat);

            x = workmat*diff_o2_.cons_full() - x;

            resp_.f(x.segment(diff_o2_.cons_start(), diff_o2_.nb_nodes()),
                    x.segment(diff_co2_.cons_start(), diff_co2_.nb_nodes()));



        };

        /* Construction of the function
        *
        * IN: A #vector# describing the evaluation of the (!) reaction (!) equation at each gridpoint
        * OUT: void. Upon completion, the #vector# has been re-initialised (as the memory is recycled) and
        *              contains the existing reactions between components, ignoring the diffusion model
        *
        */
        void f_react_only(vec_type & x){

            x.setZero();
            resp_.f(x.segment(diff_o2_.cons_start(), diff_o2_.nb_nodes()),
                    x.segment(diff_co2_.cons_start(), diff_co2_.nb_nodes()));
        };


        /* Construction of the Jacobian
        *
        * IN: A #sparse_matrix# describing the Jacobian of the (non-linear) equation
        * OUT: void. Upon completion, the #sparse_matrix# has been re-initialised (as the memory is recycled) and
        *              contains the diffusion of each of the components as well as the existing reactions between them
        *
        */
        void J(mat_type & Jmat){
            diff_o2_.setSparsityPattern(Jmat);

            diff_o2_.J(Jmat);
            diff_co2_.J(Jmat);
            resp_.J(Jmat);
        };

        void suppress_nonlinearity(d_type alpha){   resp_.suppress_nonlinearity(alpha);}
        Eigen::Ref<vec_type> cons(){                return diff_o2_.cons_full();}
        int size(){                                 return diff_o2_.nb_nodes()+diff_co2_.nb_nodes();}


    private:
        pear::diffusion<d_type, vec_type, mat_type> diff_o2_;
        pear::diffusion<d_type, vec_type, mat_type> diff_co2_;
        pear::respiration<d_type, vec_type, mat_type> resp_;

    };

}

#endif