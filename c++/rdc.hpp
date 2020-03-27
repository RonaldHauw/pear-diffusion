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

        // dedicated constructor for the pear problem
        rdc(pear::diffusion<d_type, vec_type, mat_type> diff_o2,
                pear::diffusion<d_type, vec_type, mat_type> diff_co2,
                pear::respiration<d_type, vec_type, mat_type> resp)
                : diff_o2_(diff_o2)
                , diff_co2_(diff_co2)
                , resp_(resp)
                {
                    std::cout<<"Reaction diffusion equation composed."<<std::endl;
                }

        void f(vec_type & x, mat_type & workmat){

            x.setZero();

            diff_o2_.f(x.segment(diff_o2_.cons_start(), diff_o2_.nb_nodes()), workmat);
            diff_co2_.f(x.segment(diff_co2_.cons_start(), diff_co2_.nb_nodes()), workmat);
            resp_.f(x.segment(diff_o2_.cons_start(), diff_o2_.nb_nodes()),
                    x.segment(diff_co2_.cons_start(), diff_co2_.nb_nodes()));
        };

        void f_react_only(vec_type & x, mat_type & workmat){

            x.setZero();
            resp_.f(x.segment(diff_o2_.cons_start(), diff_o2_.nb_nodes()),
                    x.segment(diff_co2_.cons_start(), diff_co2_.nb_nodes()));
        };

        void J_diff_only(mat_type & Jmat){
            Jmat.setZero();
            diff_o2_.J(Jmat.block(diff_o2_.cons_start(), diff_o2_.cons_start(), diff_o2_.nb_nodes(), diff_o2_.nb_nodes()));
            diff_co2_.J(Jmat.block(diff_co2_.cons_start(), diff_co2_.cons_start(), diff_co2_.nb_nodes(), diff_co2_.nb_nodes()));
        };

        void J_react_only(mat_type & Jmat){
            Jmat.setZero();
            resp_.J(Jmat);
        };



        void J(mat_type & Jmat){
            Jmat.setZero();
            diff_o2_.J(Jmat.block(diff_o2_.cons_start(), diff_o2_.cons_start(), diff_o2_.nb_nodes(), diff_o2_.nb_nodes()));
            diff_co2_.J(Jmat.block(diff_co2_.cons_start(), diff_co2_.cons_start(), diff_co2_.nb_nodes(), diff_co2_.nb_nodes()));
            resp_.J(Jmat);
        };

        void suppress_nonlinearity(d_type alpha){
            resp_.suppress_nonlinearity(alpha);
        }


        Eigen::Ref<vec_type> cons(){
            return diff_o2_.cons_full();
        }

        int size(){
            return diff_o2_.nb_nodes()+diff_co2_.nb_nodes();
        }


    private:
        pear::diffusion<d_type, vec_type, mat_type> diff_o2_;
        pear::diffusion<d_type, vec_type, mat_type> diff_co2_;
        pear::respiration<d_type, vec_type, mat_type> resp_;
    };



}




#endif