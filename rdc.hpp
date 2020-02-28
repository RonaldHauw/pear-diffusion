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
#include "eigen-3.3.7/Eigen/Dense"


namespace pear {

    template <typename d_type, typename vec_type, typename mat_type>
    class rdc{
    public:

        // dedicated constructor for the pear problem
        rdc(    pear::diffusion<d_type, vec_type, mat_type> diff_o2,
                pear::diffusion<d_type, vec_type, mat_type> diff_co2,
                pear::respiration<d_type, vec_type, mat_type> resp
                )
                : diff_o2_(diff_o2)
                , diff_co2_(diff_co2)
                , resp_(resp)
        {std::cout<<"Reaction diffusion equation composed."<<std::endl;}

        void f(vec_type & x){
            mat_type K;
            K.resize(diff_o2_.nb_nodes(), diff_o2_.nb_nodes());
            diff_o2_.f(x, K);
        };

        void f(vec_type & x){
            mat_type K;
            K.resize(diff_o2_.nb_nodes(), diff_o2_.nb_nodes());
            diff_o2_.f(x.segment(0, diff_o2_.nb_nodes() ), K );
            diff_co2_.f(x.segment(diff_o2_.nb_nodes()-1, diff_co2_.nb_nodes() ), K );
        };

        void J(mat_type & J){
            J.setZero();
            diff_o2_.J(J);
        };

        void J(mat_type & J){
            J.setZero();
            diff_o2_.J(K.block(0,0,diff_o2_.nb_nodes(),diff_o2_.nb_nodes()) );
            diff_co2_.J(K.block(diff_o2_.nb_nodes()-1,diff_o2_.nb_nodes()-1,diff_co2_.nb_nodes(),diff_co2_.nb_nodes()) );
        };

        void set_cons(vec_type & x){
            diff_o2_.set_cons(x);
        }

        int size(){
            return diff_o2_.nb_nodes();//+diff_co2_.nb_nodes()
        }

    private:
        pear::diffusion<d_type, vec_type, mat_type> diff_o2_;
        pear::diffusion<d_type, vec_type, mat_type> diff_co2_;
        pear::respiration<d_type, vec_type, mat_type> resp_;
    };



}




#endif