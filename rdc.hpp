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
        rdc(    pear::diffusion<d_type, vec_type> diff_o2,
                pear::diffusion<d_type, vec_type> diff_co2,
                pear::respiration<d_type, vec_type> resp
                )
                : diff_o2_(diff_o2)
                , diff_co2_(diff_co2)
                , resp_(resp)
        {std::cout<<"Reaction diffusion equation composed."<<std::endl;}


        void get_cons(vec_type & x){
            x <<  diff_o2_.get_cons(),  diff_co2_.get_cons(); ; // concatanation
        };


        void f(vec_type & x){
            //
        };

        void J(mat_type & J){
            //
        };

        void r(vec_type & x){}




    private:
        pear::diffusion<d_type, vec_type> diff_o2_;
        pear::diffusion<d_type, vec_type> diff_co2_;
        pear::respiration<d_type, vec_type> resp_;
    };



}




#endif