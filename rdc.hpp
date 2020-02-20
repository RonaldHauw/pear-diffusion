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

    template <typename d_type, typename vec_type>
    class rdc{
    public:

        // dedicated constructor for the pear problem
        rdc(pear::component<d_type, vec_type> & o2, pear::component<d_type, vec_type> & co2, pear::diffusion<d_type, vec_type> diff_o2, pear::diffusion<d_type, vec_type> diff_co2,
        pear::respiration_o2<d_type, vec_type> react_o2, pear::respiration_co2<d_type, vec_type> react_co2)
                : o2_(o2)
                , co2_(co2)
                , diff_o2_(diff_o2)
                , diff_co2_(diff_co2)
                , react_o2_(react_o2)
                , react_co2_(react_co2)
        {
            std::cout<<"Reaction diffusion equation composed."<<std::endl;
        }


        Eigen::Matrix<d_type, Eigen::Dynamic, 1> eval(){

            return o2_.get_concentration();
            //return Eigen::Matrix<d_type, Eigen::Dynamic, 1>::Random(o2_.length(), 1);
        };


        Eigen::Matrix<d_type, 6, 6> jacobian(){
            Eigen::Matrix<d_type, 6, 6> test;
            return test;
        };

        void set_concentrations(Eigen::Matrix<d_type, Eigen::Dynamic, 1> x){
            o2_.set_concentration(x);
        };


    private:
        pear::component<d_type, vec_type> & o2_;
        pear::component<d_type, vec_type> & co2_;
        pear::respiration_o2<d_type, vec_type> react_o2_;
        pear::respiration_co2<d_type, vec_type> react_co2_;
        pear::diffusion<d_type, vec_type> diff_o2_;
        pear::diffusion<d_type, vec_type> diff_co2_;
    };



}




#endif