#ifndef COMPONENTS_HPP
#define COMPONENTS_HPP

#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "grid.hpp"
#include "eigen-3.3.7/Eigen/Dense"



namespace pear {

    template <typename d_type, typename vec_type>
    class component{
    public:

        component(std::string name, pear::grid<d_type> & grid, vec_type & c)
                : name_(name)
                , concentration_(c)
                , grid_(grid)
        {
            std::cout<<"Initialised component: "<<name_<<std::endl;
        }

        std::string name(){
            return name_;
        }

        void set_concentration(Eigen::Matrix<d_type, Eigen::Dynamic, 1> & c){
            concentration_ = c;
        };

        Eigen::Matrix<d_type, Eigen::Dynamic, 1> get_concentration(){
            return concentration_;
        }


        int length(){
            return grid_.length();
        };

    private:
        std::string name_;
        vec_type & concentration_;
        pear::grid<d_type> & grid_;

    };



}




#endif