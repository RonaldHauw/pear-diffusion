#ifndef COMPONENTS_HPP
#define COMPONENTS_HPP
// hogere orde benadering driehoekjes
// reeksontwikkeling niet lineariteit in basis functies (termen vallen weg)
// discretisatie step te groot (convergentie plot)
// minteken verkeerd ergens
// boundary integraal juist?
// expressie templates eigen aan/uit?

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

        component(std::string name, pear::grid<d_type> & grid, vec_type & c, int start, int stop, int stride)
                : name_(name)
                , concentration_(c)
                , grid_(grid)
                , start_(start)
                , stop_(stop)
                , stride_(stride)
        {
            std::cout<<"Initialised component: "<<name_<<std::endl;
        }

        std::string name(){
            return name_;
        }



        vec_type & get_cons(){
            return concentration_(seq(start_, stop_, stride_), 1);
        }

        d_type & concentration(int i){
            return concentration_(start_ + i*stride_, 1);
        }

        int length(){
            return grid_.length();
        };

    private:
        std::string name_;
        vec_type & concentration_;
        pear::grid<d_type> & grid_;
        int start_;
        int stop_;
        int stride_;

    };



}




#endif