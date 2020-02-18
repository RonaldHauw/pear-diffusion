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


namespace pear {

    template <typename d_type>
    class component{
    public:

        component(std::string name, pear::grid<d_type> grid)
                : name_(name)
                , concentration_(std::vector<d_type>(grid.length()))
        {
            std::cout<<"Initialised component: "<<name_<<std::endl;

        }

        std::string name(){
            return name_;
        }

    private:
        std::string name_;
        std::vector<d_type> concentration_;

    };



}




#endif