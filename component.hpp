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


namespace pear {

    class component{
    public:

        component(std::string name)
                : name_(name)
        {
            std::cout<<"Initialised component: "<<name_<<std::endl;
        }

        std::string name(){
            return name_;
        }

    private:
        std::string name_;

    };



}




#endif