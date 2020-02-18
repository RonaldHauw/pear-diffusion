#ifndef NLSOLVER_HPP
#define NLSOLVER_HPP

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
#include "rdc.hpp"

namespace pear {

    template <typename d_type, typename f_type>
    class nlsolver{
    public:


        nlsolver(f_type f)
        :f_(f)

        {
            std::cout<<"Non-linear solver coupled with abstract function."<<std::endl;
        }



    private:
        f_type f_;
    };



}




#endif