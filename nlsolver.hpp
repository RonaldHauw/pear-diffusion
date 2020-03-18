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
#include "eigen-3.3.7/Eigen/Dense"


namespace pear {

    template <typename d_type, typename f_type>
    class nlsolver{
    public:


        nlsolver(f_type f)
        :f_(f)
        //,lasolve_(lasolver)
        {
            std::cout<<"Non-linear solver coupled with abstract function."<<std::endl;
        }


        int solve(){
            for (int i = 1; i<10; i++) {
                Eigen::Matrix<d_type, Eigen::Dynamic, Eigen::Dynamic> J = f_.jacobian();
                Eigen::Matrix<d_type, Eigen::Dynamic, 1> fx = f_.eval();
                Eigen::Matrix<d_type, Eigen::Dynamic, 1> x = J.colPivHouseholderQr().solve(fx);
                f_.set_concentrations(x);
            }
            return 1;
        }



    private:
        f_type f_;
        //lasolve_type lasolve_;
    };



}




#endif