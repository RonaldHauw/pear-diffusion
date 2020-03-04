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
#include "eigen/Eigen/Dense"


namespace pear {

    template <typename d_type, typename f_type, typename vec_type, typename mat_type>
    class nlsolver{
    public:


        nlsolver(f_type f)
        :f_(f)
        {
            std::cout<<"Non-linear solver coupled with abstract function."<<std::endl;
        }

        int solve(){

            std::cout<<"pear::nlsolver.solve(): allocating work memory: "<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()<<", "<<f_.size()<<")"<<std::endl;
            std::cout<<"       - vec_type of size  "<<f_.size()<<std::endl;
            mat_type J; J.resize(f_.size(), f_.size());
            vec_type f; f.resize(f_.size(), 1);
            vec_type dx; dx.resize(f_.size(), 1);

            for (int i = 1; i<4; i++) {
                // reset memory to zero
                J.setZero();
                f.setZero();
                // load the Jacobian and rhs
                f_.J(J);
                f_.f(f);
                // solve the linear system
                f_.cons() = f_.cons() - J.fullPivLu().solve(f) ;
            }
            return 1;
        }



    private:
        f_type f_;
    };



}




#endif