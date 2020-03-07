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

        int solve(int maxit = 10, d_type steplength = 1.){

            std::cout<<"pear::nlsolver.solve(): allocating work memory: "<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()<<", "<<f_.size()<<")"<<std::endl;
            std::cout<<"       - vec_type of size  "<<f_.size()<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()/2<<", "<<f_.size()/2<<")"<<std::endl;
            mat_type workmat;
            workmat.resize(f_.size()/2, f_.size()/2);
            mat_type J; J.resize(f_.size(), f_.size());
            vec_type f; f.resize(f_.size(), 1);

            for (int i = 1; i<maxit; i++) {
                // load the Jacobian and rhs
                f_.J(J);
                f_.f(f, workmat);
                std::cout<<"iterations = "<<i<<"  residual = "<<f.norm()<<std::endl;

                // solve the linear system
                f_.cons() = f_.cons() - steplength*J.fullPivLu().solve(f) ; // shortened step length
            }
            return 1;
        }

// observation: CO2 keeps increasing in solution after O2 reaches zero
// step length stabilizes
// corresponds to Matlab
// parsing
//
    private:
        f_type f_;
    };



}




#endif