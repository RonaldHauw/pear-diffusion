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

        int solve(int maxit = 10, d_type steplength = 1., d_type alpha = 0.){

            std::cout<<"pear::nlsolver.solve(): allocating work memory: "<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()<<", "<<f_.size()<<")"<<std::endl;
            std::cout<<"       - vec_type of size  "<<f_.size()<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()/2<<", "<<f_.size()/2<<")"<<std::endl;
            mat_type workmat;
            workmat.resize(f_.size()/2, f_.size()/2);
            mat_type J; J.resize(f_.size(), f_.size());
            mat_type J2; J2.resize(f_.size(), f_.size());
            vec_type f; f.resize(f_.size(), 1);
            vec_type f2; f2.resize(f_.size(), 1);

            d_type cur_alpha = 0.0;
            f_.suppress_nonlinearity(cur_alpha);
            for (int i = 0; i < 1./alpha; i++){
                // prediction step
                f_.f_react_only(f2, workmat);  // f stores H
                f_.J_diff_only(J); // J stores K
                f_.J_react_only(J2); // J2 stores dHdc
                J = J+cur_alpha*J2; // J stores dGammadc
                f_.cons() = f_.cons() + alpha*J.fullPivLu().solve(f2);
                cur_alpha += alpha; // take a step
                f_.suppress_nonlinearity(cur_alpha);

                for (int j = 1; j < maxit; j++ ){
                    if (f.norm()>1e-9){
                        steplength = 0.3;
                    } else if(f.norm()>1e-8) {
                        steplength = 0.2;
                    } else {
                        steplength = 0.99;
                    }
                    f_.f(f, workmat);
                    f_.J(J);
                    f_.cons() = f_.cons() - steplength*J.fullPivLu().solve(f) ; // shortened step length
                    std::cout<<"iterations = "<<i<<"  newton residual = "<<f.norm()<<std::endl;
                    if (f.norm() < 1e-20) {
                        break;
                    }
                }
            }
            return 1;
        }

// observation: CO2 keeps increasing in solution after O2 reaches zero
// step length stabilizes
// corresponds to Matlab
// parsing
//
// boundary conditions violated
// add C_amb?
// dubbele randpunten zijn een issue



    private:
        f_type f_;
    };



}




#endif