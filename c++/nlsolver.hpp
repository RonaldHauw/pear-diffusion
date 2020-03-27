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

// todo: still two hyperparameters in the code (! )
// todo: check dubbele randpunten

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
            vec_type f3; f3.resize(f_.size(), 1);

            d_type cur_alpha = 0.0;
            f_.suppress_nonlinearity(cur_alpha);

            d_type res = 1.0;
            int i = 0;
            while (cur_alpha < 1.){
                i++;


                // prediction step
                f_.f_react_only(f2, workmat);  // f stores H
                f_.J(J);
                f = J.fullPivLu().solve(-f2); // search direction

                f2 = f_.cons(); // save concentrations

                // backtracking on prediction step length
                steplength = 1.-cur_alpha;
                for (int k = 0; k < 100; k++){
                    f_.cons() = f2+steplength*f; // take a test step
                    f_.suppress_nonlinearity(cur_alpha+steplength);
                    f_.f(f3, workmat); // residual
                    if (f3.norm()>1e-11) { // if residual is small enough, keep step
                        steplength *= 0.5;
                    } else {
                        std::cout<<"prediction steplength = "<<steplength<<std::endl;
                        break;
                    }
                }
                cur_alpha += steplength;

                // newton solver
                for (int j = 1; j < maxit; j++ ){

                    f_.f(f, workmat);
                    f_.J(J);
                    f2 = J.fullPivLu().solve(f);  // direction

                    res = f.norm();

                    steplength = 1.;
                    f = f_.cons(); // store current concentrations

                    // backtracking linesearch for residual decrease
                    for (int k = 0; k < 100; k++){
                        f_.cons() = f-steplength*f2;
                        f_.f(f3, workmat); // residual
                        if (f3.norm()>res) {
                            steplength *= 0.5;
                        } else {
                            std::cout<<"steplength = "<<steplength<<std::endl;
                            break;
                        }
                    }


                    res = f3.norm();

                    std::cout<<"iterations = "<<i<<"  newton residual = "<<res<<std::endl;
                    if (res < 5e-19) {
                        break;
                    }
                }
            }
            return 1;
        }



    private:
        f_type f_;
    };



}




#endif