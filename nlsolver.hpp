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
#include "eigen/Eigen/SparseCore"
#include "eigen/Eigen/SparseLU"


namespace pear {

    template <typename d_type, typename f_type, typename vec_type, typename mat_type>
    class nlsolver{
    public:


        nlsolver(f_type f, pear::grid<d_type, mat_type> & grid)
        :f_(f)
        , grid_(grid)
        {
            std::cout<<"Non-linear solver coupled with abstract function."<<std::endl;
        }

        int solve(int maxit = 10, d_type steplength = 1., d_type alpha = 0.){

            std::cout<<"pear::nlsolver.solve(): allocating work memory: "<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()<<", "<<f_.size()<<")"<<std::endl;
            std::cout<<"       - vec_type of size  "<<f_.size()<<std::endl;
            std::cout<<"       - mat_type of size ("<<f_.size()/2<<", "<<f_.size()/2<<")"<<std::endl;

            mat_type workmat(f_.size()/2, f_.size()/2);
            mat_type J(f_.size(), f_.size());
            mat_type J2(f_.size(), f_.size());

            grid_.setSparsityPattern(J);
            grid_.setSparsityPattern(J2);

            vec_type f; f.resize(f_.size(), 1);
            vec_type f2; f2.resize(f_.size(), 1);
            vec_type f3; f3.resize(f_.size(), 1);

            d_type cur_alpha = 0.0;
            f_.suppress_nonlinearity(cur_alpha);

            d_type res = 1.0;
            Eigen::SparseLU<Eigen::SparseMatrix<double> > linsolver;

            for (int i = 0; i < 1./alpha; i++){

                // prediction step
                f_.f_react_only(f2);  // f stores H
                f_.J_diff_only(J); // J stores K
                std::cout<<"J done"<<std::endl;
                f_.J_react_only(J2); // J2 stores dHdc
                std::cout<<"J2 done"<<std::endl;
                J = J+cur_alpha*J2; // J stores dGammadc
                J.makeCompressed();

                if (i == 0){ linsolver.analyzePattern(J); }
                linsolver.compute(J);
                f_.cons() = f_.cons() + alpha*linsolver.solve(f2);
                cur_alpha += alpha; // take a step
                f_.suppress_nonlinearity(cur_alpha);

                for (int j = 1; j < maxit; j++ ){

                    f_.f(f, workmat);
                    f_.J(J);
                    J.makeCompressed();

                    if (i == 0){ linsolver.analyzePattern(J); }
                    linsolver.compute(J);
                    f2 = linsolver.solve(f2);  // direction

                    res = f.norm();

                    steplength = 1.;
                    f = f_.cons(); // store current concentrations

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
        pear::grid<d_type, mat_type> & grid_;
    };

}




#endif