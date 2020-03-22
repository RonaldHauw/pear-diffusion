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

            mat_type J(f_.size(), f_.size());
            mat_type J2(f_.size(), f_.size());

            grid_.setSparsityPattern(J);
            grid_.setSparsityPattern(J2);

            std::cout<<J<<std::endl;
            std::cout<<J2<<std::endl;

            vec_type f; f.resize(f_.size(), 1);
            vec_type f2; f2.resize(f_.size(), 1);
            vec_type f3; f3.resize(f_.size(), 1);

            d_type cur_alpha = 0.0;
            f_.suppress_nonlinearity(cur_alpha);

            d_type res = 1.0;
            Eigen::SparseLU<Eigen::SparseMatrix<double> > linsolver;

            std::cout<<"nlsolve: checkpoint 1"<<std::endl;

            for (int i = 0; i < 1./alpha; i++){

                // prediction step
                f_.f_react_only(f2);  // f stores H
                std::cout<<"f2 done"<<std::endl;
                f_.J_diff_only(J); // J stores K
                std::cout<<"J done"<<std::endl;
                f_.J_react_only(J2); // J2 stores dHdc
                std::cout<<"J2 done"<<std::endl;
                J = J+cur_alpha*J2; // J stores dGammadc
                J.makeCompressed();

                std::cout<<"nlsolve: checkpoint 2"<<std::endl;

                if (i == 0){ linsolver.analyzePattern(J); }
                linsolver.compute(J);
                f_.cons() = f_.cons() + alpha*linsolver.solve(f2);
                cur_alpha += alpha; // take a step
                f_.suppress_nonlinearity(cur_alpha);

                std::cout<<"nlsolve: checkpoint 3"<<std::endl;

                for (int j = 1; j < maxit; j++ ){

                    grid_.setSparsityPattern(J);
                    grid_.setSparsityPattern(J2);
                    std::cout<<"set sparsity pattern"<<std::endl;
                    f_.f(f, J2);
                    std::cout<<"set f"<<std::endl;
                    f_.J(J);
                    std::cout<<"set J"<<std::endl;
                    J.makeCompressed();

                    if (i == 0){ linsolver.analyzePattern(J); }
                    linsolver.compute(J);
                    f2 = linsolver.solve(f2);  // direction

                    std::cout<<"nlsolve: checkpoint 4"<<std::endl;

                    res = f.norm();

                    steplength = 1.;
                    f = f_.cons(); // store current concentrations

                    std::cout<<"nlsolve: checkpoint 5"<<std::endl;

                    for (int k = 0; k < 100; k++){
                        f_.cons() = f-steplength*f2;
                        f_.f(f3, J2); // residual
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