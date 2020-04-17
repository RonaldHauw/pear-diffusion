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
#include "eigen/Eigen/IterativeLinearSolvers"

// todo: still two hyperparameters in the code (! )
// todo: check dubbele randpunten

namespace pear {

    template <typename d_type, typename f_type, typename vec_type, typename mat_type>
    class nlsolver{
    public:

        /* Contructor for the non-linear solver
         *
         * IN:  - A #rdc_equation# describing the non-linear equation to be solved
         *      - A #grid# describing the mesh of the problem
         * OUT :   the private variables  which correspond to the above mentioned quantities
         */
        nlsolver(f_type f, pear::grid<d_type, mat_type> & grid)
        :f_(f)
        , grid_(grid)
        {
            std::cout<<"//NON-LINEAR//"<<std::endl;
        }

        /* Solutions the non-linear equation
         *
         * IN:      - A #int# describing the maximal number of iterations to be done. Default is 10.
         *          - A #double# describing the step length to be taken. Default is 1.
         *          - A #double# indicating the start value of alpha. Default is 0.
         * OUT :    A #int#=1 describing the . Upon completion, the solutions are exported in a .txt file.
         *              The residuals have been printed.
         *
         */
        int solve(int maxit = 10, d_type steplength = 1., d_type res_pred = 1e-12, d_type res_new = 1e-19){

            std::cout<<"        Working memory allocated: "<<std::endl;
            std::cout<<"            mat_type of size ("<<f_.size()<<", "<<f_.size()<<")"<<std::endl;
            std::cout<<"            vec_type of size  "<<f_.size()<<std::endl;
            std::cout<<"            mat_type of size ("<<f_.size()/2<<", "<<f_.size()/2<<")"<<std::endl;

            // Initialisation system
            mat_type J(f_.size(), f_.size());
            mat_type J_work(f_.size(), f_.size());

            J.reserve(f_.size()*16);        // ! Presumption! Highly unlikely that a node has more than 16 neighbours
            J_work.reserve(f_.size()*16);

            grid_.setSparsityPattern(J);
            grid_.setSparsityPattern(J_work);

            vec_type direction; direction.resize(f_.size(), 1);
            vec_type residual; residual.resize(f_.size(), 1);
            vec_type workvec; workvec.resize(f_.size(), 1);

            //Initialisation homotopy continuation
            d_type cur_alpha = 0.0;
            f_.suppress_nonlinearity(cur_alpha);

            // Initialisation non-linear solver
            d_type res = 1.0;
            Eigen::ConjugateGradient<mat_type> linsolver;
            linsolver.analyzePattern(J);
            int i = 0;

            std::cout<<"        Homotopy progress --> "<<100*(cur_alpha)<<"%"<<std::endl;
            // Homotopy continuation
            while (cur_alpha < 1.){
                i++;

                // Direction of prediction
                f_.f_react_only(residual);
                f_.J(J);
                linsolver.factorize(J);
                direction = linsolver.solve(residual);
                workvec = f_.cons();

                // Backtracking on prediction step length
                steplength = 1.-cur_alpha;
                for (int k = 0; k < 6; k++){
                    f_.cons() = workvec-steplength*direction;        // Try a test step...
                    f_.suppress_nonlinearity(cur_alpha+steplength);
                    f_.f(residual, J_work);                          // ... compute the residual ...
                    if (residual.norm()>res_pred) {                     // ... if the residual remains too large, halve the steplength
                        std::cout<<"        prediction residual = "<<residual.norm()<<"  steplength = "<<steplength<<std::endl;
                        steplength *= 0.5;
                    } else {
                        break;
                    }
                }
                cur_alpha += steplength;
                std::cout<<"                          --> "<<round(100*cur_alpha)<<"%"<<std::endl;

                for (int j = 1; j < maxit; j++ ){

                    f_.f(residual, J_work);
                    f_.J(J);
                    res = residual.norm();
                    linsolver.factorize(J);
                    direction = linsolver.solve(residual);

                    steplength = 1.;
                    workvec = f_.cons();
                    // Backtracking on line-search step length
                    for (int k = 0; k < 100; k++){
                        f_.cons() = workvec-steplength*direction;
                        f_.f(residual, J_work);
                        if (residual.norm()>res) {
                            std::cout<<"        newton residual = "<<residual.norm()<<"    steplength = "<<steplength<<std::endl;
                            steplength *= 0.5;
                        } else {
                            break;
                        }
                    }

                    if (residual.norm() < res_new) {
                        std::cout<<"        with "<<i<<" Newton iterations until a residual norm of "<<residual.norm()<<std::endl;
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