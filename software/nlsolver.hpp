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
        int solve(int maxit = 10, d_type steplength = 1., d_type res_pred = 1e-14, d_type res_new = 5e-16){

            std::cout<<"        Working memory allocated: "<<f_.size()*2*16+f_.size()*3<<" elements of d_type:"<<std::endl;
            std::cout<<"           - 2 times a mat_type of size ("<<f_.size()<<", "<<f_.size()<<")"<<std::endl;
            std::cout<<"             of which 16 elements per row are preallocated;"<<std::endl;
            std::cout<<"           - 3 times a vec_type of size "<<f_.size()<<";"<<std::endl;

            std::cout<<"// START SOLVE //"<<std::endl;

            mat_type J(f_.size(), f_.size());
            mat_type J_work(f_.size(), f_.size());

            J.reserve(f_.size()*16);
            J_work.reserve(f_.size()*16);


            grid_.setSparsityPattern(J);
            grid_.setSparsityPattern(J_work);


            vec_type workvec; workvec.resize(f_.size(), 1);
            vec_type direction; direction.resize(f_.size(), 1);
            vec_type residual; residual.resize(f_.size(), 1);

            d_type cur_alpha = 0.0;
            f_.suppress_nonlinearity(cur_alpha);

            d_type res = 1.0;
            Eigen::SparseLU<mat_type> linsolver;
            f_.J(J);
            linsolver.analyzePattern(J);

            int i = 0;

            while (cur_alpha < 1.){
                i++;

                // Direction of prediction
                f_.f_react_only(residual);
                f_.J(J);
                linsolver.factorize(J);
                direction = linsolver.solve(residual);


                // Backtracking on prediction step length
                workvec = f_.cons();
                steplength = 1.-cur_alpha;
                for (int k = 0; k < 6; k++){
                    f_.cons() = workvec-steplength*direction;           // take a test step
                    f_.suppress_nonlinearity(cur_alpha+steplength);
                    f_.f(residual, J_work);                             // residual
                    if (residual.norm()/f_.cons().norm()>res_pred) {    // if residual is small enough, keep step
                        steplength *= 0.5;
                    } else {
                        break;
                    }
                }


                cur_alpha += steplength;


                for (int j = 1; j < maxit; j++ ){

                    display_progress( cur_alpha,  residual.norm()/f_.cons().norm(),  res_new, res_pred);

                    f_.f(residual, J_work);
                    f_.J(J);

                    linsolver.factorize(J);
                    direction = linsolver.solve(residual);  // direction

                    res = residual.norm()/f_.cons().norm();

                    steplength = 1.;
                    workvec = f_.cons(); // store current concentrations

                    // backtracking linesearch for residual decrease
                    for (int k = 0; k < 6; k++){
                        f_.cons() = workvec-steplength*direction;
                        f_.f(residual, J_work); // residual
                        if (residual.norm()/f_.cons().norm()>res) {
                            steplength *= 0.5;
                        } else {
                            break;
                        }
                    }

                    res = residual.norm()/f_.cons().norm();

                    if (res < res_new) {
                        display_progress( cur_alpha,  residual.norm()/f_.cons().norm(),  res_new, res_pred);

                        break;
                    }
                }
            }
            std::cout<<std::endl;
            return 1;
        }


        void display_progress(d_type cur_alpha, d_type newton_res, d_type res_new, d_type res_pred){
            int bar_width = 20;
            d_type hom_progress = cur_alpha;
            d_type newton_progress = (-log(res_pred)+log(newton_res))/(-log(res_pred)+log(res_new));
            //std::cout<<"\\u1b[1F"<<"\\u1B[0K"<<"\rhomotopy: [";
            std::cout<<"\rhomotopy: [";

            d_type hom_pos = bar_width * hom_progress;
            for (int k = 0; k < bar_width; k++){
                if (k<hom_pos) {
                    std::cout << "=";
                } else {
                    std::cout << " ";
                }
            }
            std::cout<<"] "<< round(hom_progress*100) << "%     newton: [";

            d_type new_pos = bar_width * newton_progress;
            for (int k = 0; k < bar_width; k++){
                if (k<new_pos) {
                    std::cout << "=";
                } else {
                    std::cout << " ";
                }
            }
            int newton_perc = 100;
            if (round(newton_progress*100) < 100){newton_perc = round(newton_progress*100);};
            std::cout<<"] "<< newton_perc<<"%";
            fflush(stdout);
        }

    private:
        f_type f_;
        pear::grid<d_type, mat_type> & grid_;
    };

}




#endif
