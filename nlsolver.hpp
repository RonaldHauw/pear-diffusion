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

    template <typename d_type, typename f_type, typename vec_type, typename mat_type>
    class nlsolver{
    public:


        nlsolver(f_type f)
        :f_(f)
        {
            std::cout<<"Non-linear solver coupled with abstract function."<<std::endl;
        }

        int solve(){
            mat_type J; J.resize(f_.size(), f_.size());
            vec_type f; f.resize(f_.size(), 1);

            for (int i = 1; i<3; i++) {
                J.setZero(); f_.J(J); f_.f(f);
                std::cout<<J<<std::endl;
                // Jabobian likely singular
                vec_type x; x.resize(f_.size(), 1);
                x = J.fullPivLu().solve(f);
                //std::cout<<"x"<<x<<std::endl;
                std::cout<<" correctness: "<<(J*x-f).norm()<<"  "<<f.norm()<<std::endl;
                f_.set_cons(x);
                // set concentrations doesnt work

                std::cout<<" rows"<<J.rows()<<"   cols "<<J.cols()<<std::endl;
                Eigen::SelfAdjointEigenSolver<mat_type> eigensolver(J);
                if (eigensolver.info() != Eigen::Success) abort();
                std::cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
                std::cout << "Here's a matrix whose columns are eigenvectors of J \n"
                     << "corresponding to these eigenvalues:\n"
                     << eigensolver.eigenvectors() << std::endl;
            }
            return 1;
        }



    private:
        f_type f_;
        //lasolve_type lasolve_;
    };



}




#endif