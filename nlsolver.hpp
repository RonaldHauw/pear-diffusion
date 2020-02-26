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
        //,lasolve_(lasolver)
        {
            std::cout<<"Non-linear solver coupled with abstract function."<<std::endl;
        }


        int solve(){
            mat_type J; J.resize(f_.size(), f_.size());
            vec_type f; f.resize(f_.size(), 1);
            for (int i = 1; i<10; i++) {
                J.setZero(); f_.J(J); f_.f(f);
                std::cout<<J.rows()<<" "<<J.cols()<<std::endl;
                std::cout<<f.rows()<<std::endl;

                vec_type x = J.colPivHouseholderQr().solve(f);
                f_.set_cons(x);
            }
            return 1;
        }



    private:
        f_type f_;
        //lasolve_type lasolve_;
    };



}




#endif