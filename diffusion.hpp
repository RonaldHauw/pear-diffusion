#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP

#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "grid.hpp"

namespace pear {

    template <typename d_type>
     class diffusion{
     public:

         diffusion(pear::component comp, pear::grid<d_type> & grid, d_type sigma)
         : comp_(comp)
         , sigma_(sigma)
         , grid_(grid)
         {
             std::cout<<"Component: "<<comp_.name()<<" will diffuse with sigma = "<<sigma_<<std::endl;
         }

         /* evaluate
          *
          * returns the diffusion at each grid point given the concentrations c in a vector d.
          *
          * can be interpreted as the calculation of the matrix vector product K_c*c = d, where K_c the diffusion
          * matrix for component c.
          */
         void evaluate(std::vector<d_type> c, std::vector<d_type> d){
             std::cout<<"diffusion.evaluate is empty"<<std::endl;
         }

         /* d_matrix
          *
          * returns the diffusion matrix. Useful for debugging purposes.
          */
         void d_matrix(std::vector<d_type> c, std::vector<d_type> d){
             std::cout<<"diffusion.evaluate is empty"<<std::endl;
         }

     private:
         pear::component comp_;
         d_type sigma_;
         pear::grid<d_type> & grid_;

     };



}




#endif