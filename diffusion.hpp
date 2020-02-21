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

         diffusion(pear::component<d_type> & comp, pear::grid<d_type> & grid, d_type sigma_r, d_type sigma_z)
         : comp_(comp)
         , sigma_r_(sigma_r)
         , sigma_z_(sigma_z)
         , grid_(grid)
         {
             std::cout<<"Component: "<<comp_.name()<<" will diffuse with sigma_r = "<<sigma_r_<< " and sigma_z = " <<sigma_z_<< std::endl;
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
         pear::component<d_type> & comp_;
         d_type sigma_r_;
         d_type sigma_z_;
         pear::grid<d_type> & grid_;

     };



}




#endif