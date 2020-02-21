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

    template <typename d_type, typename vec_type>
     class diffusion{
     public:

         diffusion(pear::component<d_type, vec_type> & comp, pear::grid<d_type> & grid, d_type sigma)
         : comp_(comp)
         , sigma_(sigma)
         , grid_(grid)
         {
             std::cout<<"Component: "<<comp_.name()<<" will diffuse with sigma = "<<sigma_<<std::endl;
         }


         void J(){

         };

         void f(){};

         void r(){};


         vec_type & get_cons(){
             return comp_.get_cons();
         }

         int length(){
             return grid_.length();
         };



     private:
         pear::component<d_type, vec_type> & comp_;
         d_type sigma_;
         pear::grid<d_type> & grid_;

     };



}




#endif