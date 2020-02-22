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


         diffusion(pear::component<d_type, vec_type> & comp, pear::grid<d_type> & grid, d_type sigma_r, d_type sigma_z)
         : comp_(comp)
         , sigma_r_(sigma_r)
         , sigma_z_(sigma_z)
         , grid_(grid)
         {
             std::cout<<"Component: "<<comp_.name()<<" will diffuse with sigma_r = "<<sigma_r_<< " and sigma_z = " <<sigma_z_<< std::endl;
         }


         void J(){


             for (int t = 0; t<grid_.elements_length(); t++){
                //elem_nodes = grid_.element(t); // vector of length 3

                //elem_node_1_coords = grid_.node(elem_nodes(1)); // vector of length 2
                //elem_node_2_coords = grid_.node(elem_nodes(2));
                //elem_node_3_coords = grid_.node(elem_nodes(3));

                //d_type r1 = elem_node_1_coords(1);
                //d_type r2 = elem_node_2_coords(1);
                //d_type r3 = elem_node_3_coords(1);

                //d_type z1 = elem_node_1_coords(2);
                //d_type z2 = elem_node_2_coords(2);
                //d_type z3 = elem_node_3_coords(2);

                //d_type omega = ((r2-r1)*(z3-z1)-(r3-r1)*(z2-z1))*0.5;

                //d_type sum_r = r1+r2+r3;

                //d_type C_12_1 = 1.0/6 * 1.0/2/omega * (z1-z3)*(z3-z2) ;
                //d_type C_12_2 = 1.0/6 * 1.0/2/omega * (r1-r3)*(r3-r2) ;

                //d_type C_23_1 = 1.0/6 * 1.0/2/omega * (z2-z1)*(z1-z3) ;
                //d_type C_23_2 = 1.0/6 * 1.0/2/omega * (r2-r1)*(r1-r3) ;

                //d_type C_13_1 = 1.0/6 * 1.0/2/omega * (z1-z2)*(z2-z3) ;
                //d_type C_13_2 = 1.0/6 * 1.0/2/omega * (r1-r2)*(r2-r3) ;



             };

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
         d_type sigma_r_;
         d_type sigma_z_;
         pear::grid<d_type> & grid_;

     };



}




#endif