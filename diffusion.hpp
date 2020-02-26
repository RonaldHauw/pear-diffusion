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

    template <typename d_type, typename vec_type, typename mat_type>
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


         void J(mat_type K) const {


             for (int t = 1; t<grid_.nb_elements(); t++){
                 std::vector<d_type> elem_nodes = grid_.element(t); // vector of length 3 [node_1_n, node_2_n, node_3_n]

                 d_type r1   = grid_.node(elem_nodes[0])[0];     d_type z1 = grid_.node(elem_nodes[0])[1];
                 d_type r2   = grid_.node(elem_nodes[1])[0];     d_type z2 = grid_.node(elem_nodes[1])[1];
                 d_type r3   = grid_.node(elem_nodes[2])[0];     d_type z3 = grid_.node(elem_nodes[2])[1];

                 d_type omega = ((r2-r1)*(z3-z1)-(r3-r1)*(z2-z1))*0.5;

                 d_type sum_r = r1+r2+r3;

                 d_type C_12_1 = 1.0/6 * 1.0/2/omega * (z1-z3)*(z3-z2) ;
                 d_type C_12_2 = 1.0/6 * 1.0/2/omega * (r1-r3)*(r3-r2) ;

                 d_type C_23_1 = 1.0/6 * 1.0/2/omega * (z2-z1)*(z1-z3) ;
                 d_type C_23_2 = 1.0/6 * 1.0/2/omega * (r2-r1)*(r1-r3) ;

                 d_type C_13_1 = 1.0/6 * 1.0/2/omega * (z1-z2)*(z2-z3) ;
                 d_type C_13_2 = 1.0/6 * 1.0/2/omega * (r1-r2)*(r2-r3) ;

                 d_type C_11_1 = 1.0/6 * 1.0/2/omega * (z2-z3)*(z2-z3);
                 d_type C_11_2 = 1.0/6 * 1.0/2/omega * (r2-r3)*(r2-r3);

                 d_type C_22_1 = 1.0/6 * 1.0/2/omega * (z1-z3)*(z1-z3);
                 d_type C_22_2 = 1.0/6 * 1.0/2/omega * (r1-r3)*(r1-r3);

                 d_type C_33_1 = 1.0/6 * 1.0/2/omega * (z1-z2)*(z1-z2);
                 d_type C_33_2 = 1.0/6 * 1.0/2/omega * (r1-r2)*(r1-r2);

                 K(elem_nodes[1], elem_nodes[2]) =  K(elem_nodes[1], elem_nodes[2]) + sigma_r_*C_12_1+sigma_z_*C_12_2 * sum_r;
                 K(elem_nodes[2], elem_nodes[1]) =  K(elem_nodes[2], elem_nodes[1]) + sigma_r_*C_12_1+sigma_z_*C_12_2 * sum_r;

                 K(elem_nodes[2], elem_nodes[3]) =  K(elem_nodes[2], elem_nodes[3]) + sigma_r_*C_23_1+sigma_z_*C_23_2 * sum_r;
                 K(elem_nodes[3], elem_nodes[2]) =  K(elem_nodes[3], elem_nodes[2]) + sigma_r_*C_23_1+sigma_z_*C_23_2 * sum_r;

                 K(elem_nodes[1], elem_nodes[3]) =  K(elem_nodes[1], elem_nodes[3]) + sigma_r_*C_13_1+sigma_z_*C_13_2 * sum_r;
                 K(elem_nodes[3], elem_nodes[1]) =  K(elem_nodes[3], elem_nodes[1]) + sigma_r_*C_13_1+sigma_z_*C_13_2 * sum_r;

                 K(elem_nodes[1], elem_nodes[1]) =  K(elem_nodes[1], elem_nodes[1]) + sigma_r_*C_11_1+sigma_z_*C_11_2 * sum_r;

                 K(elem_nodes[2], elem_nodes[2]) =  K(elem_nodes[2], elem_nodes[2]) + sigma_r_*C_22_1+sigma_z_*C_22_2 * sum_r;

                 K(elem_nodes[3], elem_nodes[3]) =  K(elem_nodes[3], elem_nodes[3]) + sigma_r_*C_33_1+sigma_z_*C_33_2 * sum_r;

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