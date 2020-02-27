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


         diffusion(pear::component<d_type, vec_type> & comp, pear::grid<d_type> & grid, std::vector<d_type> param)
         : comp_(comp)
         , sigma_r_(param[0])
         , sigma_z_(param[1])
         , r_(param[2])
         , C_amb_(param[3])
         , grid_(grid)
         {
             std::cout<<"Component: "<<comp_.name()<<" will diffuse with sigma_r = "<<sigma_r_<< " and sigma_z = " <<sigma_z_<< std::endl;
             std::cout<<"  C_amb = "<<C_amb_<<"   r = "<<r_<<std::endl;
         }


         void J(mat_type & K) const {


             for (int t = 1; t<grid_.nb_elements(); t++){
                 std::vector<int> elem_nodes = grid_.element(t);

                 d_type r1   = grid_.node(elem_nodes[0])[0];     d_type z1 = grid_.node(elem_nodes[0])[1];
                 d_type r2   = grid_.node(elem_nodes[1])[0];     d_type z2 = grid_.node(elem_nodes[1])[1];
                 d_type r3   = grid_.node(elem_nodes[2])[0];     d_type z3 = grid_.node(elem_nodes[2])[1];

                 d_type omega = ((r2-r1)*(z3-z1)-(r3-r1)*(z2-z1))*0.5;

                 d_type sum_r = r1+r2+r3;

                 d_type C_12_1 = 1./omega * (z1-z3)*(z3-z2) / 12. ;
                 d_type C_12_2 = 1./omega * (r1-r3)*(r3-r2) / 12. ;

                 d_type C_23_1 = 1./omega * (z2-z1)*(z1-z3) / 12. ;
                 d_type C_23_2 = 1./omega * (r2-r1)*(r1-r3) / 12. ;


                 d_type C_13_1 = 1./omega * (z1-z2)*(z2-z3) / 12. ;
                 d_type C_13_2 = 1./omega * (r1-r2)*(r2-r3) / 12. ;

                 d_type C_11_1 = 1./omega * (z2-z3)*(z2-z3) / 12.;
                 d_type C_11_2 = 1./omega * (r2-r3)*(r2-r3) / 12.;

                 d_type C_22_1 = 1./omega * (z1-z3)*(z1-z3) / 12.;
                 d_type C_22_2 = 1./omega * (r1-r3)*(r1-r3) / 12.;

                 d_type C_33_1 = 1./omega * (z1-z2)*(z1-z2) / 12.;
                 d_type C_33_2 = 1./omega * (r1-r2)*(r1-r2) / 12.;

                 K(elem_nodes[0]-1, elem_nodes[1]-1) =  K(elem_nodes[0]-1, elem_nodes[1]-1) + sigma_r_*C_12_1+sigma_z_*C_12_2 * sum_r;
                 K(elem_nodes[1]-1, elem_nodes[0]-1) =  K(elem_nodes[1]-1, elem_nodes[0]-1) + sigma_r_*C_12_1+sigma_z_*C_12_2 * sum_r;

                 K(elem_nodes[1]-1, elem_nodes[2]-1) =  K(elem_nodes[1]-1, elem_nodes[2]-1) + sigma_r_*C_23_1+sigma_z_*C_23_2 * sum_r;
                 K(elem_nodes[2]-1, elem_nodes[1]-1) =  K(elem_nodes[2]-1, elem_nodes[1]-1) + sigma_r_*C_23_1+sigma_z_*C_23_2 * sum_r;

                 K(elem_nodes[0]-1, elem_nodes[2]-1) =  K(elem_nodes[0]-1, elem_nodes[2]-1) + sigma_r_*C_13_1+sigma_z_*C_13_2 * sum_r;
                 K(elem_nodes[2]-1, elem_nodes[0]-1) =  K(elem_nodes[2]-1, elem_nodes[0]-1) + sigma_r_*C_13_1+sigma_z_*C_13_2 * sum_r;

                 K(elem_nodes[0]-1, elem_nodes[0]-1) =  K(elem_nodes[0]-1, elem_nodes[0]-1) + sigma_r_*C_11_1+sigma_z_*C_11_2 * sum_r;

                 K(elem_nodes[1]-1, elem_nodes[1]-1) =  K(elem_nodes[1]-1, elem_nodes[1]-1) + sigma_r_*C_22_1+sigma_z_*C_22_2 * sum_r;

                 K(elem_nodes[2]-1, elem_nodes[2]-1) =  K(elem_nodes[2]-1, elem_nodes[2]-1) + sigma_r_*C_33_1+sigma_z_*C_33_2 * sum_r;

             };

         };

         void f(vec_type & f_vector, mat_type & K) const{
             // Outer boundary
             for (int t = 0; t<grid_.nb_outer_edges(); t++) {
                 std::vector<int> edge_nodes = grid_.outer_edge(t);

                 d_type r1     = grid_.node(edge_nodes[0])[0];   d_type z1 = grid_.node(edge_nodes[0])[1];
                 d_type r2     = grid_.node(edge_nodes[1])[0];   d_type z2 = grid_.node(edge_nodes[1])[1];
                 d_type length = sqrt((r1-r2)*(r1-r2) + (z1-z2)*(z1-z2));

                 f_vector( edge_nodes[0] )   = f_vector( edge_nodes[0] )   + r_ * C_amb_ * (2*r1+r2) * length / 6. ;
                 f_vector( edge_nodes[1] )   = f_vector( edge_nodes[1] )   + r_ * C_amb_ * (r1+2*r2) * length / 6. ;

             };

             // Inner boundary: in case of non-trivial Neumann conditions on the inner booundary, insert similar loop here

             // Add the matrix vector product
             K.setZero(); J(K);
             f_vector = f_vector+ K*comp_.concentrations();
             std::cout<<"concentrations"<<comp_.concentrations()<<std::endl;

         };



         void set_cons(vec_type & x){
             comp_.concentrations() = x;
         }

         int nb_nodes(){
             return comp_.nb_nodes();
         }


     private:

         pear::component<d_type, vec_type> & comp_;
         d_type sigma_r_;
         d_type sigma_z_;
         d_type r_;
         d_type C_amb_;
         pear::grid<d_type> & grid_;

     };



}




#endif