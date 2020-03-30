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
#include "component.hpp"
#include "eigen/Eigen/Dense"


namespace pear {

    template <typename d_type, typename vec_type, typename mat_type>
     class diffusion{
     public:

        /* Constructor for diffusion
         *
         * IN :    -A #component# describing the component to be diffused
         *         -A #grid# describing the mesh of the problem
         *         -A #vector# (of length 4X1) describing the different physical parameters necessary in the diffusion
         *                 model; sigma_r, sigma_z, r, c_amb
         * OUT :   the private variables  which correspond to the above mentioned quantities
         */
         diffusion(pear::component<d_type, vec_type, mat_type> & comp, pear::grid<d_type, mat_type> & grid, std::vector<d_type> param)
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

         /* Function evaluation of the diffusion model
         *
         * IN: A #reference(vector)# describing the value of the diffusion equation for this particular component
         * OUT: void. Upon completion, the contribution from the diffusion equation has been added to the #vector#
         *                 referenced
         *
         */
         void f(Eigen::Ref<vec_type> f_vector) const{

             // Outer boundary
             for (int t = 1; t<grid_.nb_outer_edges()+1; t++) {
                 std::vector<int> edge_nodes = grid_.outer_edge(t);

                 d_type r1     = grid_.node(edge_nodes[0])[0];   d_type z1 = grid_.node(edge_nodes[0])[1];
                 d_type r2     = grid_.node(edge_nodes[1])[0];   d_type z2 = grid_.node(edge_nodes[1])[1];
                 d_type length = sqrt((r1-r2)*(r1-r2) + (z1-z2)*(z1-z2));

                 f_vector( edge_nodes[0]-1 ) += r_ * C_amb_ * (2*r1+r2) * length / 6. ;
                 f_vector( edge_nodes[1]-1 ) += r_ * C_amb_ * (r1+2*r2) * length / 6. ;

             };
             // Inner boundary: in case of non-trivial Neumann conditions on the inner booundary, insert similar loop here
         };

         /* Jacobian of the diffusion model
         *
         * IN: A #sparse_matrix# describing the Jacobian of the (non-linear) equation
         * OUT: void. Upon completion, the contribution from the diffusion to the Jacobian of the (non-linear) equation
         *                 from this particular component (K_i) has been added to the #sparse_matrix#
         *
         */
         void J(mat_type & K) const {

             int i = comp_.start();

             for (int t = 1; t<grid_.nb_elements()+1; t++){

                 std::vector<int> elem_nodes = grid_.element(t);

                 d_type r1   = grid_.node(elem_nodes[0])[0];     d_type z1 = grid_.node(elem_nodes[0])[1];
                 d_type r2   = grid_.node(elem_nodes[1])[0];     d_type z2 = grid_.node(elem_nodes[1])[1];
                 d_type r3   = grid_.node(elem_nodes[2])[0];     d_type z3 = grid_.node(elem_nodes[2])[1];

                 d_type omega = ((r2-r1)*(z3-z1)-(r3-r1)*(z2-z1)) * 0.5;
                 d_type sum_r = r1+r2+r3;

                 d_type C_12_1 = (z1-z3)*(z3-z2) / 12. /omega ;
                 d_type C_12_2 = (r1-r3)*(r3-r2) / 12. /omega ;

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


                 K.coeffRef(elem_nodes[0]-1 + i, elem_nodes[1]-1 + i) += (sigma_r_*C_12_1+sigma_z_*C_12_2) * sum_r;
                 K.coeffRef(elem_nodes[1]-1 + i, elem_nodes[0]-1 + i) += (sigma_r_*C_12_1+sigma_z_*C_12_2) * sum_r;

                 K.coeffRef(elem_nodes[1]-1 + i, elem_nodes[2]-1 + i) += (sigma_r_*C_23_1+sigma_z_*C_23_2) * sum_r;
                 K.coeffRef(elem_nodes[2]-1 + i, elem_nodes[1]-1 + i) += (sigma_r_*C_23_1+sigma_z_*C_23_2) * sum_r;

                 K.coeffRef(elem_nodes[0]-1 + i, elem_nodes[2]-1 + i) += (sigma_r_*C_13_1+sigma_z_*C_13_2) * sum_r;
                 K.coeffRef(elem_nodes[2]-1 + i, elem_nodes[0]-1 + i) += (sigma_r_*C_13_1+sigma_z_*C_13_2) * sum_r;

                 K.coeffRef(elem_nodes[0]-1 + i, elem_nodes[0]-1 + i) += (sigma_r_*C_11_1+sigma_z_*C_11_2) * sum_r;
                 K.coeffRef(elem_nodes[1]-1 + i, elem_nodes[1]-1 + i) += (sigma_r_*C_22_1+sigma_z_*C_22_2) * sum_r;

                 K.coeffRef(elem_nodes[2]-1 + i, elem_nodes[2]-1 + i) += (sigma_r_*C_33_1+sigma_z_*C_33_2) * sum_r;

             };
             for (int t = 1; t<grid_.nb_outer_edges()+1; t++) {
                 std::vector<int> edge_nodes = grid_.outer_edge(t);

                 d_type r1   = grid_.node(edge_nodes[0])[0];     d_type z1 = grid_.node(edge_nodes[0])[1];
                 d_type r2   = grid_.node(edge_nodes[1])[0];     d_type z2 = grid_.node(edge_nodes[1])[1];

                 d_type len = sqrt( (r2-r1)*(r2-r1) + (z2-z1)*(z2-z1) );

                 d_type parallel_term_1     = 1./12. * len * ( 3*r1 +   r2 ) ;
                 d_type parallel_term_2     = 1./12. * len * (   r1 + 3*r2 ) ;
                 d_type cross_term          = 1./12. * len * (   r1 +   r2 ) ;


                 K.coeffRef( edge_nodes[0]-1+i, edge_nodes[0]-1+i )  += r_ * parallel_term_1 ;
                 K.coeffRef( edge_nodes[0]-1+i, edge_nodes[1]-1+i )  += r_ * cross_term      ;
                 K.coeffRef( edge_nodes[1]-1+i, edge_nodes[0]-1+i )  += r_ * cross_term      ;
                 K.coeffRef( edge_nodes[1]-1+i, edge_nodes[1]-1+i )  += r_ * parallel_term_2 ;

             }

         };

         void set_cons(vec_type & x){           comp_.concentrations() = x;}
         Eigen::Ref<vec_type> cons(){           return comp_.cons();}
         Eigen::Ref<vec_type> cons_full(){      return comp_.cons_full();}
         int cons_start(){                      return comp_.cons_start(); };
         int cons_stop(){                       return comp_.cons_stop(); };
         int cons_stride(){                     return comp_.cons_stride(); };

         int nb_nodes(){
             return comp_.nb_nodes();
         }

         void setSparsityPattern(mat_type & m){
             grid_.setSparsityPattern(m);
         }

     private:

         pear::component<d_type, vec_type, mat_type> & comp_;
         d_type sigma_r_;
         d_type sigma_z_;
         d_type r_;
         d_type C_amb_;
         pear::grid<d_type, mat_type> & grid_;

     };



}




#endif