#ifndef REACTION_HPP
#define REACTION_HPP

#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "grid.hpp"
#include "eigen/Eigen/Dense"
#include "component.hpp"
#include "diffusion.hpp"
#include "respiration.hpp"

namespace pear {


    template <typename d_type, typename vec_type, typename mat_type>
    class reaction{
    public:

        /* Constructor for reaction
         *
         * IN :    -Two #component#'s describing the components to react
         *         -A #grid# describing the mesh of the problem
         *         -A #respiration# which contains the reaction dynamics
         * OUT :   the private variables  which correspond to the above mentioned quantities
         */
        reaction(
                pear::component<d_type, vec_type, mat_type> & co2,
                pear::component<d_type, vec_type, mat_type> & o2,
                pear::grid<d_type, mat_type> & grid,
                pear::respiration<d_type> resp)
                : o2_(o2)
                , co2_(co2)
                , grid_(grid)
                , resp_(resp)
        {
            std::cout<<"// RESPIRATION //"<<std::endl<<"        Between O2 and CO2"<<std::endl;
        }

        // Respiration dynamics: functions
        d_type R_u(d_type C_u, d_type C_v){
            return resp_.R_u(C_u, C_v); };
        d_type R_v(d_type C_u, d_type C_v){
            return resp_.R_v(C_u, C_v);};

        // Respiration dynamics: derivatives
        d_type dR_u_u(d_type C_u, d_type C_v){
            return resp_.dR_u_u(C_u, C_v);};
        d_type dR_u_v(d_type C_u, d_type C_v){
            return resp_.dR_u_v(C_u, C_v); };
        d_type dR_v_u(d_type C_u, d_type C_v){
            return resp_.dR_v_u(C_u, C_v); };
        d_type dR_v_v(d_type C_u, d_type C_v){
            return resp_.dR_v_v(C_u, C_v); };

        /* Function evaluation of the reaction model
        *
        * IN: Two #reference(vector)# describing the value of the reaction equation for each of the particular component
        * OUT: void. Upon completion, the contribution from the reaction equation has been added to the #vector#'s
        *                 referenced
        *
        */
        void f(Eigen::Ref<vec_type>  H_o2, Eigen::Ref<vec_type>  H_co2) {
            for (int t = 1; t < grid_.nb_elements()+1; t++) {

                std::vector<int> elem_nodes = grid_.element(t);

                d_type r1 = grid_.node(elem_nodes[0])[0];   d_type z1 = grid_.node(elem_nodes[0])[1];
                d_type r2 = grid_.node(elem_nodes[1])[0];   d_type z2 = grid_.node(elem_nodes[1])[1];
                d_type r3 = grid_.node(elem_nodes[2])[0];   d_type z3 = grid_.node(elem_nodes[2])[1];

                d_type rm = (r1+r2+r3) / 3.;
                d_type co2m = (o2_.concentration(elem_nodes[0])+o2_.concentration(elem_nodes[1])+o2_.concentration(elem_nodes[2]))/3.;
                d_type cco2m = (co2_.concentration(elem_nodes[0])+co2_.concentration(elem_nodes[1])+co2_.concentration(elem_nodes[2]))/3.;


                d_type area = (r2*z3 - r3*z2) - (r1*z3-r3*z1) + (r1*z2-z1*r2) ;

                /* uncomment for quadrature formula on edges
                H_o2(elem_nodes[0]-1) += r1 * R_u(o2_.concentration(elem_nodes[0]), co2_.concentration(elem_nodes[0])) * area / 6. ;
                H_o2(elem_nodes[1]-1) += r2 * R_u(o2_.concentration(elem_nodes[1]), co2_.concentration(elem_nodes[1])) * area / 6. ;
                H_o2(elem_nodes[2]-1) += r3 * R_u(o2_.concentration(elem_nodes[2]), co2_.concentration(elem_nodes[2])) * area / 6. ;


                H_co2(elem_nodes[0]-1) -= r1 * R_v(o2_.concentration(elem_nodes[0]), co2_.concentration(elem_nodes[0])) * area / 6. ;
                H_co2(elem_nodes[1]-1) -= r2 * R_v(o2_.concentration(elem_nodes[1]), co2_.concentration(elem_nodes[1])) * area / 6. ;
                H_co2(elem_nodes[2]-1) -= r3 * R_v(o2_.concentration(elem_nodes[2]), co2_.concentration(elem_nodes[2])) * area / 6. ;
                */

                // quadrature formula in center
                H_o2(elem_nodes[0]-1) += rm * R_u(co2m, cco2m) * area / 6. ;
                H_o2(elem_nodes[1]-1) += rm * R_u(co2m, cco2m) * area / 6. ;
                H_o2(elem_nodes[2]-1) += rm * R_u(co2m, cco2m) * area / 6. ;


                H_co2(elem_nodes[0]-1) -= rm * R_v(co2m, cco2m) * area / 6. ;
                H_co2(elem_nodes[1]-1) -= rm * R_v(co2m, cco2m) * area / 6. ;
                H_co2(elem_nodes[2]-1) -= rm * R_v(co2m, cco2m) * area / 6. ;

            };
        };

        /* Jacobian of the diffusion model
        *
        * IN: A #sparse_matrix# describing the Jacobian of the (non-linear) equation
        * OUT: void. Upon completion, the contribution from the reaction to the Jacobian of the (non-linear) equation
        *                 from this particular component (H_i) has been added to the #sparse_matrix#
        *
        */
        void J(mat_type & dH) {

            /*  uncomment for quadrature formula in edges

            for (int t = 1; t < grid_.nb_nodes()+1; t++) {
                std::vector<int> elements = grid_.elements_for_node(t);
                d_type area = 0.;
                for (int e = 1; e < static_cast<signed int>(elements.size())+1; e++) {

                        std::vector<int> nodes = grid_.element(elements[e - 1]);

                        d_type r1 = grid_.node(nodes[0])[0]; d_type z1 = grid_.node(nodes[0])[1];
                        d_type r2 = grid_.node(nodes[1])[0]; d_type z2 = grid_.node(nodes[1])[1];
                        d_type r3 = grid_.node(nodes[2])[0]; d_type z3 = grid_.node(nodes[2])[1];
                        area += (r2 * z3 - r3 * z2) - (r1 * z3 - r3 * z1) + (r1 * z2 - z1 * r2);

                };

                dH.coeffRef(t-1, t-1)                                            +=  area * grid_.node(t)[0] * dR_u_u(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH.coeffRef(t-1, t + grid_.nb_nodes() -1)                        +=  area * grid_.node(t)[0] * dR_u_v(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH.coeffRef(t + grid_.nb_nodes() -1, t-1)                        += -area * grid_.node(t)[0] * dR_v_u(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
                dH.coeffRef(t + grid_.nb_nodes() -1, t + grid_.nb_nodes() -1)    += -area * grid_.node(t)[0] * dR_v_v(o2_.concentration(t), co2_.concentration(t) ) / 6. ;
            };
             */

            // quadrature formula in center
            for (int t = 1; t<grid_.nb_elements()+1; t++) {

                std::vector<int> elem_nodes = grid_.element(t);

                d_type r1 = grid_.node(elem_nodes[0])[0];
                d_type z1 = grid_.node(elem_nodes[0])[1];
                d_type r2 = grid_.node(elem_nodes[1])[0];
                d_type z2 = grid_.node(elem_nodes[1])[1];
                d_type r3 = grid_.node(elem_nodes[2])[0];
                d_type z3 = grid_.node(elem_nodes[2])[1];

                d_type area = ((r2*z3 - r3*z2) - (r1*z3-r3*z1) + (r1*z2-z1*r2))/2. ;

                d_type rm = (r1+r2+r3) / 3.;
                d_type co2m = (o2_.concentration(elem_nodes[0])+o2_.concentration(elem_nodes[1])+o2_.concentration(elem_nodes[2]))/3.;
                d_type cco2m = (co2_.concentration(elem_nodes[0])+co2_.concentration(elem_nodes[1])+co2_.concentration(elem_nodes[2]))/3.;



                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[0]-1+ o2_.cons_start())   += area/9. * rm * dR_u_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[1]-1+ o2_.cons_start())   += area/9. * rm * dR_u_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[2]-1+ o2_.cons_start())   += area/9. * rm * dR_u_u( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[1]-1+ o2_.cons_start(), elem_nodes[0]-1+ o2_.cons_start())   += area/9. * rm * dR_u_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[1]-1+ o2_.cons_start(), elem_nodes[1]-1+ o2_.cons_start())   += area/9. * rm * dR_u_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[1]-1+ o2_.cons_start(), elem_nodes[2]-1+ o2_.cons_start())   += area/9. * rm * dR_u_u( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[2]-1+ o2_.cons_start(), elem_nodes[0]-1+ o2_.cons_start())     += area/9. * rm * dR_u_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[2]-1+ o2_.cons_start(), elem_nodes[1]-1+ o2_.cons_start())     += area/9. * rm * dR_u_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[2]-1+ o2_.cons_start(), elem_nodes[2]-1+ o2_.cons_start())     += area/9. * rm * dR_u_u( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[0]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[1]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[2]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[0]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[1]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[2]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[0]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[1]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ o2_.cons_start(), elem_nodes[2]-1 + co2_.cons_start())   += area/9. * rm * dR_u_v( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[0]-1+ co2_.cons_start(), elem_nodes[0]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ co2_.cons_start(), elem_nodes[1]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ co2_.cons_start(), elem_nodes[2]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[1]-1+ co2_.cons_start(), elem_nodes[0]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[1]-1+ co2_.cons_start(), elem_nodes[1]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[1]-1+ co2_.cons_start(), elem_nodes[2]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[2]-1+ co2_.cons_start(), elem_nodes[0]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[2]-1+ co2_.cons_start(), elem_nodes[1]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[2]-1+ co2_.cons_start(), elem_nodes[2]-1+ o2_.cons_start())    += -area/9. * rm * dR_v_u( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[0]-1+ co2_.cons_start(), elem_nodes[0]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ co2_.cons_start(), elem_nodes[1]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[0]-1+ co2_.cons_start(), elem_nodes[2]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[1]-1+ co2_.cons_start(), elem_nodes[0]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[1]-1+ co2_.cons_start(), elem_nodes[1]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[1]-1+ co2_.cons_start(), elem_nodes[2]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;

                dH.coeffRef(elem_nodes[2]-1+ co2_.cons_start(), elem_nodes[0]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[2]-1+ co2_.cons_start(), elem_nodes[1]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;
                dH.coeffRef(elem_nodes[2]-1+ co2_.cons_start(), elem_nodes[2]-1 + co2_.cons_start())  += -area/9. * rm * dR_v_v( co2m, cco2m ) ;


            }

        };

        void suppress_nonlinearity(d_type alpha){
            resp_.suppress_nonlinearity(alpha);
        };

    private:
        pear::component<d_type, vec_type, mat_type> & o2_;
        pear::component<d_type, vec_type, mat_type> & co2_;
        pear::grid<d_type, mat_type> & grid_;
        pear::respiration<d_type> resp_;
    };




}



#endif