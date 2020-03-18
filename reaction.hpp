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


namespace pear {


    template <typename d_type, typename vec_type>
    class respiration{
    public:

        respiration(pear::component<d_type, vec_type> & co2,
                pear::component<d_type, vec_type> & o2,
                pear::grid<d_type> & grid,
                d_type p1,
                d_type p2,
                d_type p3,
                d_type p4,
                d_type p5,
                d_type p6
        )
                : o2_(o2)
                , co2_(co2)
                , grid_(grid)
                , p1_(p1)
                , p2_(p2)
                , p3_(p3)
                , p4_(p4)
                , p5_(p5)
                , p6_(p6)
        {
            std::cout<<"Respiration between O2 and CO2"<<std::endl;
        }

        void eval_co2(){
            std::cout<<"empty function"<<std::endl;
        }

        void eval_o2(){
            std::cout<<"empty function"<<std::endl;
        }

        void J_co2(){
            std::cout<<"empty function"<<std::endl;
        }

        void J_o2(){
            std::cout<<"empty function"<<std::endl;
        }


    private:
        pear::component<d_type, vec_type> & o2_;
        pear::component<d_type, vec_type> & co2_;
        pear::grid<d_type> & grid_;
        d_type p1_, p2_, p3_, p4_, p5_, p6_;

    };




}



#endif