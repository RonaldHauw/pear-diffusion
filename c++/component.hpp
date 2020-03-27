#ifndef COMPONENTS_HPP
#define COMPONENTS_HPP
// hogere orde benadering driehoekjes
// reeksontwikkeling niet lineariteit in basis functies (termen vallen weg)
// discretisatie step te groot (convergentie plot)
// minteken verkeerd ergens
// boundary integraal juist?
// expressie templates eigen aan/uit?

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
#include "eigen/Eigen/Eigen"
#include "eigen/Eigen/Core"


namespace pear {

    template <typename d_type, typename vec_type>
    class component{
    public:

        component(std::string name, pear::grid<d_type> & grid, vec_type & c, int start, int stop, int stride)
                : name_(name)
                , concentration_(c)
                , grid_(grid)
                , start_(start)
                , stop_(stop)
                , stride_(stride)
        {
            std::cout<<"Initialised component: "<<name_<<std::endl;
        }

        std::string name(){
            return name_;
        }



        Eigen::Ref<vec_type> cons()  {
            return concentration_.segment(start_, grid_.nb_nodes());
        }

        Eigen::Ref<vec_type> cons_full()  {
            return concentration_;
        }

        int cons_start(){ return start_; };
        int cons_stop(){ return stop_; };
        int cons_stride(){ return stride_; };

        d_type & concentration(int i){
            return concentration_(start_ + (i-1)*stride_);
        }

        int nb_nodes(){
            return grid_.nb_nodes();
        };

        void set_initial(d_type c_amb){
            this->cons().setOnes();
            this->cons() *= c_amb;
        }

    private:
        std::string name_;
        vec_type & concentration_;
        pear::grid<d_type> & grid_;
        int start_;
        int stop_;
        int stride_;

    };



}




#endif