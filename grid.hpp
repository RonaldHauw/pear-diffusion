#ifndef GRID_HPP
#define GRID_HPP

#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>


namespace pear {

    template <typename d_type>
    class grid{
    public:


        /* constructor for grid
         *
         * loads coordinates from a text file in the following format:
         * int grid_number d_type x_coordinate d_type y_coordinate
         *
         * for example:
         * 1 0.0 0.0
         * 2 0.5 0.5
         */
        grid(const char * file_name)
        : file_name_(file_name)
        {
            // read file
            std::ifstream myfile(file_name_); // read file

            d_type x, y; int n, N;

            if (myfile.is_open()){
                // count number of grid points
                N = 0;
                while(myfile>>n>>x>>y){
                    N++;
                }
                N_ = N;

                grid_coordinates_ = std::vector<d_type>(2*N_);

                // reopen
                myfile.close();
                std::ifstream myfile(file_name_);
                // save coordinates
                while(myfile>>n>>x>>y){
                    grid_coordinates_[2*n] = x;
                    grid_coordinates_[2*n+1] = y;
                }

                std::cout<<"Grid loaded from: "<<file_name_<<" length: "<<N_<<std::endl;
            }

            else {std::cout << "Unable to open file, check the executable folder";}
            // reorder grid points

        }

        void get_coords(d_type & x, d_type & y, int i) const {
            x = grid_coordinates_(2*i);
            y = grid_coordinates_(2*i+1);
        }

        int length() const {
            return N_;
        }



    private:
        std::string file_name_;
        std::vector<d_type> grid_coordinates_; //  length 2N_: x1 y1 x2 y2 x3 y3 for sequantial mem access
        int N_; // nb of grid points

    };



}




#endif