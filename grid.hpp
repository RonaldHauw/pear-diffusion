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
        grid(std::string file_name)
        : file_name_(file_name)
        {
            // NODES
            std::ifstream myfile(file_name_ + "_Nodes.txt"); // read file
            d_type x, y; int n, N;
            if (myfile.is_open()){
                // count number of grid points
                N = 0; while(myfile>>n>>x>>y){N++;}; nb_nodes_ = N;
                nodes_ = std::vector<d_type>(2*nb_nodes_);
                // reopen
                myfile.close();
                std::ifstream myfile(file_name_+"_Nodes.txt");
                // save coordinates
                while(myfile>>n>>x>>y){
                    nodes_[2*n] = x;
                    nodes_[2*n+1] = y;
                }

                std::cout<<"Grid nodes loaded from: "<<file_name_<<" length: "<<nb_nodes_<<std::endl;
                // ELEMENTS
            }

            else {std::cout << "Unable to open file, check the executable folder";}

            // ELEMENTS
            std::ifstream myfile2(file_name_ + "_Elements.txt"); // read file
            if (myfile2.is_open()){
                int n1, n2, n3;
                // count number of grid points
                N = 0; while(myfile2>>n>>n1>>n2>>n3){N++;}; nb_elements_ = N;
                elements_ = std::vector<d_type>(3*nb_elements_);
                // reopen
                myfile2.close();
                std::ifstream myfile2(file_name_+"_Nodes.txt");
                // save coordinates
                while(myfile2>>n>>n1>>n2>>n3){
                    elements_[3*n] = n1;
                    elements_[3*n+1] = n2;
                    elements_[3*n+2] = n3;
                }
                std::cout<<"Grid elements loaded from: "<<file_name_<<" length: "<<nb_elements_<<std::endl;
            }

            else {std::cout << "Unable to open file, check the executable folder";}

            // ELEMENTS
            std::ifstream myfile2(file_name_ + "_Elements.txt"); // read file
            if (myfile2.is_open()){
                int n1, n2, n3;
                // count number of grid points
                N = 0; while(myfile2>>n>>n1>>n2>>n3){N++;}; nb_elements_ = N;
                elements_ = std::vector<d_type>(3*nb_elements_);
                // reopen
                myfile2.close();
                std::ifstream myfile2(file_name_+"_Nodes.txt");
                // save coordinates
                while(myfile2>>n>>n1>>n2>>n3){
                    elements_[3*n] = n1;
                    elements_[3*n+1] = n2;
                    elements_[3*n+2] = n3;
                }
                std::cout<<"Grid elements loaded from: "<<file_name_<<" length: "<<nb_elements_<<std::endl;
            }

            else {std::cout << "Unable to open file, check the executable folder";}
            // reorder grid points

        }

        void get_coords(d_type & x, d_type & y, int i) const {
            x = nodes_(2*i);
            y = nodes_(2*i+1);
        }

        int length() const {
            return nb_nodes_;
        }



    private:
        std::string file_name_;
        std::vector<d_type> nodes_; //  length 2N_: x1 y1 x2 y2 x3 y3 for sequantial mem access
        std::vector<d_type> edges_;
        std::vector<d_type> elements_;
        std::vector<d_type> boundaries_;
        int nb_nodes_; // nb of grid points
        int nb_elements_; // nb of elements

    };



}




#endif