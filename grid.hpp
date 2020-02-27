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
         * Loads coordinates from a text file in the following format:

         * (int) grid_number //space// (d_type) x_coordinate //space// (d_type) y_coordinate
         *
         * for example:
         * 1 0.0 0.0
         * 2 0.5 0.5
         */
        grid(std::string file_name)
        : file_name_(file_name)
        {
            // NODES

            std::ifstream nodes_file(file_name_ + "_Nodes.txt");
            d_type x, y; int n, N, n1, n2, n3;
            if (nodes_file.is_open()){
                // Count number of grid points to pre-allocate memory
                N = 0; while(nodes_file>>n>>x>>y){N++;}; nb_nodes_ = N;
                nodes_ = std::vector<d_type>(2*nb_nodes_);

                // Reopen
                nodes_file.close();
                std::ifstream nodes_file(file_name_+"_Nodes.txt");

                // Save coordinates
                while(nodes_file>>n>>x>>y){
                    nodes_[2*(n-1)] = x;
                    nodes_[2*(n-1)+1] = y;

                }
                std::cout<<"Grid nodes loaded from: "<<file_name_<<" length: "<<nb_nodes_<<std::endl;

            }
            else {std::cout << "Unable to open file, check the executable folder";}

            // ELEMENTS

            std::ifstream elements_file(file_name_ + "_Elements.txt");
            if (elements_file.is_open()){
                // Count number of elements to pre-allocate memory
                N = 0; while(elements_file>>n>>n1>>n2>>n3){N++;}; nb_elements_ = N;
                elements_ = std::vector<d_type>(3*nb_elements_);
                // Reopen
                elements_file.close();
                std::ifstream elements_file(file_name_+"_Elements.txt");

                // Save coordinates
                while(elements_file>>n>>n1>>n2>>n3){
                    elements_[3*(n-1)] = n1;
                    elements_[3*(n-1)+1] = n2;
                    elements_[3*(n-1)+2] = n3;
                }
                std::cout<<"Grid elements loaded from: "<<file_name_<<" length: "<<nb_elements_<<std::endl;
            }
            else {std::cout << "Unable to open file, check the executable folder";}

            // OUTER EDGES

            std::ifstream outer_file(file_name_ + "_OuterEdges.txt");
            if (outer_file.is_open()){
                // Count number of outer edges to pre-allocate memory
                N = 0; while(outer_file>>n>>n1>>n2>>n3){N++;}; nb_outer_edges_ = N;
                outer_edges_ = std::vector<d_type>(2*nb_outer_edges_);
                // Reopen
                outer_file.close();
                std::ifstream outer_file(file_name_+"_OuterEdges.txt");
                // Save coordinates
                while(outer_file>>n>>n1>>n2){
                    outer_edges_[2*(n-1)] = n1;
                    outer_edges_[2*(n-1)+1] = n2;
                }
                std::cout<<"Outer edges loaded from: "<<file_name_<<" length: "<<nb_outer_edges_<<std::endl;
            }

            else {std::cout << "Unable to open file, check the executable folder";}

            // INNER EDGES

            std::ifstream inner_file(file_name_ + "_InnerEdges.txt"); // read file
            if (inner_file.is_open()){
                // Count number of inner edges to pre-allocate memory
                N = 0; while(inner_file>>n>>n1>>n2){N++;}; nb_inner_edges_ = N;
                inner_edges_ = std::vector<d_type>(2*nb_inner_edges_);
                // Reopen
                inner_file.close();
                std::ifstream inner_file(file_name_+"_OuterEdges.txt");
                // Save coordinates
                while(inner_file>>n>>n1>>n2){
                    inner_edges_[2*(n-1)] = n1;
                    inner_edges_[2*(n-1)+1] = n2;
                }
                std::cout<<"Inner edges loaded from: "<<file_name_<<" length: "<<nb_inner_edges_<<std::endl;
            }
            else {std::cout << "Unable to open file, check the executable folder";};
        }

        std::vector<d_type> node(int i) const {
            std::vector<d_type> node_coord(2);
            node_coord[0] = nodes_[2*(i-1)];
            node_coord[1] = nodes_[2*(i-1)+1];
            return node_coord;
        }

        std::vector<d_type> element(int i) const {
            std::vector<d_type> elem_nodes(3);
            elem_nodes[0] = elements_[3*(i-1)];
            elem_nodes[1] = elements_[3*(i-1)+1];
            elem_nodes[2] = elements_[3*(i-1)+2];
            return elem_nodes;
        }

        std::vector<d_type> outer_edge(int i) const {
            std::vector<d_type> edge_nodes(2);
            edge_nodes[0] = outer_edges_[2*(i-1)];
            edge_nodes[1] = outer_edges_[2*(i-1)+1];
            return edge_nodes;
        }

        std::vector<d_type> inner_edge(int i) const {
            std::vector<d_type> edge_nodes(2);
            edge_nodes[0] = inner_edges_[2*(i-1)];
            edge_nodes[1] = inner_edges_[2*(i-1)+1];
            return edge_nodes;
        }

        int nb_elements() const {
            return nb_elements_;
        }

        int nb_nodes() const {
            return nb_nodes_;
        }

        int nb_outer_edges() const {
            return nb_outer_edges_;
        }

        int nb_inner_edges() const {
            return nb_inner_edges_;
        }


    private:
        std::string file_name_;
        std::vector<d_type> nodes_; //  length 2N_: x1 y1 x2 y2 x3 y3 for sequential mem access
        std::vector<d_type> inner_edges_;
        std::vector<d_type> elements_;
        std::vector<d_type> outer_edges_;
        int nb_nodes_;
        int nb_elements_;
        int nb_outer_edges_;
        int nb_inner_edges_;

    };



}




#endif