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
#include "eigen/Eigen/Eigen"


namespace pear {

    template <typename d_type, typename mat_type>
    class grid{

        typedef Eigen::Triplet<d_type> T;

    public:

        /* Constructor for grid
         *
         * Loads the mesh description and stores them in the datastructures compatible with later C++ use
         *
         * IN:      file directory with .txt files, in the format of the input described below
         * OUT:     - 4 #vector#'s describing the nodes, edges and elements of the mesh
         *          - 4 #scalar#'s describing the length of each of the above vectors
         *          - 1 #TripletList# describing the sparsity pattern of the future matrices
         *
         * NODES : _Nodes.txt
         * (int) grid_number //space// (d_type) x_coordinate //space// (d_type) y_coordinate
         *
         * ELEMENTS : _Elements.txt
         * (int) element_number //space// (int) first_node //space// (int) second_node //space// (int) third_node
         *
         * OUTER EDGES : _OuterEdges.txt
         * (int) edge_number //space// (int) first_node //space// (int) second_node
         *
         * INNER EDGES : _InnerEdges.txt
         * (int) edge_number //space// (int) first_node //space// (int) second_node
         *
         */
        grid(std::string file_name)
        : file_name_(file_name)
        {

            std::cout<<"//GRID//"<<std::endl;

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

                    // Form sparsity pattern
                    n -= 1;
                    // Upper left block
                    tripletList_.push_back(T(n,n,0.0));

                    // Lower left block
                    tripletList_.push_back(T(n+nb_nodes_,n,0.0));

                    // Upper right block
                    tripletList_.push_back(T(n,n+nb_nodes_,0.0));

                    // Lower right block
                    tripletList_.push_back(T(n+nb_nodes_,n+nb_nodes_,0.0));

                }
                std::cout<<"        "<<nb_nodes_<<" nodes loaded from: "<<file_name_<<std::endl;
                nodes_file.close();

            }
            else {std::cout << "Unable to open file, check the executable folder";}

            // ELEMENTS

            std::ifstream elements_file(file_name_ + "_Elements.txt");
            if (elements_file.is_open()){
                // Count number of elements to pre-allocate memory
                N = 0; while(elements_file>>n>>n1>>n2>>n3){N++;}; nb_elements_ = N;
                elements_ = std::vector<int>(3*nb_elements_);

                // Reopen
                elements_file.close();
                std::ifstream elements_file(file_name_+"_Elements.txt");

                while(elements_file>>n>>n1>>n2>>n3){
                    // Save coordinates
                    elements_[3*(n-1)] = n1;
                    elements_[3*(n-1)+1] = n2;
                    elements_[3*(n-1)+2] = n3;

                    // Form sparsity pattern
                    n1 -= 1;
                    n2 -= 1;
                    n3 -= 1;

                    // Upper left block
                    tripletList_.push_back(T(n1,n2,0.0));
                    tripletList_.push_back(T(n2,n1,0.0));
                    tripletList_.push_back(T(n1,n3,0.0));
                    tripletList_.push_back(T(n3,n1,0.0));
                    tripletList_.push_back(T(n2,n3,0.0));
                    tripletList_.push_back(T(n3,n2,0.0));
                    tripletList_.push_back(T(n1,n1,0.0));
                    tripletList_.push_back(T(n2,n2,0.0));
                    tripletList_.push_back(T(n3,n3,0.0));


                    // Lower left block
                    tripletList_.push_back(T(n1+nb_nodes_,n2,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n1,0.0));
                    tripletList_.push_back(T(n1+nb_nodes_,n3,0.0));
                    tripletList_.push_back(T(n3+nb_nodes_,n1,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n3,0.0));
                    tripletList_.push_back(T(n3+nb_nodes_,n2,0.0));
                    tripletList_.push_back(T(n1+nb_nodes_,n1,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n2,0.0));
                    tripletList_.push_back(T(n3+nb_nodes_,n3,0.0));


                    // Upper right block
                    tripletList_.push_back(T(n1,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2,n1+nb_nodes_,0.0));
                    tripletList_.push_back(T(n1,n3+nb_nodes_,0.0));
                    tripletList_.push_back(T(n3,n1+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2,n3+nb_nodes_,0.0));
                    tripletList_.push_back(T(n3,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n1,n1+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n3,n3+nb_nodes_,0.0));

                    // Lower right block
                    tripletList_.push_back(T(n1+nb_nodes_,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n1+nb_nodes_,0.0));
                    tripletList_.push_back(T(n1+nb_nodes_,n3+nb_nodes_,0.0));
                    tripletList_.push_back(T(n3+nb_nodes_,n1+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n3+nb_nodes_,0.0));
                    tripletList_.push_back(T(n3+nb_nodes_,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n1+nb_nodes_,n1+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n3+nb_nodes_,n3+nb_nodes_,0.0));

                }
                std::cout<<"        "<<nb_elements_<<" elements loaded from: "<<file_name_<<std::endl;
                elements_file.close();
            }
            else {std::cout << "Unable to open file, check the executable folder";}

            // OUTER EDGES

            std::ifstream outer_file(file_name_ + "_OuterEdges.txt");
            if (outer_file.is_open()){
                // Count number of outer edges to pre-allocate memory
                N = 0; while(outer_file>>n>>n1>>n2){N++;}; nb_outer_edges_ = N;
                outer_edges_ = std::vector<int>(2*nb_outer_edges_);
                // Reopen
                outer_file.close();
                std::ifstream outer_file(file_name_+"_OuterEdges.txt");
                // Save coordinates
                while(outer_file>>n>>n1>>n2){
                    outer_edges_[2*(n-1)] = n1;
                    outer_edges_[2*(n-1)+1] = n2;

                    n1 -= 1;
                    n2 -= 1;

                    // Upper left block
                    tripletList_.push_back(T(n1,n2,0.0));
                    tripletList_.push_back(T(n2,n1,0.0));

                    // Lower left block
                    tripletList_.push_back(T(n1+nb_nodes_,n2,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n1,0.0));



                    // Upper right block
                    tripletList_.push_back(T(n1,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2,n1+nb_nodes_,0.0));


                    // Lower right block
                    tripletList_.push_back(T(n1+nb_nodes_,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n1+nb_nodes_,0.0));



                }
                std::cout<<"        "<<nb_outer_edges_<<" outer edges loaded from: "<<file_name_<<std::endl;
                outer_file.close();
            }

            else {std::cout << "Unable to open file, check the executable folder";}

            // INNER EDGES

            std::ifstream inner_file(file_name_ + "_InnerEdges.txt"); // read file
            if (inner_file.is_open()){
                // Count number of inner edges to pre-allocate memory
                N = 0; while(inner_file>>n>>n1>>n2){N++;}; nb_inner_edges_ = N;
                inner_edges_ = std::vector<int>(2*nb_inner_edges_);
                // Reopen
                inner_file.close();
                std::ifstream inner_file(file_name_+"_InnerEdges.txt");
                // Save coordinates
                while(inner_file>>n>>n1>>n2){
                    inner_edges_[2*(n-1)] = n1;
                    inner_edges_[2*(n-1)+1] = n2;

                    n1 -= 1;
                    n2 -= 1;

                    // Upper left block
                    tripletList_.push_back(T(n1,n2,0.0));
                    tripletList_.push_back(T(n2,n1,0.0));

                    // Lower left block
                    tripletList_.push_back(T(n1+nb_nodes_,n2,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n1,0.0));


                    // Upper right block
                    tripletList_.push_back(T(n1,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2,n1+nb_nodes_,0.0));


                    // Lower right block
                    tripletList_.push_back(T(n1+nb_nodes_,n2+nb_nodes_,0.0));
                    tripletList_.push_back(T(n2+nb_nodes_,n1+nb_nodes_,0.0));
                }
                std::cout<<"        "<<nb_inner_edges_<<" inner edges loaded from: "<<file_name_<<std::endl;
                inner_file.close();
            }
            else {std::cout << "Unable to open file, check the executable folder";};
        }

        /* Coordinates of each node
         *
         * IN:      #int# describing the number of the node
         * OUT:     #vector<double># (of length 2X1) describing the (x, y) coordinates of the node
         *
         */
        std::vector<d_type> node(int i) const{
            std::vector<d_type> node_coord(2);
            node_coord[0] = nodes_[2*(i-1)];
            node_coord[1] = nodes_[2*(i-1)+1];
            return node_coord;
        }

        /* Elements associated to each node
         *
         * IN:      #int# describing the number of the node
         * OUT:     #vector<int># (of variable length NX1) describing the elements to which the node belongs
         *
         */
        std::vector<int> elements_for_node(int node) const{
            std::vector<int> edge_nodes(15);         // Pre-allocation to upper-limit

            int N = 0;                                  // Number of elements found
            for (int i = 1; i < nb_elements_ + 1 ; i++) {   // Loop on all the elements

                if (elements_[3*(i-1)] == node || elements_[3*(i-1)+1] == node || elements_[3*(i-1)+2] == node){
                    ++N;                                                            // Update count
                    if ( N > static_cast<signed int>(edge_nodes.size()) ) { edge_nodes.resize(edge_nodes.size()*2); }        // Resize if necessary
                    edge_nodes[N-1] = i;
                }

            }
            edge_nodes.resize(N);                   // Final resize
            return edge_nodes;
        }

        /* Nodes of each element
         *
         * IN:      #int# describing the number of the element
         * OUT:     #vector<int># (of length 3X1) describing the numbers of the nodes/vertices of each element
         *
         */
        std::vector<int> element(int i) const{
            std::vector<int> elem_nodes(3);
            elem_nodes[0] = elements_[3*(i-1)];
            elem_nodes[1] = elements_[3*(i-1)+1];
            elem_nodes[2] = elements_[3*(i-1)+2];
            return elem_nodes;
        }

        /* Nodes of each outer edge
         *
         * IN:      #int# describing the number of the outer edge
         * OUT:     #vector<int># (of length 2X1) describing the numbers of the nodes/vertices of each edge
         *
         */
        std::vector<int> outer_edge(int i) const{
            std::vector<int> edge_nodes(2);
            edge_nodes[0] = outer_edges_[2*(i-1)];
            edge_nodes[1] = outer_edges_[2*(i-1)+1];
            return edge_nodes;
        }

        /* Nodes of each inner edge
         *
         * IN:      #int# describing the number of the inner edge
         * OUT:     #vector<int># (of length 2X1) describing the numbers of the nodes/vertices of each edge
         *
         */
        std::vector<int> inner_edge(int i) const{
            std::vector<int> edge_nodes(2);
            edge_nodes[0] = inner_edges_[2*(i-1)];
            edge_nodes[1] = inner_edges_[2*(i-1)+1];
            return edge_nodes;
        }

        int nb_elements() const {       return nb_elements_;}
        int nb_nodes() const {          return nb_nodes_;}
        int nb_outer_edges() const {    return nb_outer_edges_;}
        int nb_inner_edges() const {    return nb_inner_edges_;}

        /* Sparsity pattern
         *
         * IN:      #sparse_matrix# to be patterned
         * OUT:     void. Upon output, the #sparse_matrix# has been allocated the entries later required;
         *              - the diagonal (each node to itself)
         *              - neighbouring nodes in the same element
         *              - neighbouring nodes on the same edge
         *              Warning: the sparsity pattern depends on the quadrature rule
         *
         */
        void setSparsityPattern(mat_type & m) const {
            m.setFromTriplets(tripletList_.begin(), tripletList_.end());
        }

    private:
        std::string file_name_;
        std::vector<d_type> nodes_; //  length 2N_: x1 y1 x2 y2 x3 y3 for sequential mem access
        std::vector<int> inner_edges_;
        std::vector<int> elements_;
        std::vector<int> outer_edges_;
        std::vector<T> tripletList_;
        int nb_nodes_;
        int nb_elements_;
        int nb_outer_edges_;
        int nb_inner_edges_;

    };



}




#endif