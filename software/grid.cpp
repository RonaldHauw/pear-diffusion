#ifndef GRID_CPP
#define GRID_CPP

#include "grid.hpp"


namespace pear {

    template <typename d_type, typename mat_type>
    std::vector<d_type> grid<d_type, mat_type>::node(int i) const {
        std::vector<d_type> node_coord(2);
        node_coord[0] = nodes_[2*(i-1)];
        node_coord[1] = nodes_[2*(i-1)+1];
        return node_coord;
    }

    template <typename d_type, typename mat_type>
    std::vector<int> grid<d_type, mat_type>::elements_for_node(int node) const {
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

    template <typename d_type, typename mat_type>
    std::vector<int> grid<d_type, mat_type>::element(int i) const {
            std::vector<int> elem_nodes(3);
            elem_nodes[0] = elements_[3*(i-1)];
            elem_nodes[1] = elements_[3*(i-1)+1];
            elem_nodes[2] = elements_[3*(i-1)+2];
            return elem_nodes;
    }

    template <typename d_type, typename mat_type>
    std::vector<int> grid<d_type, mat_type>::outer_edge(int i) const {
            std::vector<int> edge_nodes(2);
            edge_nodes[0] = outer_edges_[2*(i-1)];
            edge_nodes[1] = outer_edges_[2*(i-1)+1];
            return edge_nodes;
    }


    template <typename d_type, typename mat_type>
    std::vector<int> grid<d_type, mat_type>::inner_edge(int i) const {
            std::vector<int> edge_nodes(2);
            edge_nodes[0] = inner_edges_[2*(i-1)];
            edge_nodes[1] = inner_edges_[2*(i-1)+1];
            return edge_nodes;
    }

}




#endif