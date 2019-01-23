#include"graphHOST.h"
#include"fstream"
#include "iostream"
#include"vector"
using namespace std;

GraphHOST::GraphHOST(char* filename, char* filename_w, int type) {

    ifstream finput;
    finput.open(filename, fstream::in | fstream::binary);

    // Read number of nodes on 4 bytes
    finput.read((char *) &nb_nodes, 4);
    assert(finput.rdstate() == ios::goodbit);

    // Read cumulative degree sequence: 8 bytes for each node
    // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
    degrees.resize(nb_nodes);
    finput.read((char *) &degrees[0], nb_nodes * 8);

    // Read links: 4 bytes for each link (each link is counted twice)
    nb_links = degrees[nb_nodes - 1];
    links.resize(nb_links);
    finput.read((char *) (&links[0]), (long) nb_links * 4);

    // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
    weights.resize(0);

    if (type == WEIGHTED) {
        ifstream finput_w;
        finput_w.open(filename_w, fstream::in | fstream::binary);
        weights.resize(nb_links);
        finput_w.read((char *) &weights[0], (long) nb_links * 4);

    }

    // Compute total weight
    total_weight = 0;
    for (unsigned int i = 0; i < nb_nodes; i++) {
        total_weight += (double) weighted_degree(i);
    }
    if (type == UNWEIGHTED) {
        std::cout << std::endl << "UNWEIGHTED" << std::endl;
    } else {
        std::cout << std::endl << "WEIGHTED" << std::endl;
    }
    std::cout << " total_weight = " << total_weight << std::endl;
}

GraphHOST::GraphHOST() {
    nb_nodes = 0;
    nb_links = 0;
    total_weight = 0;
}

void
GraphHOST::display() {



    for (unsigned int node = 0; node < nb_nodes; node++) {
        pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
        //thrust::pair<thrust::host_vector<unsigned int>::iterator, thrust::host_vector<float>::iterator > p = neighbors(node);
        /*
          if (node >= 31220 && node <= 31223)
             cout << "node: " << node << " : nr_neighbor: " << nb_neighbors(node) << endl;

         if (node == 34693)
             cout << "node: " << node << " : nr_neighbor: " << nb_neighbors(node) << endl;
         */
        cout << node << " : ";
        for (unsigned int i = 0; i < nb_neighbors(node); i++) {
            if (true) {
                if (weights.size() != 0)
                    cout << " (" << *(p.first + i) << " " << *(p.second + i) << ")";
                else
                    cout << " " << *(p.first + i);
            }
        }
        cout << endl;

    }
}
