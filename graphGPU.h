
#ifndef GRAPHGPU_H
#define	GRAPHGPU_H

#include "thrust/device_vector.h"
#include "iostream"
#include"commonconstants.h"

struct GraphGPU {
    // Graph 
    unsigned int nb_nodes;
    unsigned long nb_links;
    double total_weight;

    int type;

    //thrust::device_vector<unsigned long> degrees;

    thrust::device_vector<int> indices;

    thrust::device_vector<unsigned int> links;
    thrust::device_vector<float> weights;
    thrust::device_vector<int> colors;
    void greedyColoring(unsigned int wrpSz);

    //unsigned int nb_neighbors(unsigned int node);
    //double weighted_degree(unsigned int node);


};


#endif	/* GRAPHGPU_H */

