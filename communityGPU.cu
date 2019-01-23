/*

    Copyright (C) 2016, University of Bergen

    This file is part of Rundemanen - CUDA C++ parallel program for
    community detection

    Rundemanen is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rundemanen is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rundemanen.  If not, see <http://www.gnu.org/licenses/>.
    
    */

#include <algorithm>
#include <iostream>
#include "communityGPU.h"

#include <functional>
#include"numeric"

Community::Community(const GraphHOST& input_graph, int nb_pass, double min_mod) {


    //Graph
    g.nb_links = input_graph.nb_links;
    g.nb_nodes = input_graph.nb_nodes;
    g.type = UNWEIGHTED;



    //Copy degree array into indices with an extra zero(0) at the beginning
    g.indices = thrust::device_vector<int>(input_graph.nb_nodes + 1, 0);
    thrust::copy(input_graph.degrees.begin(), input_graph.degrees.end(), g.indices.begin() + 1); // 0 at first position

    /********************Gather Graph Statistics***************/
    std::vector< int> vtxDegs;

    vtxDegs.resize(input_graph.degrees.size());

    std::adjacent_difference(input_graph.degrees.begin(), input_graph.degrees.end(), vtxDegs.begin());

    int totNbrs = std::accumulate(vtxDegs.begin(), vtxDegs.end(), 0);
    int maxDeg = *std::max_element(vtxDegs.begin(), vtxDegs.end());

    double sumSquareDiff = 0;
    double avgDeg = (double) totNbrs / g.nb_nodes;

    for (int i = 0; i < vtxDegs.size(); i++) {
        double delta = ((double) vtxDegs[i] - avgDeg);
        sumSquareDiff += delta*delta;
    }

    double standardDeviation = sqrt(sumSquareDiff / input_graph.nb_nodes);

    std::cout << "MaxDeg = " << maxDeg << " AvgDeg = " << avgDeg << " STD = "
            << standardDeviation << " STD2AvgRatio = " << standardDeviation / avgDeg << std::endl;

    std::cout << "totNbrs =" << totNbrs << " #links =" << input_graph.nb_links << std::endl;

    if (input_graph.nb_nodes < 10) {
        std::cout << std::endl;
        for (int i = 0; i < vtxDegs.size(); i++) {
            std::cout << vtxDegs[i] << " ";
        }
        std::cout << std::endl;
    }
    /**********************************/

    //copy all edges
    g.links.resize(g.nb_links);
    g.links = input_graph.links;

    //copy all weights
    g.weights.resize(input_graph.weights.size());
    g.weights = input_graph.weights;

    std::cout << std::endl << "Copied  " << g.weights.size() << " weights" << std::endl;

    g.total_weight = input_graph.total_weight;


    if (input_graph.weights.size() > 0) {
        g.type = WEIGHTED;
        std::cout << " Setting type to WEIGHTED" << std::endl;
    } else {
        std::cout << "Type is already set to UNWEIGHTED" << std::endl;
    }

    //Community
    community_size = g.nb_nodes;
    min_modularity = min_mod;

    std::cout << std::endl << "(Dev Graph) " << " #Nodes: " << g.nb_nodes << "  #Links: " << g.nb_links / 2 << "  Total_Weight: " << g.total_weight / 2 << std::endl;
    std::cout << "community_size: " << community_size << std::endl;
    // seriously !!
}
