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
#include"hostconstants.h"

void Community::preProcess() {



    //Compute degree of each node
    thrust::device_vector<int> sizesOfNhoods(g.indices.size() - 1, 0);


    thrust::transform(g.indices.begin() + 1, g.indices.end(), g.indices.begin(),
            sizesOfNhoods.begin(), thrust::minus<int >());

    //Find all degree 1 vertices
    IsInRange<int, int> filter_SNL_1(1, 1);

    int nrC_SNL_1 = thrust::count_if(thrust::device, sizesOfNhoods.begin(),
            sizesOfNhoods.end(), filter_SNL_1);

    std::cout << "#vertices (SNL=1):" << nrC_SNL_1 << std::endl;

    //Lets copy Identities of all communities  in g_next.links

    g_next.links.resize(community_size, 0);
    thrust::sequence(g_next.links.begin(), g_next.links.end(), 0);

    //Use g_next.indices to copy community ids with  SLN =1
    g_next.indices.resize(community_size, -1);


    //Collet all  degree 1 vertices in g_next.indices
    thrust::copy_if(thrust::device, g_next.links.begin(), g_next.links.end(),
            sizesOfNhoods.begin(), g_next.indices.begin(), filter_SNL_1);



    /*
    void reduceGraph(int* indices, unsigned int* links, float* weights, int gType,
            int* uniDegvrts, unsigned int nrUniDegVrts, unsigned int mark,
            int* vtsForPostProcessing);
     */


    int mark = g.nb_nodes * 2;
    thrust::device_vector<int> vtsForPostProcessing(nrC_SNL_1, -1);

    unsigned int nrBlk = (nrC_SNL_1 + NR_THREAD_PER_BLOCK - 1) / NR_THREAD_PER_BLOCK;

    //initialization of n2c
    n2c.resize(community_size);
    thrust::sequence(n2c.begin(), n2c.end(), 0);
    if(nrC_SNL_1>0)
    reduceGraph << <nrBlk, NR_THREAD_PER_BLOCK>>>(
            thrust::raw_pointer_cast(g.indices.data()),
            thrust::raw_pointer_cast(g.links.data()),
            thrust::raw_pointer_cast(g.weights.data()), g.type,
            thrust::raw_pointer_cast(g_next.indices.data()), nrC_SNL_1, mark,
            thrust::raw_pointer_cast(vtsForPostProcessing.data()),
            thrust::raw_pointer_cast(n2c.data()));

    /*
    // let's process each vertex in  vtsForPostProcessing with a warp
    nrBlk = (nrC_SNL_1 + (NR_THREAD_PER_BLOCK / PHY_WRP_SZ) - 1) / (NR_THREAD_PER_BLOCK / PHY_WRP_SZ);

    void editEdgeList(int* indices, unsigned int* links, float* weights, int gType,
            int* uniDegvrts, unsigned int nrUniDegVrts, unsigned int mark,
            int* vtsForPostProcessing);

    editEdgeList << <nrBlk, NR_THREAD_PER_BLOCK>>>(
            thrust::raw_pointer_cast(g.indices.data()),
            thrust::raw_pointer_cast(g.links.data()),
            thrust::raw_pointer_cast(g.weights.data()), g.type,
            thrust::raw_pointer_cast(g_next.indices.data()), nrC_SNL_1, mark,
            thrust::raw_pointer_cast(vtsForPostProcessing.data()));
     */


    if (0) {
        thrust::host_vector<int> uniDegVertices = g_next.indices;
        std::cout << std::endl;
        for (int i = 0; i < nrC_SNL_1; i++) {
            std::cout << uniDegVertices[i] << " ";
        }
        std::cout << std::endl;
    }

    if (0) {
        thrust::host_vector<unsigned int> gnlinks = g.links;
        thrust::host_vector<int> gnIndices = g.indices;
        thrust::device_vector<float> gnWeights = g.weights;
        for (unsigned int i = 0; i < g.nb_nodes; i++) {

            unsigned int startNbr = gnIndices[i];
            unsigned int endNbr = gnIndices[i + 1];
            //thrust::sort(gnlinks.begin() + startNbr, gnlinks.begin() + endNbr);

           if(i<10){
	    std::cout << i << ":";

            for (unsigned int j = startNbr; j < endNbr; j++) {

                float edgeWt = 1;
                if (g.type == WEIGHTED)
                    edgeWt = gnWeights[j];
                std::cout << " " << gnlinks[j] << "(" << edgeWt << ")";
            }
            std::cout << std::endl;
	    }
        }
    }

    n2c.clear();
    vtsForPostProcessing.clear();
    sizesOfNhoods.clear();
    g_next.links.clear();
    g_next.indices.clear();

}
