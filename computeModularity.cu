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

struct my_modularity_functor_2 {
    double m2;

#ifdef RUNONGPU

    __host__ __device__
#endif

    my_modularity_functor_2(double _m2) : m2(_m2) {
    }

#ifdef RUNONGPU

    __host__ __device__
#endif

    double operator()(const float& x, const float& y) {
        return ((double) x / m2 - ((double) y / m2)*((double) y / m2));
    }

};

double Community::modularity(thrust::device_vector<float>& tot, thrust::device_vector<float>& in) { // put i=j in equation (1)

    float q = 0.;
    float m2 = (float) g.total_weight;

    std::cout << "m2: " << m2 << std::endl;

    /*
    thrust::host_vector<float> in_ = in;
    thrust::host_vector<float> tot_ = tot;
    //std::cout << "QQ:" << std::endl;
    for (int i = 0; i < community_size; i++) {
        if (tot_[i] > 0) {
            //std::cout << in_[i] << " " << tot_[i] << " " << m2 << " : ";
            q += in_[i] / m2 - (tot_[i] / m2)*(tot_[i] / m2);
        }
    }
    //std::cout << std::endl;
    return q;
     */

    /*
    double q = 0.;
    double m2 = (double) g.total_weight;

    for (int i = 0; i < size; i++) {
        if (tot[i] > 0)
            q += (double) in[i] / m2 - ((double) tot[i] / m2)*((double) tot[i] / m2);
    }

    return q;
     */

    bool hostPrint = false;

    if (hostPrint) {
        std::cout << std::endl << " Inside  modularity() " << std::endl;

        thrust::copy(tot.begin(), tot.end(), std::ostream_iterator<float>(std::cout, " "));
        std::cout << std::endl;

        thrust::copy(in.begin(), in.end(), std::ostream_iterator<float>(std::cout, " "));
        std::cout << std::endl;
    }

    int sc;
    //std::cout << community_size << " |in|:" << in.size() << " |tot|:" << tot.size() << std::endl;

    sc = 0; //std::cin>>sc;
    thrust::device_vector<double> result_array(community_size, 0.0);
    sc = 0; //std::cin>>sc;

    thrust::transform(thrust::device, in.begin(), in.end(), tot.begin(), result_array.begin(), my_modularity_functor_2(m2));

    q = thrust::reduce(thrust::device, result_array.begin(), result_array.end(), (double) 0, thrust::plus<double>());
    result_array.clear();
    return q;
}

