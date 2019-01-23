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

#include <fstream>
#include <algorithm>
#include <iostream>
#include "communityGPU.h"

void Community::set_new_graph_as_current() {

    g.indices.clear();
    g.links.clear();
    g.weights.clear();

    g.indices.resize(g_next.nb_nodes + 1);
    g.links.resize(g_next.nb_links);
    g.weights.resize(g_next.nb_links);

    /*
    {
        std::cout << std::endl << "[Before Exchange] g.size(): " << g.nb_nodes << std::endl;
    }
     */
    g = g_next;

    //deep copy or sallow copy ?
    {

        g_next.indices.clear();
        g_next.links.clear();
        g_next.weights.clear();
    }
    /*
    {
        std::cout << std::endl << "[After Exchange] g.size(): " << g.nb_nodes << " #E: " << g.links.size() << std::endl;
    }
     */

    //std::cout << "...............Set New Graph as Current.................................." << std::endl;

    int sc;
    sc = 0; //std::cin>>sc;
    if (sc > 0) {
        print_vector(g.indices, "g.indices( after exchanging): ");
        print_vector(g.weights, "g.weights( after exchanging): ");
        print_vector(g.links, "g.links( after exchanging): ");
    }

    community_size = g.nb_nodes;

}

void Community::readPrimes(std::string filename) {

    std::ifstream filePtr(filename.c_str());

    if (filePtr.is_open()) {

        int nrPrimes, aPrimeNum;

        filePtr>>nrPrimes;

        assert(nrPrimes > 0);

        std::cout << "Reading " << nrPrimes << " prime numbers." << std::endl;

        //Read primes in host memory
        hostPrimes = new int [nrPrimes];
        int index = 0;
        while (filePtr >> aPrimeNum) {

            hostPrimes[index++] = aPrimeNum;
            if (index >= nrPrimes)
                break;
            //std::cout << aPrimeNum << " ";
        }

        std::cout << std::endl;

        assert(nrPrimes == index);
        nb_prime = nrPrimes;

        //Copy prime numbers to device memory
        devPrimes.resize(nb_prime);
        thrust::copy(hostPrimes, hostPrimes + nb_prime, devPrimes.begin());

    } else {
        std::cout << "Can't open file containing prime numbers." << std::endl;
    }
    return;
}

