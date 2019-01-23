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

#include <string>
#include <iostream>

#include "thrust/device_vector.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <algorithm>

#include <chrono>
#include<time.h>
#include <sys/time.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include"fstream"
#include "iostream"
#include "graphHOST.h"
#include "graphGPU.h"
#include "communityGPU.h"
#include"list"


int main(int argc, char** argv) {


	char* file_w = NULL;
	int type = UNWEIGHTED;

	ofstream logFile;
	string logFileName = "Log/louvain_method_gpu_runtime_and_modularity.csv";
	ifstream infile(logFileName);
	bool existing_file = infile.good();
	logFile.open(logFileName, ios_base::out | ios_base::app | ios_base::ate);
	if (!existing_file) {
		logFile << "GraphName" << "," << "Total Time" << "," << "Modularity" << std::endl;
	}


	std::cout << "#Args: " << argc << std::endl;
	for (int i = 0; i < argc; i++) {
		std::cout << i << " : " << argv[i] << std::endl;
	}

	if (argc == 3) {
		file_w = argv[2];
		type = WEIGHTED;
		if (type == WEIGHTED)
			std::cout << "Weighted Graph \n";
	}

	if (file_w)
		std::cout << "inputGraph: " << argv[1] << " Corresponding Weight: " << file_w << std::endl;
	else if (argc==2)
		std::cout << "inputGraph: " << argv[1] << std::endl;
	else 
		std::cout<<"No input graph provided, creating a sample graph"<<std::endl;

	// Read Graph in  host memory
	//GraphHOST input_graph(argv[1], file_w, type);

	//Create a graph in host memory
	GraphHOST input_graph; // Sample graph
	if (1) {

		input_graph.nb_nodes = 7;
		input_graph.degrees.resize(input_graph.nb_nodes);
		int statIndices[] = {3, 4, 5, 8, 10, 13, 14};
		std::copy(statIndices, statIndices + (input_graph.nb_nodes), input_graph.degrees.begin());

		input_graph.nb_links = 14;
		input_graph.links.resize(input_graph.nb_links);
		unsigned int edges[] = {1, 2, 3, 0, 0, 0, 4, 5, 3, 5, 3, 4, 6, 5};
		std::copy(edges, edges + input_graph.nb_links, input_graph.links.begin());

	}

	input_graph.display();

	double threshold = 0.000001;
	if(argc==4) threshold = atof(argv[2]);
	double binThreshold = 0.01;
	if(argc==4) binThreshold=atof(argv[3]);
	//binThreshold=threshold;
	//Copy Graph to Device
	Community dev_community(input_graph, -1, threshold);
	double cur_mod = -1.0, prev_mod = 1.0;
	bool improvement = false;

	std::cout << "threshold: " << threshold << " binThreshold: " << binThreshold << std::endl;

	//Read Prime numbers
	dev_community.readPrimes("fewprimes.txt");

	cudaStream_t *streams = NULL;
	int n_streams = 8;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//clock_t clk_decision, clk_contraction;
	std::vector<clock_t> clkList_decision;
	std::vector<clock_t> clkList_contration;

	clock_t t1, t2, t3;
	t1 = clock();

	time_t time_begin, time_end;
	time(&time_begin);
	struct timespec start_comm, end_comm;
	clock_gettime(CLOCK_MONOTONIC, &start_comm);

	/*			
				dev_community.preProcess();
				dev_community.gatherStatistics(true); // MUST be "true" to filter out isolated vertices at the beginning
				dev_community.compute_next_graph();
				dev_community.set_new_graph_as_current();
	 */
	bool TEPS = true;
	bool islastRound = false;
	int szSmallComm = 100000;
	bool isGauss =true;// false;

	if(isGauss)
		std::cout<<"\n Update method:  Gauss–Seidel (in batch) \n";
	else
		std::cout<<"\n Update method: Jacobi\n";

	int stepID = 1;

	do {

		std::cout << "---------------Calling method for modularity optimization------------- \n";
		t2 = clock();
		prev_mod = cur_mod;

		cur_mod = dev_community.one_levelGaussSeidel(cur_mod, islastRound,
				szSmallComm, binThreshold, isGauss &&(dev_community.community_size > szSmallComm),
				streams, n_streams, start, stop);

		t2 = clock() - t2;

		clkList_decision.push_back(t2); // push the clock for the decision

		std::cout<< "step: " <<stepID <<", Time for modularity optimization: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;
		stepID++;
		if (TEPS == true) {
			std::cout<<binThreshold<<"_"<<threshold<< " #E:" <<  dev_community.g.nb_links << "  TEPS: " << dev_community.g.nb_links / (((float) t2) / CLOCKS_PER_SEC) << std::endl;
			TEPS = false;
		}

		//break;
		std::cout << "Computed modularity: " << cur_mod << " ( init_mod = " << prev_mod << " ) " << std::endl;

		if ((cur_mod - prev_mod) > threshold) {

			t2 = clock();   
			t3 = t2;
			dev_community.gatherStatistics();
			t2 = clock() - t2;
			//std::cout << "T_gatherStatistics: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;

			t2 = clock();
			dev_community.compute_next_graph(streams, n_streams, start, stop);
			t2 = clock() - t2;
			std::cout << "Time to compute next graph: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;

			//clkList_contration.push_back(t2); // push the clock for the contraction

			t2 = clock();
			dev_community.set_new_graph_as_current();
			t2 = clock() - t2;
			t3 = clock() -t3;

			clkList_contration.push_back(t3); // push the clock for the contraction

			//std::cout << "T_new_graph_as_current: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;
			// break;

			//std::cout << "\n Back to Main \n";
			//int sc;
			//std::cin>>sc;


		} else {
			if (islastRound == false) {
				islastRound = true;
			} else {
				break;
			}
		}
		//improvement = false;
	} while (true);

	std::cout<< "#phase: "<<stepID<<std::endl;

	clock_gettime(CLOCK_MONOTONIC, &end_comm);
	double elapsed_time = ((end_comm.tv_sec*1000 + (end_comm.tv_nsec/1.0e6)) - (start_comm.tv_sec*1000 + (start_comm.tv_nsec/1.0e6)));

	time(&time_end);
	logFile<<argv[1]<<","<<elapsed_time<<","<<prev_mod<<std::endl;

	t2 = clock();
	float diff = ((float) t2 - (float) t1);
	float seconds = diff / CLOCKS_PER_SEC;

	if( argc ==1){
		std::cout <<  binThreshold<<"_"<<threshold<<" Running Time: " << seconds << " ;  Final Modularity: "
			<< prev_mod  << std::endl;
	}else{

		std::cout <<  binThreshold<<"_"<<threshold<<" Running Time: " << seconds << " ;  Final Modularity: "
			<< prev_mod << " inputGraph: " << argv[1] << std::endl;
	}


	std::cout << "#Record(clk_optimization): " << clkList_decision.size()
		<< " #Record(clk_contraction):" << clkList_contration.size() << std::endl;

	int nrPhase = std::min(clkList_decision.size(), clkList_contration.size());

	std::ofstream ofs ("time.txt", std::ofstream::out);

	//------------------------------------To Plot------------------------------------//
	/*  
	    int nrphaseToPlot = 40;
	    for (int i = 0; i < std::min(nrphaseToPlot, nrPhase); i++) {
	    std::cout << (float) clkList_decision[i] / CLOCKS_PER_SEC << " ";
	    }
	    std::cout << std::endl;

	    for (int i = 0; i < std::min(nrphaseToPlot, nrPhase); i++) {
	    std::cout << (float) clkList_contration[i] / CLOCKS_PER_SEC << " ";
	    }
	    std::cout << std::endl;
	 */


	//----------------------------------------------------//

	float t_decision = 0, t_contraction = 0;

	for (int i = 0; i < clkList_decision.size(); i++) {

		t_decision += (float) clkList_decision[i] / CLOCKS_PER_SEC;
		if(i<nrPhase) ofs<< (float) clkList_decision[i] / CLOCKS_PER_SEC<<" ";
		else std::cout<<  (float) clkList_decision[i] / CLOCKS_PER_SEC<<" -> "<<std::endl;

	}

	ofs<<"\n";

	for (int i = 0; i < clkList_contration.size(); i++) {
		t_contraction += (float) clkList_contration[i] / CLOCKS_PER_SEC;
		if(i<nrPhase) ofs<< (float) clkList_contration[i] / CLOCKS_PER_SEC<<" ";

	}

	ofs<<"\n";
	ofs.close();

	std::cout<< " Optimization and contraction time  ratio:"
		<< (100 * t_decision)/(t_decision + t_contraction) << " " << (100 * t_contraction)/(t_decision+t_contraction) << std::endl;


	/*
	   for (int i = nrPhase; i < clkList_decision.size(); i++) {

	   std::cout << clkList_decision[i] << " "
	   << (float) clkList_decision[i] / CLOCKS_PER_SEC << std::endl;
	   }
	 */

	std::cout << "(graph):      #V  " << dev_community.g.nb_nodes << " #E   " << dev_community.g.nb_links << std::endl;
	std::cout << "(new graph)  #V  " << dev_community.g_next.nb_nodes << " #E  " << dev_community.g_next.nb_links << std::endl;
	/* 
	   for (int i = 0; i < n_streams; i++) {
	//CHECK(cudaStreamDestroy(streams[i]));
	}
	 */

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//free(streams);

	return 0;
}


