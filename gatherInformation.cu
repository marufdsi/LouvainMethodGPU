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


#include"fstream"
#include"communityGPU.h"

template < class T>
__global__
void filter_entries_by_threshold(T *source, T* dest, T threshold, int nr_old_communities, int* locations) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < nr_old_communities) {
        // If source[tid] contains something > threshold  and 
        // locations[i] is greater than zero(0)
        T data = source[tid];
        int index = locations[tid];

        if (data > threshold && index >= 0) { //zero based indexing
            //printf("\n %d %d %d\n", tid, index, data);
            dest[index] = source[tid];
        }

        tid += blockDim.x * gridDim.x;
    }
}

void Community::gatherStatistics(bool isPreprocess) {

    bool hostPrint = false;
    int sc;
    sc = 0; //std::cin>>sc;
    hostPrint = (sc > 1);

    thrust::device_vector<int> renumber(community_size, 0);

/*
    thrust::host_vector<int> hn2c = n2c;
    std::ofstream ofs;
    ofs.open ("n2c.txt", std::ofstream::out | std::ofstream::app);
    for(int i=0; i< hn2c.size(); i++) {

	ofs<<i<<":"<<hn2c[i]<<" ";

	if(hn2c[i]==320)
		std::cout<<"####" << i  <<"; (cid before renumbering) "<<hn2c[i]<<std::endl;

	if( (i+1)%100 == 0)
		ofs<<"\n";		
    }

    ofs<<"\n";
    ofs.close();
*/
    if (hostPrint) {
        print_vector(renumber, " Size of Communities: ");
        print_vector(n2c, "n2c: ");
    }

    int load_per_blk = CHUNK_PER_WARP * (NR_THREAD_PER_BLOCK / PHY_WRP_SZ);
    int nr_of_block = (community_size + load_per_blk - 1) / load_per_blk;


    //Count the size of each new community and store the sizes in renumber

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    if(isPreprocess ==false ){
    get_size_of_communities << < nr_of_block, NR_THREAD_PER_BLOCK >>>(
            thrust::raw_pointer_cast(renumber.data()), thrust::raw_pointer_cast(n2c.data()), g.nb_nodes);
    }else{
    get_size_of_communities_NEW << < nr_of_block, NR_THREAD_PER_BLOCK >>>(
            thrust::raw_pointer_cast(renumber.data()), thrust::raw_pointer_cast(n2c.data()), g.nb_nodes,thrust::raw_pointer_cast(g.indices.data()));
    }
    report_time(start, stop, "get_size_of_communities");


    if (hostPrint) {
        print_vector(renumber, "Size of Communities: (by Atomic Add) ");
    }


    n2c_new.resize(community_size, 0);// DO NOT clear it here, it will be cleared in gatherStatistics
    // 1 to indicate community of size greater or equal 1, 

    assert(renumber.size() == n2c_new.size());
  
    //int myZero = 0;
    //if(isPreprocess == true)
	//myZero = 1;

    thrust::transform(thrust::device, renumber.begin(), renumber.end(),
            n2c_new.begin(), IsGreaterThanZero<int>(0));

    //NOTE:n2c_new contains 0s and 1s only

    if (hostPrint) {
        print_vector(n2c_new, "Before Prefix Sum: ");
    }

    thrust::inclusive_scan(thrust::device, n2c_new.begin(), n2c_new.end(), n2c_new.begin());

    int new_nb_comm = n2c_new.back();

    if (hostPrint) {
        print_vector(n2c_new, "After Prefix Sum: ");
    }


    thrust::transform(thrust::device, renumber.begin(), renumber.end(),
            n2c_new.begin(), n2c_new.begin(), Community_ID_By_Prefix_Sum<int>());

    // After Transform, n2c_new contains mapping from old_CId to new_Cid
    {

        //(renumber): 0 4 0 0 0, 0 4 0 0 0,......( #member per Community)
        // (n2c_new): 0 1 1 1 1, 1 2 2 2 2,...............   

        //( n2c_new): X 0 X X X, X 1 X X X,.......( New CId, X = -1 in my case)

        //n2c_new now holds renumbered CIDs *******
    }

    /*
    thrust::host_vector<int> hn2c_new = n2c_new;
   
    	for(int i=0; i< hn2c_new.size(); i++) {
    		if(hn2c_new[i] == 63 || ( i>=320 && i<=329) || (i >=2945 && i <=2949) )
			std::cout<<"node:"<<i<<" ; new_CID "<<hn2c_new[i]<<std::endl;
	}
    */
    
    if (hostPrint) {
        print_vector(n2c_new, "CID by Prefix Sum: ");
    }

   //std::cout << "#NewCommunity: " << new_nb_comm << std::endl;

    //for next phase
    pos_ptr_of_new_comm.resize(new_nb_comm + 1, 0);



    if (hostPrint) {
        print_vector(pos_ptr_of_new_comm, "Pos ptrs: ");
    }

    // copy all non zero entries that corresponds to #nodes per communities 

    //Note: renumber still contains #nodes per community
    //n2c_new contains renumbered cIds of communities having size >=1

    cudaEventRecord(start, 0);

    filter_entries_by_threshold << < nr_of_block, NR_THREAD_PER_BLOCK >>>(
            thrust::raw_pointer_cast(renumber.data()), thrust::raw_pointer_cast(pos_ptr_of_new_comm.data()) + 1,
            (int) 0, community_size, thrust::raw_pointer_cast(n2c_new.data()));

    report_time(start, stop, "filter_entries_by_threshold");

    // copied all non zero entries that corresponds to #nodes per communities into pos_ptr_of_new_comm

    if (hostPrint) {
        print_vector(pos_ptr_of_new_comm, "Pos ptrs: ");
    }

    // after prefix sum pos_ptr points to start of each community where 
    //nodes of same community are placed consecutively

    thrust::inclusive_scan(thrust::device, pos_ptr_of_new_comm.begin(),
            pos_ptr_of_new_comm.end(), pos_ptr_of_new_comm.begin(),
            thrust::plus<int>());


    if (hostPrint) {
        print_vector(pos_ptr_of_new_comm, "Pos ptrs: ");

    }

    g_next.nb_nodes = new_nb_comm;

    renumber.clear();

}
