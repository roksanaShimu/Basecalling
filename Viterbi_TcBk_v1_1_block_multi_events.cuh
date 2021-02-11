/*
 * Viterbi_TcBk_v1_multi_events.cuh
 *
 *  Created on: Jun 20, 2017
 *      Author: roksana
 */

#ifndef VITERBI_TCBK_V1_1_BLOCK_MULTI_EVENTS_CUH_
#define VITERBI_TCBK_V1_1_BLOCK_MULTI_EVENTS_CUH_

#include<iostream>
#include<cuda.h>
#include <ctime>
#include <math.h>
#include <cuda_runtime.h>
#include <vector>
#include <thrust/device_vector.h>
#include <inttypes.h>

#define __BCALL_I__
#include "bcall-i.h"
#include "bcall.h"
#include "Kmer.cuh"
#include "Predefine_Values.cuh"
#include "Device_Memory_Allocation.cuh"
#include "Viterbi_TcBk_v1_functions.cuh"

using namespace std;


__global__ void kernel(int32_t *devPtr_pm, int32_t *devPtr_log_pr, int32_t *alpha, uint16_t *beta, int32_t *d_amax, unsigned int *d_jmax, int32_t *d_mutex, int ith_event, int64_t *d_temp_sum_pmax, int start_event, int end_event){

	unsigned id = threadIdx.x ;//+ blockIdx.x * threads_per_block;   //id is j from the serial code

	unsigned ith_tcbk=blockIdx.x;

	unsigned alpha_read=ith_tcbk*2; unsigned alpha_write=ith_tcbk*2+1;  	
	int32_t ln_pe_list[number_of_iteration];
	int32_t p_max=d_pmax[blockIdx.x];

	//if(ith_event<34)if (threadIdx.x==0)printf("blockIdx.x=%d and p_max=%d\n",blockIdx.x,p_max);
	
	//if(ith_event<34)if (threadIdx.x==0)printf("ith_tcbk= %d and d_temp_sum_pmax=%" PRIu64 "\n",ith_tcbk,d_temp_sum_pmax[ith_tcbk]);
	int number_of_events=end_event-start_event+1;
	__syncthreads();

	unsigned counter=0;
	for(int j=start_event; j<=end_event;j++){	
		for (int k=0;k<number_of_iteration;k++){
			unsigned temp_id=id+k*threads_per_block;
			if( temp_id < N_STATES){
				//TD DO: try to move alpha to shared memory so that it can save time accessing global memory
				int32_t ln_pt=dev_ln_ptransition(devPtr_log_pr, temp_id, alpha, beta, alpha_read, counter, ith_tcbk, number_of_events, ith_event);
				//if(ith_event==3 && threadIdx.x==2)printf("ith_event= %d ith_tcbk=%d temp_id=%d ln_pt=%" PRIu32 " \n", ith_event,ith_tcbk, temp_id, ln_pt);

				ln_pe_list[k] = dev_ln_pemission(devPtr_pm, p_max,  ln_pt, temp_id, counter, ith_tcbk);

				alpha[temp_id+alpha_write*N_STATES] = ln_pe_list[k];
			
				//if(ith_event==1 && k<3 && j<start_event+2) if(threadIdx.x==0) printf("ith_tcbk=%d j=%d k=%d and ln_pt=%d \n",ith_tcbk,j, k,ln_pt);
		
			}
		}

		/*if(ith_event==1){
			if(threadIdx.x==0){
				for (int k=0;k<number_of_iteration;k++){
					printf("k=%d and ln_pe=%d \n",k,ln_pe_list[k]);

				}
			}
		}*/

		
		__syncthreads();

		int32_t max_ln_pe=ln_pe_list[0];   	//in thread now
		unsigned max_id=id;		//in thread now
		for(int k=1; k<number_of_iteration; k++ ){
			if (max_ln_pe<ln_pe_list[k]){
				 max_ln_pe=ln_pe_list[k];
				 max_id=id+k*threads_per_block;
			}
		}


		/*if(ith_event==1){
			if(threadIdx.x==0){
				printf("max_ln_pe_in_thread=%d and max_temp_id_in_thread=%d \n",max_ln_pe,max_id);
			}
		}*/
		__syncthreads();



		//find max ln_pe and the corresponding id
		__shared__ int cache[threads_per_block]; __shared__ unsigned shared_max_id_available_in_thread;
		cache[threadIdx.x]=max_ln_pe;
		unsigned max_id_available_in_thread=threadIdx.x;
		__syncthreads();


		// reduction
		unsigned int i = blockDim.x/2;
		while(i != 0){
			if(threadIdx.x < i){
				if(cache[threadIdx.x]< cache[threadIdx.x + i]){
					cache[threadIdx.x]=cache[threadIdx.x + i];
					max_id_available_in_thread=threadIdx.x + i;
	
				}
			}

			__syncthreads();
			i /= 2;
		}
		__syncthreads();
	
		p_max=cache[0];

		if (threadIdx.x==0){
			d_amax[ith_tcbk]=cache[0];
			shared_max_id_available_in_thread=max_id_available_in_thread;
			d_temp_sum_pmax[ith_tcbk]=d_temp_sum_pmax[ith_tcbk]+cache[0];
		}
		__syncthreads();

		if (threadIdx.x==shared_max_id_available_in_thread){
			d_jmax[ith_tcbk]=max_id_available_in_thread;
		}
		__syncthreads();



		/*if(j>7359){
			if (threadIdx.x==0){
				printf("i=%d and amax=%d\n",j, *d_amax);
			}
		}*/
		unsigned temp=alpha_read;  
		alpha_read=alpha_write;
		alpha_write=temp;
		
		//alpha_read=!alpha_read; alpha_write=!alpha_write;
		counter++;
	}
	

}




//void preprocess_for_parallel(struct bcall_tcbk * fix, int32_t alpha[N_STATES], int32_t p_max[], int i, int32_t * devPtr_pm, int32_t * devPtr_log_pr, int32_t * dev_alpha,  uint16_t *dev_beta, uint16_t *h_beta, int32_t *h_amax, int32_t *d_amax, unsigned int *h_jmax, unsigned int *d_jmax, int32_t *d_mutex, int64_t *h_temp_sum_pmax, int64_t *d_temp_sum_pmax, const int n_tcbk  ){


void preprocess_for_parallel(vector<bcall_tcbk>& Bcall_TcBk, int32_t alpha[N_STATES], int32_t p_max[], int i, int32_t * devPtr_pm, int32_t * devPtr_log_pr, int32_t * dev_alpha,  uint16_t *dev_beta, uint16_t *h_beta, int32_t *h_amax, int32_t *d_amax, unsigned int *h_jmax, unsigned int *d_jmax, int32_t *d_mutex, int64_t *h_temp_sum_pmax, int64_t *d_temp_sum_pmax ){
	
	const int n_tcbk=Bcall_TcBk.size();

	struct fix_event * e;

	int j, k;


	int start_event=i; int end_event;
	if(start_event== (Bcall_TcBk[0].n_events-1) ){
		end_event=start_event;	
	}else if((start_event+number_of_events_per_kernel-1)<Bcall_TcBk[0].n_events){
		end_event=start_event+number_of_events_per_kernel-1;
	}else{
		end_event=Bcall_TcBk[0].n_events-1;
	}

	const int number_of_events=end_event-start_event+1;

	//if (start_event>7359)	cout<<"start_event="<<start_event<<"    end_event= "<<end_event<<" number_of_events="<<number_of_events<<endl;


	//struct bcall_tcbk * fix=&Bcall_TcBk[0];

	//event
	int32_t event_host[4* number_of_events*n_tcbk];


	k=0;
for(int Bcall_iter=0; Bcall_iter<n_tcbk;Bcall_iter++){
	for(j=start_event; j<=end_event; j++){

		e = &Bcall_TcBk[Bcall_iter].event[j];
		event_host[k]=e->stdv; k++;
		event_host[k]=e->stdv_inv; k++;
		event_host[k]=e->mean; k++;
		event_host[k]=e->log_stdv_1p5; k++;
	}
}


	/* adjust the maximum probability for this round */

for(int Bcall_iter=0; Bcall_iter<n_tcbk;Bcall_iter++){
	//h_amax[Bcall_iter] = INT32_MIN; h_jmax[Bcall_iter]=N_STATES;
	h_temp_sum_pmax[Bcall_iter]=0;

}
	cudaMemset(d_mutex, 0, sizeof(int32_t)*n_tcbk);
	//cudaMemset(d_temp_sum_pmax, 0, sizeof(int64_t)*n_tcbk);



	//copy event_host to constant memory
	cudaMemcpyToSymbol(d_ev, event_host, 4*sizeof(int32_t)*number_of_events*n_tcbk);
	cudaMemcpyToSymbol(d_pmax, p_max, sizeof(int32_t)*n_tcbk);


	//cudaMemcpy(d_amax, h_amax, sizeof(int32_t), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_jmax, h_jmax, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_alpha,alpha, sizeof(int32_t)*N_STATES*2*n_tcbk, cudaMemcpyHostToDevice);
	
	cudaMemcpy(d_temp_sum_pmax, h_temp_sum_pmax, sizeof(int64_t)*n_tcbk, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	kernel<<< number_of_blocks*n_tcbk, threads_per_block >>>(devPtr_pm, devPtr_log_pr, dev_alpha, dev_beta, d_amax, d_jmax, d_mutex, i, d_temp_sum_pmax, start_event,end_event);
	cudaDeviceSynchronize();

	cudaMemcpy(alpha,dev_alpha, sizeof(int32_t)*N_STATES*2*n_tcbk, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_beta,dev_beta, sizeof(uint16_t)*N_STATES*number_of_events*n_tcbk, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_amax,d_amax, sizeof(int32_t)*n_tcbk, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_jmax,d_jmax, sizeof(int32_t)*n_tcbk, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_temp_sum_pmax,d_temp_sum_pmax, sizeof(int64_t)*n_tcbk, cudaMemcpyDeviceToHost);


	for(int Bcall_iter=0; Bcall_iter<n_tcbk; Bcall_iter++){ 
	int counter=0;
	for(k=start_event; k<=end_event; k++){
		for (j = 0; j < N_STATES; ++j) {
			Bcall_TcBk[Bcall_iter].beta[k % TRACEBACK_LEN][j]=h_beta[j+counter*N_STATES  + Bcall_iter*number_of_events*N_STATES];
		}
		counter++;
	}
    }

	if(number_of_events%2!=0){
	    for(int Bcall_iter=0;Bcall_iter<n_tcbk; Bcall_iter++){
		for (j = 0; j < N_STATES; ++j) {
			alpha[j+Bcall_iter*2*N_STATES]=alpha[j+N_STATES +Bcall_iter*2*N_STATES ];
			
		}
	    }
	}


}


void tcbk_sequence_fill(vector<bcall_tcbk>& Bcall_TcBk)
{

	const int n_tcbk= Bcall_TcBk.size();
	cout<<"n_tcbk="<<n_tcbk<<endl;

	/* Pr[ MLSS producing e_1 ... e_i, with S_i == j ] */
	int32_t alpha[N_STATES*2*n_tcbk];
	struct fix_event * e;
	unsigned int j_max[n_tcbk];
	int32_t a_max[n_tcbk];
	int32_t p_max[n_tcbk];
	unsigned int i;
	unsigned int j, k;
	/* Store the maximum probability in each round */
	
	
	//struct bcall_tcbk * fix;

	for(int Bcall_iter=0; Bcall_iter<n_tcbk; Bcall_iter++){

		p_max[Bcall_iter] = INT32_MIN;
		a_max[Bcall_iter] = INT32_MIN;
		j_max[Bcall_iter] = N_STATES;

		//fix= &Bcall_TcBk[Bcall_iter];
	
	
		//fix->sum_log_p_max = 0;
		Bcall_TcBk[Bcall_iter].sum_log_p_max = 0;

	
		//e = &fix->event[0];
		e = &Bcall_TcBk[Bcall_iter].event[0];

		for (j = 0; j < N_STATES; ++j) {
			//struct fix_pm_state * pm = &fix->map[j].pm;
			struct fix_pm_state * pm = &Bcall_TcBk[Bcall_iter].map[j].pm;
			int32_t ln_pe;

			ln_pe = ln_pemission(pm, e, p_max[Bcall_iter],  0);

			if (ln_pe > a_max[Bcall_iter]) {
				a_max[Bcall_iter] = ln_pe;
				j_max[Bcall_iter] = j;
			}

			alpha[j+Bcall_iter*2*N_STATES] = ln_pe;
			//fix->beta[0][j] = N_STATES;
			Bcall_TcBk[Bcall_iter].beta[0][j] = N_STATES;
		}


		/* Save p_max for next round */
		p_max[Bcall_iter] = a_max[Bcall_iter];
		/* Accumulate max probability */
		//fix->sum_log_p_max += a_max[Bcall_iter];
		Bcall_TcBk[Bcall_iter].sum_log_p_max += a_max[Bcall_iter];
	
		//cout<<"after event no 1 and bcall_iter= "<<Bcall_iter<<"   sum_log_p_max = "<<Bcall_TcBk[Bcall_iter].sum_log_p_max<<endl;
	
	}
	//--------------------------------------------------------------------------
	//allocate event independent variables
	//--------------------------------------------------------------------------
	//pm
	int32_t pm_rearranged[7*N_STATES*n_tcbk]; size_t pm_rearranged_size= sizeof(int32_t) * 7*N_STATES*n_tcbk;
	for(int Bcall_iter=0;Bcall_iter<n_tcbk;Bcall_iter++){
		for (j = 0; j < N_STATES; ++j) {
			//struct fix_pm_state * pm = &fix->map[j].pm;
			struct fix_pm_state * pm = &Bcall_TcBk[Bcall_iter].map[j].pm;


			unsigned group_no=j/warp_size;
			unsigned start_index=group_no*warp_size*7 +j%warp_size + Bcall_iter*7*N_STATES;

			pm_rearranged[start_index]=pm->sd_mean;
			pm_rearranged[start_index+1*warp_size]=pm->sd_mean_inv;
			pm_rearranged[start_index+2*warp_size]=pm->sd_lambda_p5;
			pm_rearranged[start_index+3*warp_size]=pm->level_mean;
			pm_rearranged[start_index+4*warp_size]=pm->level_stdv_inv_2;
			pm_rearranged[start_index+5*warp_size]=pm->log_sd_lambda_p5;
			pm_rearranged[start_index+6*warp_size]=pm->log_level_stdv_2pi;

			/*if((j<2) || (j==N_STATES-1)){
				cout<<"iter="<<Bcall_iter<< " j= "<<j<<" pm: "<< pm_rearranged[start_index] <<"  "<< pm_rearranged[start_index+1*warp_size]<<"  "<<pm_rearranged[start_index+5*warp_size]<<"  "<<pm_rearranged[start_index+6*warp_size]<<endl;

			}*/

		}
	}

	//use page lock memory for pm
	//int32_t * devPtr_pm=assign_page_locked_memory_int32_t(pm_rearranged, pm_rearranged_size);

	//instead of page lock memory, global memory is used
	int32_t * devPtr_pm;
	cudaMalloc((void**)&devPtr_pm, pm_rearranged_size) ; // device
	cudaMemcpy(devPtr_pm,pm_rearranged, pm_rearranged_size, cudaMemcpyHostToDevice);

	//log_pr
	int32_t log_pr_rearranged[21*N_STATES*n_tcbk]; size_t log_pr_rearranged_size= sizeof(int32_t) *21* N_STATES *n_tcbk;
	for(int Bcall_iter=0;Bcall_iter<n_tcbk;Bcall_iter++){
		
		for (j = 0; j < N_STATES; ++j) {

			unsigned group_no=j/warp_size;
			unsigned start_index=group_no*warp_size*21 +j%warp_size + Bcall_iter*21*N_STATES;

			for(k=0;k<21;k++){
				//log_pr_rearranged[start_index+ k*warp_size]=fix->map[j].log_pr[k];
				log_pr_rearranged[start_index+ k*warp_size]=Bcall_TcBk[Bcall_iter].map[j].log_pr[k];
			}

			/*if((j==32)||(j==67) || (j==N_STATES-1)){
				cout<<"iter="<<Bcall_iter<< " j= "<<j<<" logpr: "<< log_pr_rearranged[start_index] <<"  "<< log_pr_rearranged[start_index+ 6*warp_size]<<"  "<<log_pr_rearranged[start_index+ 10*warp_size]<<"  "<<log_pr_rearranged[start_index+ 20*warp_size]<<endl;

			}*/
		}
	
	}


	//use page-locked memory
	//int32_t * devPtr_log_pr=assign_page_locked_memory_int32_t(log_pr_rearranged, log_pr_rearranged_size);
	//instead of page lock memory, global memory is used
	int32_t * devPtr_log_pr;
	cudaMalloc((void**)&devPtr_log_pr, log_pr_rearranged_size) ; // device
	cudaMemcpy(devPtr_log_pr,log_pr_rearranged, log_pr_rearranged_size, cudaMemcpyHostToDevice);


	//alpha
	int32_t *dev_alpha;
	cudaMalloc((void**)&dev_alpha, sizeof(int32_t)*N_STATES*2*n_tcbk) ; // device
	//keep alpha in global memory


	//beta needs a write only memory : so global memory
	uint16_t * dev_beta;
	cudaMalloc((void**)&dev_beta, sizeof(uint16_t)*N_STATES*number_of_events_per_kernel*n_tcbk) ; // device
	uint16_t * h_beta;
	h_beta=(uint16_t*)malloc(sizeof(uint16_t)*N_STATES*number_of_events_per_kernel*n_tcbk);



	// a_max and j_max;
	int32_t *h_amax;
	h_amax = (int32_t*)malloc(sizeof(int32_t)*n_tcbk); //host
	int32_t *d_amax;
	cudaMalloc((void**)&d_amax, sizeof(int32_t)*n_tcbk);//device

	

	unsigned int *h_jmax;
	h_jmax = (unsigned int*)malloc(sizeof(unsigned int)*n_tcbk); //host
	unsigned int *d_jmax;
	cudaMalloc((void**)&d_jmax, sizeof(unsigned int)*n_tcbk);//device

	int32_t *d_mutex;
	cudaMalloc((void**)&d_mutex, sizeof(int32_t)*n_tcbk);

	int64_t *h_temp_sum_pmax;
	h_temp_sum_pmax = (int64_t*)malloc(sizeof(int64_t)*n_tcbk); //host
	int64_t *d_temp_sum_pmax;
	cudaMalloc((void**)&d_temp_sum_pmax, sizeof(int64_t)*n_tcbk);//device

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//fix= &Bcall_TcBk[0];
	//for (i = 1; i < fix->n_events; i=i+number_of_events_per_kernel) {
	for (i = 1; i < Bcall_TcBk[0].n_events; i=i+number_of_events_per_kernel) {

		for(int Bcall_iter=0;Bcall_iter<n_tcbk;Bcall_iter++){		

			if ((i % TRACEBACK_CNT) == 1){
				//fixpt_traceback(fix, i, j_max[Bcall_iter]); // need to update beta
				fixpt_traceback(Bcall_TcBk, Bcall_iter, i, j_max[Bcall_iter]); // need to update beta
			}

		}

		//preprocess_for_parallel(fix, alpha, p_max, i, devPtr_pm, devPtr_log_pr, dev_alpha, dev_beta, h_beta, h_amax, d_amax, h_jmax, d_jmax, d_mutex, h_temp_sum_pmax, d_temp_sum_pmax, n_tcbk );
		preprocess_for_parallel(Bcall_TcBk, alpha, p_max, i, devPtr_pm, devPtr_log_pr, dev_alpha, dev_beta, h_beta, h_amax, d_amax, h_jmax, d_jmax, d_mutex, h_temp_sum_pmax, d_temp_sum_pmax);


		
		for(int Bcall_iter=0;Bcall_iter<n_tcbk;Bcall_iter++){
			a_max[Bcall_iter]=h_amax[Bcall_iter]; j_max[Bcall_iter]=h_jmax[Bcall_iter];



			/* Save p_max for next round */
			p_max[Bcall_iter] = a_max[Bcall_iter];


			//if (i<100) cout<<"i= "<<i<< "  Bcall_iter= "<<Bcall_iter<<"   a_max= "<<a_max[Bcall_iter]<< " *h_temp_sum_pmax= "<<h_temp_sum_pmax[Bcall_iter]<<endl;
		
			Bcall_TcBk[Bcall_iter].sum_log_p_max += h_temp_sum_pmax[Bcall_iter]; //a_max;
		}
	

	}
	for(int Bcall_iter=0;Bcall_iter<n_tcbk;Bcall_iter++){
		fixpt_traceback(Bcall_TcBk, Bcall_iter, i, j_max[Bcall_iter]); // need to update beta
	}
	//cout<<"beta[0][0]= "<<fix->beta[0][0]<<"  beta[0][1]= "<<  fix->beta[0][1]<<"  beta[0][200]= "<<  fix->beta[0][200]<<endl;


	//cout<<"state seq [0]= "<<fix->state_seq[0]<<"  state seq [1]= "<<  fix->state_seq[1]<<"  state seq [200]= "<<  fix->state_seq[200]<<endl;

	
	//delete (pm_rearranged); //delete [] devPtr_pm;
	//free(log_pr_rearranged); free(devPtr_log_pr);

	cudaFree(dev_alpha); cudaFree(dev_beta); cudaFree(d_amax); 
	cudaFree(d_jmax); cudaFree(d_mutex); cudaFree (d_temp_sum_pmax);
}



#endif /* VITERBI_TCBK_V1_MULTI_EVENTS_CUH_ */
