/*
 ============================================================================
 Name        : thread_read.cu
 Author      : Roksana Hossain
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <thread>
#include <dirent.h>
#include <vector>
#include <string>
#include<mutex>
//#include <fstream>

#include "thread_reading_binary_file.cuh"
#include "bcall.h"
#include "Pore_Model.cuh"
#include "Sequence.cuh"
//#include "Viterbi_TcBk_v1.cuh"
#include "Viterbi_TcBk_v1_1_block_multi_events.cuh"
#include "Viterbi_TcBk_v1_functions.cuh"

#include "reading_tcbk_files.cuh"


#ifdef __unix__
# include <unistd.h>
#elif defined _WIN32
# include <windows.h>
#define sleep(x) Sleep(1000 * x)
#endif

using namespace std;


//--------------------------------------------------------------
//--------------------------------------------------------------
int read_the_bin_file(string file_name, string file_location, string out_file_destination, bool trans, bool fixpt, bool tcbk){

	string temp_dest=file_location;

	file_location.append("/"); file_location.append(file_name);
	//cout<<file_location<<endl;
	//cout<<"file name is "<<file_name<<endl;
	const char * fname= file_location.c_str();



	//char prefix[PATH_MAX + 1];
	//char outfn[PATH_MAX + 1];
	struct sequence * seq;

	uint32_t begin;
	int32_t dt;
	char * cp;
	FILE * f;


	//printf("- Loading file: \"%s\"\n", fname);
	if ((seq = sequence_from_bin_file(fname)) == NULL) {
		fprintf(stderr,
			"%s: invalid sequence file \"%s\"\n", prog, fname);
		fflush(stderr);
		return -1;
	}
	//cout<<"loading completed..."<<endl;

	//------------------------------------------------------------
	//output-file

	std::mutex _mu;

	string out_file_name=file_name.replace(file_name.length()-3,file_name.length()-1,"txt");
	out_file_name=out_file_name.insert(file_name.length()-4, "-bcall_tcbk");
	//cout<<"output file name"<< out_file_name <<endl<<endl;



	out_file_destination= out_file_destination.append("/");
	out_file_destination= out_file_destination.append(out_file_name);

	temp_dest.append("/");temp_dest.append(out_file_name);

	const char * destination= temp_dest.c_str();

	//std::unique_lock<mutex> locker(_mu);
	if (!(f = fopen(destination, "w"))) {
		fprintf(stderr,
				"%s: can't create file \"%s\"\n", prog, fname);
		fflush(stderr);
		return -1;
	};


	//--------------------------------------------------------------
	fprintf(f, "%d", seq->st_cnt);
	for(int i=0; i<seq->st_cnt; i++){
		fprintf(f, " %d", seq->st[i].pm_cnt);

	}
	fprintf(f, "\n");

	for(int i=0; i<seq->st_cnt; i++){

		for(int j=0;j<seq->st[i].pm_cnt;j++){

			//cout<<"file_name="<<file_name<<endl;
			//f<<seq->st[i].id << endl<<file_name<<endl<<seq->st[i].no<<endl<<j<<endl;
			//fprintf(f,"%s\n", file_name);
			fprintf(f, "%s\n%d\n%d\n", seq->st[i].id,seq->st[i].no, j);


			struct bcall_tcbk * bcall_tcbk ;

			bcall_tcbk = tcbk_blk_alloc();


			tcbk_strand_load(bcall_tcbk, seq, i, j);

			if (tcbk) {
				if (trans) {
					//printf("- Computing transitions map...\n");
					fflush(stdout);
					tcbk_compute_transitions(bcall_tcbk, seq->st[i].pm[j].st_param);
				}
			}

			fprintf(f,"log_pr\n");
			for(int k=0;k<N_STATES;k++){
				for(int m=0; m<MAX_TRANS; m++){
					fprintf(f,"%d ", bcall_tcbk->map[k].log_pr[m]);
				}
				fprintf(f,"\n");
			}

			fprintf(f,"pm\n");
			for(int k=0;k<N_STATES;k++){
				fprintf(f,"%d", bcall_tcbk->map[k].pm.level_mean);
				fprintf(f," %d", bcall_tcbk->map[k].pm.sd_mean);
				fprintf(f," %d", bcall_tcbk->map[k].pm.sd_mean_inv);
				fprintf(f," %d", bcall_tcbk->map[k].pm.level_stdv_inv_2);
				fprintf(f," %d", bcall_tcbk->map[k].pm.sd_lambda_p5);
				fprintf(f," %d", bcall_tcbk->map[k].pm.log_level_stdv_2pi);
				fprintf(f," %d\n", bcall_tcbk->map[k].pm.log_sd_lambda_p5);
			}


			fprintf(f,"events\n");
			fprintf(f,"%d\n",bcall_tcbk->n_events);
			for(int k=0;k<bcall_tcbk->n_events;k++){
				fprintf(f,"%d", bcall_tcbk->event[k].mean);
				fprintf(f," %d", bcall_tcbk->event[k].stdv);
				fprintf(f," %d", bcall_tcbk->event[k].stdv_inv);
				fprintf(f," %d\n", bcall_tcbk->event[k].log_stdv_1p5);

			}

		}
	}

	fclose(f);
	//locker.unlock();
	int result=rename(temp_dest.c_str(), out_file_destination.c_str());

//--------------------------------------------------------------



	return 0;
}

//--------------------------------------------------------------

void get_and_read_the_bin_files(){


	char trans[256];
	bool trans_set = false;

	char pmodel[256];
	bool pmodel_set = false;

	bool fixpt = true;
	bool tcbk = true;

	string old_location="/home/roksana/Desktop/bcall_file_steps/only_bin_files";
	string new_location="/home/roksana/Desktop/bcall_file_steps/read_completed";
	string out_file_destination="/home/roksana/Desktop/bcall_file_steps/bcall_tcbk";


	/*string old_location="/home/roksana/Dropbox/bcall_file_steps/only_bin_files";
	string new_location="/home/roksana/Dropbox/bcall_file_steps/read_completed";
	string out_file_destination="/home/roksana/Dropbox/bcall_file_steps/bcall_tcbk";*/

	/*string old_location="/home/shimu/Dropbox/bcall_file_steps/only_bin_files";
	string new_location="/home/shimu/Dropbox/bcall_file_steps/read_completed";
	string out_file_destination="/home/shimu/Dropbox/bcall_file_steps/bcall_tcbk";*/
	
	
	for(int inf_loop=0; inf_loop<100;inf_loop++){

		vector<string> file_names;

		read_the_file_names(file_names, old_location, ".bin");

		//
		print_file_names(file_names);



		// read the files and after the read move the files to a different folder


		for(int i=0;i<file_names.size();i++){
			//read the file
			int r=read_the_bin_file(file_names[i], old_location,  out_file_destination, !trans_set, fixpt, tcbk);
			//remove the file
			remove_the_file(old_location, new_location, file_names, i);
		}
	}

}

//-----------------------------------------------------------------------




void get_bcall_tcbk_files_and_do_viterbi(){

	FILE * f;

	string file_location="/home/roksana/Desktop/bcall_file_steps/bcall_tcbk";
	string destination_file_location="/home/roksana/Desktop/bcall_file_steps/final_result";
	string read_complete_location="/home/roksana/Desktop/bcall_file_steps/read_completed";
	
	/*string file_location="/home/roksana/Dropbox/bcall_file_steps/bcall_tcbk";
	string destination_file_location="/home/roksana/Dropbox/bcall_file_steps/final_result";
	string read_complete_location="/home/roksana/Dropbox/bcall_file_steps/read_completed";*/

	/*string file_location="/home/shimu/Dropbox/bcall_file_steps/bcall_tcbk";
	string destination_file_location="/home/shimu/Dropbox/bcall_file_steps/final_result";
	string read_complete_location="/home/shimu/Dropbox/bcall_file_steps/read_completed";*/
	
	
	int64_t	total_number_of_output_base=0;
	
for(int inf_loop=0; inf_loop<100;inf_loop++){
	vector<string> file_names;
	read_the_file_names(file_names, file_location, "-bcall_tcbk.txt");

	//cout<<"size of file_names: "<<file_names.size()<<endl;

	for(int i=0;i<file_names.size();i++){
		cout<<"file name is "<<file_names[i]<<endl;
		
		//destination file
		string seq_name=file_names[i];seq_name.erase(seq_name.length()-15, 15);
		string destination_file_name= final_destination_file_name(file_names[i], destination_file_location);
		//cout<<"i="<<i<<"    "<<destination_file_name<<endl;

		const char * destination= destination_file_name.c_str();
		if( f = fopen(destination, "w") ){

			vector<bcall_tcbk> Bcall_TcBk;
			vector<int> strand_pm_numbers;

			//cout<<"size of bcall before ="<<Bcall_TcBk.size();
			string strand_id=read_the_tcbk_file(Bcall_TcBk, file_location, file_names[i],strand_pm_numbers);
			//cout<<"strand_id= "<<strand_id<<endl;
			//cout<<"    and after ="<<Bcall_TcBk.size()<<endl;

			/*for(int j=0; j<Bcall_TcBk.size(); j++){

				cout<<" j= "<< j<< " logpr= "<<Bcall_TcBk[j].map[0].log_pr[0]<<"  "<<Bcall_TcBk[j].map[0].log_pr[20]<<endl;

			}*/

			//do viterbi
			do_viterbi(Bcall_TcBk, strand_pm_numbers);
			
			for(int m=0;m<Bcall_TcBk.size();m++){
				cout<<"m="<<m<<"   "<<Bcall_TcBk[m].sum_log_p_max<<endl;
			}


			int counter=0;
			cout<<"strand_pm_numbers.size()= "<<strand_pm_numbers.size()<<endl;

			for(int j=0;j<strand_pm_numbers.size();j++){//iterates strand no times
				cout<<"strand no="<<j<<endl;
				char * s = NULL;
				int64_t p_max=INT64_MIN; //cout<<"INT32_MIN="<<INT32_MIN<<endl;
				for(int k=0;k<strand_pm_numbers[j];k++){
					cout<<"counter="<<counter<<"   "<<Bcall_TcBk[counter].sum_log_p_max<<endl;
					int64_t p=Bcall_TcBk[counter].sum_log_p_max;
					cout<<"pm no="<<k<<"   counter="<<counter<<" p="<<p;
					if(p>p_max){
						p_max=p;
						struct bcall_tcbk* temp_bcall=&Bcall_TcBk[counter];
						s = tcbk_fill_base_seq(temp_bcall);
						//cout<<s<<endl;
						

					}
					cout<<"  pmax="<<p_max<<endl;
					counter++;
				}
				//cout<<strand_id<<endl;
				fprintf(f, ">%s:%s:%d\n", strand_id.c_str(),seq_name.c_str(), j);
				dump_base_seq(f, s);
				//cout<<"length of s="<<strlen(s)<<endl;
				total_number_of_output_base=total_number_of_output_base+strlen(s);
			}
		}else{
			cout<<"destination file is not opening"<<endl;
		}
		fclose(f);

		remove_the_file(file_location, read_complete_location, file_names, i);
	}
}

	cout<<"total_number_of_output_base ="<<total_number_of_output_base<<endl;
}


//----------------------------------------------------------------------



int main(){

	//const int Number_of_files=1;

	//char fname[Number_of_files][PATH_MAX + 1];

	clock_t time_start = clock();
	std::thread t1(get_and_read_the_bin_files);
	
	sleep(1);
	std::thread t2(get_bcall_tcbk_files_and_do_viterbi);


	t1.join();
	t2.join();

	clock_t time_end = clock();
	printf("Time taken: %.6fs\n", (double)(time_end - time_start)/CLOCKS_PER_SEC);

	cout<<"end of main"<<endl;
	return 0;
}
