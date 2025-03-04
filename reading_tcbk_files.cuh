/*
 * reading_tcbk_files.cuh
 *
 *  Created on: Jun 15, 2017
 *      Author: shimu
 */

#ifndef READING_TCBK_FILES_CUH_
#define READING_TCBK_FILES_CUH_

#include <iostream>
#include <numeric>
#include <stdlib.h>
#include<vector>
#include <string>
#include <fstream>
//#include"Viterbi_TcBk_v1.cuh"
#include "Viterbi_TcBk_v1_1_block_multi_events.cuh"

using namespace std;


string final_destination_file_name(string a, string destination_file_location){

	string destination_file_name=a;
	destination_file_name.replace(destination_file_name.length()-14,10,"final");
	destination_file_name.insert(0,"/" );
	destination_file_name.insert(0, destination_file_location);
	return destination_file_name;
}

string read_the_tcbk_file(vector<bcall_tcbk>& Bcall_TcBk, string location, string file_name, vector <int>& strand_pm){

	//cout<<file_name<<endl;

	location.append("/");
	location=location.append(file_name);

	ifstream myfile(location);

	string strand_id="";

	if(myfile.is_open()){
		int number_of_strand;

		myfile>>number_of_strand;
		//cout<<"number_of_strand"<<number_of_strand;
		//vector<int> strand_pm;
		strand_pm.resize(number_of_strand);
		int count=0;
		for(int i=0;i<number_of_strand;i++){
			myfile>>strand_pm[i];
			//cout<<" "<<strand_pm[i];
			count=count+strand_pm[i];
		}
		//cout<<"  count="<<count<<endl;
		Bcall_TcBk.resize(count);



		int Bcall_TcBk_index=0;
		for(int i=0;i<number_of_strand;i++){
			for(int j=0; j<strand_pm[i]; j++){
				string a;
				if(strand_id==""){
					myfile>>strand_id;
					//cout<<"strand_id="<<strand_id<<endl;
				}else{
					//
					myfile>>a;
					assert(a==strand_id);
				}
				//
				int strand_no; myfile>>strand_no;
				assert(strand_no==i);

				int pm_no; myfile>>pm_no;
				assert(pm_no==j);

				 myfile>>a;
				assert(a=="log_pr");
				for(int n=0;n<N_STATES;n++){
					for(int m=0;m<MAX_TRANS;m++){
						myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].log_pr[m];
					}
				}

				myfile>>a;
				assert(a=="pm");
				for(int n=0;n<N_STATES;n++){

					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.level_mean;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.sd_mean;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.sd_mean_inv;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.level_stdv_inv_2;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.sd_lambda_p5;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.log_level_stdv_2pi;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].map[n].pm.log_sd_lambda_p5;

				}

				myfile>>a;
				assert(a=="events");
				myfile>>Bcall_TcBk[Bcall_TcBk_index].n_events;
				for(int n=0;n<Bcall_TcBk[Bcall_TcBk_index].n_events;n++){
					myfile>> Bcall_TcBk[Bcall_TcBk_index].event[n].mean;
					myfile>>Bcall_TcBk[Bcall_TcBk_index].event[n].stdv;
					myfile>> Bcall_TcBk[Bcall_TcBk_index].event[n].stdv_inv;
					myfile>> Bcall_TcBk[Bcall_TcBk_index].event[n].log_stdv_1p5;

				}

				Bcall_TcBk_index++;

			}
		}

	}else{
		cout<<"unabled to open the file"<<endl;
	}
	myfile.close();
	return strand_id;
}

void do_viterbi(vector<bcall_tcbk>& Bcall_TcBk, vector<int> strand_pm_numbers){

	/*	
	cout<<"strand_pm_numbers.size="<<strand_pm_numbers.size()<<endl;

	for(int i=0;i<strand_pm_numbers.size(); i++){
		cout<<"i = "<< i<<"   strand_pm_numbers= "<<strand_pm_numbers[i]<<endl; 

	}
	cout<<"inside to viterbi. Bcall_TcBk.size()="<<Bcall_TcBk.size()<<endl;
*/
	int counter1=0; int counter2=0;
	for(int i=0;i<strand_pm_numbers.size(); i++){
		vector<bcall_tcbk> temp_Bcall_TcBk; temp_Bcall_TcBk.resize(strand_pm_numbers[i]); 
		for(int j=0;j<strand_pm_numbers[i];j++){
			temp_Bcall_TcBk[j]=Bcall_TcBk[counter1]; 
			counter1++;
		}

		tcbk_sequence_fill(temp_Bcall_TcBk);
		for(int j=0;j<strand_pm_numbers[i];j++){
			Bcall_TcBk[counter2]=temp_Bcall_TcBk[j]; 
			counter2++;
		}

	}


/*	for(int i=0;i<Bcall_TcBk.size(); i++){
		vector<bcall_tcbk> temp_Bcall_TcBk; temp_Bcall_TcBk.resize(1);

		temp_Bcall_TcBk[0]=Bcall_TcBk[i]; 
		
		tcbk_sequence_fill(temp_Bcall_TcBk);
		Bcall_TcBk[i]=temp_Bcall_TcBk[0]; 
	}
*/

	for(int i=0;i<Bcall_TcBk.size();i++){
		//struct bcall_tcbk* bcall=&Bcall_TcBk[i];

		//tcbk_sequence_fill(Bcall_TcBk, i);
		//tcbk_sequence_fill(bcall);

		struct bcall_tcbk* bcall=&Bcall_TcBk[i];
		int p = tcbk_fill_state_seq(bcall);
	}
	for(int i=0;i<Bcall_TcBk.size();i++){
		cout<<Bcall_TcBk[i].sum_log_p_max<<endl;
	}


}


#endif /* READING_TCBK_FILES_CUH_ */
