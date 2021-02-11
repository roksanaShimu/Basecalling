/*
 * thread_reading_binary_file.cuh
 *
 *  Created on: Jun 14, 2017
 *      Author: shimu
 */

#ifndef THREAD_READING_BINARY_FILE_CUH_
#define THREAD_READING_BINARY_FILE_CUH_

#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <thread>
#include <dirent.h>
#include <vector>
#include <string>

//#include "Sequence.cuh"
//#include "bcall.h"
//#include "Pore_Model.cuh"

using namespace std;

int verbose = 0;
char * prog;


//get the bin file names from the directory
void read_the_file_names(vector<string>& file_names, string location, string extension){
	DIR *dir;
	struct dirent *ent;

	if ((dir = opendir (location.c_str())) != NULL) {
		/* print all the files and directories within directory */
		 while ((ent = readdir (dir)) != NULL) {
		    //printf ("%s\n", ent->d_name);
		    string s=ent->d_name;
		    int d=extension.length();
		    if(s.length()>d){
		    	 string s2=s.substr(s.length()-d, s.length()-1);
		    	 //cout<<s2<<endl;
		    	 if(s2==extension){
		    		 file_names.push_back(s);
		    	 }
		    }
		 }
		 closedir (dir);
	} else {
		/* could not open directory */
		perror ("");
	}


}

//----------------------------------------------------------------------
void print_file_names(vector<string>& file_names){
	//check the list of file_names
	cout<<endl<<endl;
	for(int i=0;i<file_names.size();i++){
		cout<<file_names[i]<<endl;
	}

}


//----------------------------------------------------------------------
//after reading the binary files, remove the binary file
void remove_the_file(string source, string destination, vector<string>& file_names, int i){
	source.append("/"); source.append(file_names[i]);
	//cout<<"source ="<<source<<endl;

	destination.append("/"); destination.append(file_names[i]);
	//cout<<"destination="<<destination<<endl<<endl;

	int result= rename(source.c_str(), destination.c_str());
	/* if ( result != 0 )
	    perror( "Error renaming file" );*/


}

#endif /* THREAD_READING_BINARY_FILE_CUH_ */
