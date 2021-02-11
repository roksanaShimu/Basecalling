Add libraries:

export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export PATH=/usr/local/cuda-8.0/bin${PATH:+:${PATH}}

To the directory ./Final_Basecalling

Open thread_read.cu

update line# 192 - 194  and 239-241 accordingly

save the file.

From copy all the 8 files from ./Final_Basecalling/Data_set/Data8 and paste in the folder ./Final_Basecalling/bcall_file_steps/only_bin_files

From the terminal compile the .cu file:
$nvcc thread_read.cu -use_fast_math -Xptxas -v -lcurand -std=c++11

then
$./a.out

The final results will be stored in ./Final_Basecalling/bcall_file_steps/final_result 




