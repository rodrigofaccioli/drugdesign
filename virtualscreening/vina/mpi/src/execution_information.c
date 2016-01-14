#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "futil.h"
#include "string_owner.h"
#include "execution_information.h"

#define NAME_FILE "vs_execution_information.txt"

void save_information_adaptive(const char *local_execute, const int *num_proc, 
	const int *num_dock_total, const int *nthreads){
	
	FILE *f_info;
	char *f_information = NULL;
	char *path_file_name = NULL;	
	int i;
	f_information = (char*)malloc(sizeof(char)*MAX_FILE_NAME);

	strcpy(f_information, NAME_FILE);

	path_file_name = path_join_file(local_execute, f_information);

	f_info = open_file(path_file_name, fWRITE);
	fprintf(f_info, "Process Number = %d\n", *num_proc);
	fprintf(f_info, "Threads Number = %d\n", *nthreads);	
	fprintf(f_info, "Docking Number total %d\n", *num_dock_total);

	fclose(f_info);

	free(path_file_name);
	free(f_information);
}

void save_information(const char *local_execute, const int *num_proc, 
	const int *v_num_dock, const int *nthreads){
	
	FILE *f_info;
	char *f_information = NULL;
	char *path_file_name = NULL;
	int num_dock_total;
	int i;
	f_information = (char*)malloc(sizeof(char)*MAX_FILE_NAME);

	strcpy(f_information, NAME_FILE);

	path_file_name = path_join_file(local_execute, f_information);

	f_info = open_file(path_file_name, fWRITE);
	fprintf(f_info, "Process Number = %d\n", *num_proc);
	fprintf(f_info, "Threads Number = %d\n", *nthreads);
	num_dock_total = 0;
	fprintf(f_info, "Docking Number per process\n");
	for (i = 0; i < *num_proc; i++){
		fprintf(f_info, "\t%d\t%d\n", i, v_num_dock[i]);
		num_dock_total = num_dock_total + v_num_dock[i];
	}		
	fprintf(f_info, "Docking Number total %d\n", num_dock_total);

	fclose(f_info);

	free(path_file_name);
	free(f_information);
}

void saving_time_execution(const char *local, const double *finished, const double *started,
	const time_t *f_time, const time_t *s_time){
	FILE *f_info;
	char *f_information = NULL;
	char *path_file_name = NULL;
	double diff_t;	
	f_information = (char*)malloc(sizeof(char)*MAX_FILE_NAME);	

	diff_t = difftime(*f_time, *s_time);
	strcpy(f_information, NAME_FILE);
	path_file_name = path_join_file(local, f_information);
	f_info = open_file(path_file_name, fAPPEND);	
	fprintf(f_info, "\n\nDate \n");
	fprintf(f_info,"Started on %s", asctime(gmtime(s_time)));
	fprintf(f_info, "Finished on %s", asctime(gmtime(f_time)));	
	fprintf(f_info, "Total time %f seconds \n",diff_t);		

	fprintf(f_info, "\n\nTime Execution\n");
	fprintf(f_info, "Started = %f\n", *started);
	fprintf(f_info, "Finished = %f\n", *finished);
	fprintf(f_info, "Runtime = %9.4f\n", *finished - *started);
	

	fclose(f_info);
	free(path_file_name);
	free(f_information);
}