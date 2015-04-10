#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "futil.h"
#include "string_owner.h"
#include "execution_information.h"

#define NAME_FILE "vs_execution_information.txt"

void save_information(const char *local_execute, const int *num_proc, 
	const int *num_dock, const int *num_dock_root, const int *nthreads){
	
	FILE *f_info;
	char *f_information = NULL;
	char *path_file_name = NULL;
	f_information = (char*)malloc(sizeof(char)*MAX_FILE_NAME);

	strcpy(f_information, NAME_FILE);

	path_file_name = path_join_file(local_execute, f_information);

	f_info = open_file(path_file_name, fWRITE);
	fprintf(f_info, "Process Number = %d\n", *num_proc);
	fprintf(f_info, "Threads Number = %d\n", *nthreads);
	fprintf(f_info, "Docking Number root = %d\n", *num_dock_root);
	fprintf(f_info, "Docking Number others = %d\n", *num_dock);
	fprintf(f_info, "Docking Number total %d\n", (*num_dock * (*num_proc-1) ) + *num_dock_root);

	fclose(f_info);

	free(path_file_name);
	free(f_information);
}

void saving_time_execution(const char *local, const double *finished, const double *started,
	const time_t *f_time, const time_t *s_time){
	FILE *f_info;
	char *f_information = NULL;
	char *path_file_name = NULL;
	double elapsed;	
	f_information = (char*)malloc(sizeof(char)*MAX_FILE_NAME);	

	strcpy(f_information, NAME_FILE);
	path_file_name = path_join_file(local, f_information);
	f_info = open_file(path_file_name, fAPPEND);	
	fprintf(f_info, "\n\nDate \n");
	fprintf(f_info,"Started on %s \n",asctime(gmtime(s_time)));
	fprintf(f_info, "Finished on %s \n",asctime(gmtime(f_time)));	
	fprintf(f_info, "\n\nTime Execution\n");
	fprintf(f_info, "Started = %f\n", *started);
	fprintf(f_info, "Finished = %f\n", *finished);
	fprintf(f_info, "Runtime = %9.4f\n", *finished - *started);
	

	fclose(f_info);
	free(path_file_name);
	free(f_information);
}