#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "futil.h"
#include "string_owner.h"

void save_information(const char *local_execute, const int *num_proc, 
	const int *num_dock, const int *num_dock_root){
	
	FILE *f_info;
	char *f_information = NULL;
	char *path_file_name = NULL;
	f_information = (char*)malloc(sizeof(char)*MAX_FILE_NAME);

	strcpy(f_information, "vs_execution_information.txt");

	path_file_name = path_join_file(local_execute, f_information);

	f_info = open_file(path_file_name, fWRITE);
	fprintf(f_info, "Process Number = %d\n", *num_proc);
	fprintf(f_info, "Docking Number root = %d\n", *num_dock_root);
	fprintf(f_info, "Docking Number others = %d\n", *num_dock);
	fprintf(f_info, "Docking Number total %d\n", (*num_dock * (*num_proc-1) ) + *num_dock_root);

	fclose(f_info);

	free(path_file_name);
	free(f_information);
}