#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "docking.h"
#include "messages.h"
#include "futil.h"
#include "string_owner.h"
#include "docking_log.h"

#define BASE_NAME_LOG_FILE "vs_log_"

void build_log_file_name_from_rank(char *file_name, const int *rank){
  	char *c_sufix = NULL;
  	c_sufix = (char*)malloc(sizeof(char)*MAX_SUFIX);	
	strcpy(file_name, BASE_NAME_LOG_FILE);
  	int2str(c_sufix, rank);
  	strcat(file_name, c_sufix);
  	strcat(file_name, ".log_docking");	
	free(c_sufix);

}

void initialize_log(const char *path_file_name){
	FILE *f_log;	
	f_log = open_file(path_file_name, fWRITE);
	fclose(f_log);
}

void save_log(const char *path_file_name, const docking_t *docking, 
	const time_t *f_time, const time_t *s_time){
	
	FILE *f_log;
	double diff_t;	
	char *str_f_time, *str_s_time;
	
	str_f_time = (char*)malloc(sizeof(char)*500);
	str_s_time = (char*)malloc(sizeof(char)*500);

	diff_t = difftime(*f_time, *s_time);
	f_log = open_file(path_file_name, fAPPEND);

	strcpy(str_s_time, asctime(gmtime(s_time)) );
	strcpy(str_f_time, asctime(gmtime(f_time)) );

	remove_character_enter(str_s_time);
	remove_character_enter(str_f_time);
	fprintf(f_log, "%s\t%s\t%f\t%d\t%d\t%s\t%s\n", docking->receptor,
		docking->compound,
		diff_t,
		docking->num_torsion_angle,
		docking->num_atom,
		str_s_time,
		str_f_time);
	
	fclose(f_log);

}