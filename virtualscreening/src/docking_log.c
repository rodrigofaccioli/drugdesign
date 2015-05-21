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
#define LOG_FILE_NAME_ADAPTIVE "vs_log_all_docking.log_docking"



void set_log_file_name_adaptive(char *file_name){
	strcpy(file_name, LOG_FILE_NAME_ADAPTIVE);
}

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

void save_log_by_line(const char *path_file_name, const char *line){
	FILE *f_log;
	f_log = open_file(path_file_name, fAPPEND);
	fprintf(f_log,"%s", line);
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
	
	free(str_f_time);
	free(str_s_time);
	fclose(f_log);

}

/*
* This function assigin line which will be stored
* in log file.
* It is used in docking adaptive
*/
void set_log_file_line(char *line, const docking_t *docking, 
	const time_t *f_time, const time_t *s_time){
	
	double diff_t;	
	char *s_diff_t = NULL;
	char *aux = NULL;
	char *str_f_time, *str_s_time;

	aux = (char*)malloc(sizeof(char)*20);
	
	str_f_time = (char*)malloc(sizeof(char)*500);
	str_s_time = (char*)malloc(sizeof(char)*500);

	diff_t = difftime(*f_time, *s_time);
	s_diff_t = (char*)malloc(sizeof(char)*20);
	sprintf (s_diff_t, "%f", diff_t);
	
	strcpy(str_s_time, asctime(gmtime(s_time)) );
	strcpy(str_f_time, asctime(gmtime(f_time)) );

	remove_character_enter(str_s_time);
	remove_character_enter(str_f_time);
	//Assigin line
	strcpy(line, docking->receptor);
	strcat(line,"\t");
	strcat(line, docking->compound);
	strcat(line,"\t");
	strcat(line, s_diff_t);
	strcat(line,"\t");
	int2str(aux, &docking->num_torsion_angle);
	strcat(line, aux);
	strcat(line,"\t");
	int2str(aux, &docking->num_atom);
	strcat(line, aux);
	strcat(line,"\t");
	strcat(line, str_s_time);
	strcat(line,"\t");
	strcat(line, str_f_time);
	strcat(line,"\n");

	free(aux);
	free(s_diff_t);
	free(str_f_time);
	free(str_s_time);	
}
