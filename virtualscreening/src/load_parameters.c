#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include "defines.h"
#include "load_parameters.h"
#include "futil.h"
#include "string_owner.h"
#include "messages.h"

void show_parameters(const input_parameters_t *param){
	printf("Local_Execute %s\n", param->local_execute);	
	printf("Path_Receptors %s\n", param->path_receptors);	
	printf("Path_compounds %s\n", param->path_compounds);
	printf("Path_out %s\n", param->path_out);
	printf("Path_log %s\n", param->path_log);
	printf("Config_file %s\n", param->config_vina);
	printf("Vina_program %s\n", param->vina_program);
}

static void set_parameters(input_parameters_t *param, const char *key, 
	const char *value){
	if (strncmp(key,"Local_Execute", 13) == 0){
		strcpy(param->local_execute, value);
	}else if (strncmp(key,"Path_receptor", 13) == 0){
		strcpy(param->path_receptors, value);
	}else if (strncmp(key,"Path_compounds",14) == 0){
		strcpy(param->path_compounds, value);
	}else if (strncmp(key,"Path_out",8) == 0){
		strcpy(param->path_out, value);
	}else if (strncmp(key,"Path_log",8) == 0){
		strcpy(param->path_log, value);
	}else if (strncmp(key,"Config_file",11) == 0){
		strcpy(param->config_vina, value);
	}else if (strncmp(key,"Vina_program",12) == 0){
		strcpy(param->vina_program, value);
	}else{
		fatal_error("Parameter not found\n");
	}
}

void load_parameters_from_file(input_parameters_t *in_param,
		const char *conf_file_name){
	FILE *conf=NULL;
	char *ch=NULL;
	char *line=NULL;
	char *key=NULL;
	char *value=NULL;
	int key_r;

	conf = open_file(conf_file_name, fREAD);
	line = (char*)malloc(MAX_LINE_FILE);
	key = (char*)malloc(MAX_LINE_FILE);
	value = (char*)malloc(MAX_LINE_FILE);
	key_r = 0;
	while (fgets(line, MAX_LINE_FILE, conf) != NULL){		
		ch = strtok(line, "=");
		while (ch != NULL) {
			trim(ch);
			remove_character_enter(ch);
  			if (key_r == 0){
  				strcpy(key, ch);
  				key_r = 1;  				
  			}else{
  				strcpy(value, ch);
  			}
  			ch = strtok(NULL, "=");
  		}
  		set_parameters(in_param, key, value);
  		key_r = 0;
	}
	free(key);
	free(value);
  	free(line);
  	fclose(conf);
}

void deAllocateload_parameters(input_parameters_t *param){
	free(param);
}
