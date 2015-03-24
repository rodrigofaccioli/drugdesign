#ifndef OLD_PARAMETERS_TYPE_H
#define OLD_PARAMETERS_TYPE_H

#include "defines.h"

typedef struct s_input_parameters{
	//main local execute
	char local_execute[MAX_PATH];
	//where the receptors are
	char path_receptors[MAX_PATH];
	//where the compounds are
	char path_compounds[MAX_PATH];
	//where the structures are saved
	char path_out[MAX_PATH];
	//where the log are saved
	char path_log[MAX_PATH];	
	//Vina file name config
	char config_vina[MAX_PATH_FILE_NAME];
	//Path File name of Vina
	char vina_program[MAX_PATH_FILE_NAME];
}input_parameters_t;

#endif