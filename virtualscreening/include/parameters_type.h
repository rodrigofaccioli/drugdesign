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
}input_parameters_t;

#endif