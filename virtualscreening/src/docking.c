#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "docking.h"
#include "defines.h"
#include "futil.h"
#include "string_owner.h"


void deAllocate_docking(docking_t *d){
	free(d);
}

int get_number_docking(const char *file_name){

	FILE *f_dock=NULL;
	char *line=NULL;
	int num_dock = -1;

	f_dock = open_file(file_name, fREAD);
	line = (char*)malloc(MAX_LINE_FILE);
	fgets(line, MAX_LINE_FILE, f_dock);
	num_dock = str2int(line);
  	free(line);
  	fclose(f_dock);	
  	return num_dock;
}

void set_receptor_compound(char *receptor, char *compound, char *source){
	char *ch=NULL;
    ch = strtok(source, "\t");
    ltrim(ch);
    strcpy(receptor, ch);
    ch = strtok(NULL, "\t");
    ltrim(ch);
    strcpy(compound, ch);
}