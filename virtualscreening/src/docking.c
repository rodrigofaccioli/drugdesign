#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "docking.h"
#include "defines.h"
#include "futil.h"
#include "string_owner.h"

#define BASE_NAME_FILE_DOCKING "docking_rank_"

void deAllocate_docking(docking_t *d){
	free(d);
}

/** Calculates the number of docking for each process
* v_dock is an array that contains the number of docking that each process will execute
*/
void set_number_docking(int *v_dock, const int *world_size, const int *full_dock){
  int i, j, s;
  int num_dock;  
  int num_dock_remainder;

  num_dock = (int)*full_dock  / *world_size;
  num_dock_remainder = *full_dock - (num_dock* *world_size);
  //Assign equal number of docking for all process 
  for (i = 0; i < *world_size; i++){
      v_dock[i] = num_dock;
  }
  //Assign docking number which are remainder for each process
  j=0;
  while (j < num_dock_remainder){
    s = 0;
    while ( (j < num_dock_remainder) && (s < *world_size) ){
        v_dock[s] = v_dock[s] + 1;
        j = j + 1;
        s = s + 1;
    }
  }  
}


int get_number_docking_from_file(const char *file_name){

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
    trim(ch);
    remove_character_enter(ch);
    strcpy(compound, ch);
}

void build_docking_file_name(char *f_file, const int *suf){
  char *c_sufix = NULL;
  c_sufix = (char*)malloc(sizeof(char)*MAX_SUFIX);
  strcpy(f_file, BASE_NAME_FILE_DOCKING);
  int2str(c_sufix, suf);
  strcat(f_file, c_sufix);
  strcat(f_file, ".dock");
  free(c_sufix);
}

void get_docking_file_name(char *f_file, const int *suf){  
  build_docking_file_name(f_file, suf);
}

//create a file with docking
void save_file_docking_from_array(const docking_t *v_doc, const int *num_doc, 
  const char *local_execute, const int *suf){
  FILE *f_dock=NULL;
  char *line=NULL;  
  char *f_file = NULL;
  char *path_file_name = NULL;  
  int d;

  f_file = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
  build_docking_file_name(f_file, suf);
  path_file_name = path_join_file(local_execute, f_file);

  f_dock = open_file(path_file_name, fWRITE);
  fprintf(f_dock, "%d\n", *num_doc);
  for(d = 0; d < *num_doc; d++){
    fprintf(f_dock, "%s\t%s\n", v_doc[d].receptor, v_doc[d].compound);
  }
  fclose(f_dock);
  
  free(path_file_name);
  free(f_file);
}

void load_docking_from_file(docking_t *v_doc, const int *num_dock, 
  const char *local, const int *suf){
  char *path_file_name = NULL;  
  FILE *f_dock=NULL;
  char *f_file = NULL;
  char *line=NULL;
  int d;

  line = (char*)malloc(MAX_LINE_FILE);
  
  f_file = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
  build_docking_file_name(f_file, suf);
  path_file_name = path_join_file(local, f_file);

  f_dock = open_file(path_file_name, fREAD);
  //ignoring first line of file
  fgets(line, MAX_LINE_FILE, f_dock);
  for (d = 0; d < *num_dock; d++){
      fgets(line, MAX_LINE_FILE, f_dock);
      set_receptor_compound(v_doc[d].receptor, v_doc[d].compound, line);      
  }
  fclose(f_dock);     
  free(line);
  free(f_file);
  free(path_file_name);
}