#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "vina.h"
#include "messages.h"

static char *command = NULL;
static char *program = NULL; /* program (executable) file name */
static char *opt_config = NULL;
static char *opt_path_file_config = NULL;
static char *opt_receptor = NULL;
static char *opt_path_file_receptor = NULL;
static char *opt_ligand = NULL;
static char *opt_path_file_ligand = NULL;
static char *opt_out = NULL;
static char *opt_path_file_out = NULL;
static char *opt_log = NULL;
static char *opt_path_file_log = NULL;
static char *file_log_name = NULL;
static char *file_receptor_name = NULL;
static char *file_ligand_name = NULL;
static char *file_out_name = NULL;
static char *opt_l = NULL;
static char *opt_v = NULL;
static char *opt_u = NULL;
static char *opt_a = NULL;
static char *opt_o = NULL;
static char *opt_nphs_lps = NULL;
static char *opt_hydrogens = NULL;

static inline int run_vina_program(const char *file, char *const argv[]){
  int out, err; /* file descriptors for stdout and stderr */
  int i = 0;
  strcpy(command, argv[i]);
  strcat(command, " ");
  i++;
  do{ 
    strcat(command, argv[i]);
    strcat(command, " ");
    i++;  
  }while (argv[i] != NULL);
  
  //Avoid output messages
  strcat(command, " > /dev/null 2> /dev/null ");

  system(command);

  return 1;
}

void initialize_vina_execution(){
	command = (char*)malloc(sizeof(char)*MAX_COMMAND);
	program = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
	opt_config = (char*)malloc(sizeof(char)*10);
	opt_path_file_config = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
	opt_receptor = (char*)malloc(sizeof(char)*12);
	opt_path_file_receptor  = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
	opt_ligand = (char*)malloc(sizeof(char)*10);
	opt_path_file_ligand = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
	opt_out = (char*)malloc(sizeof(char)*7);
	opt_path_file_out = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
	opt_log = (char*)malloc(sizeof(char)*7);  
	opt_path_file_log = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);

  file_log_name = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
  file_receptor_name = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
  file_ligand_name = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
  file_out_name = (char*)malloc(sizeof(char)*MAX_FILE_NAME);

  opt_l = (char*)malloc(sizeof(char)*3);
  opt_v = (char*)malloc(sizeof(char)*3);
  opt_u = (char*)malloc(sizeof(char)*3);
  opt_a = (char*)malloc(sizeof(char)*3);
  opt_o = (char*)malloc(sizeof(char)*3);
  opt_nphs_lps = (char*)malloc(sizeof(char)*10);
  opt_hydrogens = (char*)malloc(sizeof(char)*11);


  //setting option values
  strcpy(opt_config, "--config");
  strcpy(opt_receptor, "--receptor");
  strcpy(opt_ligand, "--ligand");
  strcpy(opt_out, "--out");
  strcpy(opt_log, "--log");
  strcpy(opt_l, "-l");
  strcpy(opt_v, "-v");
  strcpy(opt_u, "-U");
  strcpy(opt_a, "-A");
  strcpy(opt_o, "-o");
  strcpy(opt_nphs_lps, "nphs_lps");
  strcpy(opt_hydrogens, "hydrogens");
}

void finish_vina_execution(){
	free(command);
	free(program);
	free(opt_config);
	free(opt_path_file_config);
	free(opt_receptor);
	free(opt_path_file_receptor);
	free(opt_ligand);
	free(opt_path_file_ligand);
	free(opt_out);
	free(opt_path_file_out);
	free(opt_log);
	free(opt_path_file_log);
  free(file_log_name);
  free(file_receptor_name);
  free(file_ligand_name);
  free(file_out_name);
  free(opt_l);
  free(opt_v);
  free(opt_u);
  free(opt_a);
  free(opt_o);
  free(opt_nphs_lps);
  free(opt_hydrogens);

}

/* Setting log file name
*/
void set_log_file(char *file_name, const char *receptor, const char *ligand){
  strcpy(file_name, receptor);
  strcat(file_name, "_-_");
  strcat(file_name, ligand);
  strcat(file_name, ".log");
}

/* Setting out file name
*/
void set_out_file(char *file_name, const char *receptor, const char *ligand){
  strcpy(file_name, receptor);
  strcat(file_name, "_-_");
  strcat(file_name, ligand);
  strcat(file_name, ".pdbqt");
}

/* Setting ligand file name
*/
void set_ligand_file(char *file_name, const char *ligand){
  strcpy(file_name, ligand);
  strcat(file_name, ".pdbqt");
}

/* Setting out file name
*/
void set_receptor_file(char *file_name, const char *receptor){
  strcpy(file_name, receptor);
  strcat(file_name, ".pdbqt");  
}


/* Executes virtual screening
*/
void call_vina(const input_parameters_t *param, const docking_t *docking){  
  char *vina_args[11];

  strcpy(program, param->vina_program);
  vina_args[0] = program;
  //config
  vina_args[1] = opt_config;
  strcpy(opt_path_file_config, param->config_vina);  
  vina_args[2] = opt_path_file_config;
  //receptor
  set_receptor_file(file_receptor_name, docking->receptor);
  vina_args[3] = opt_receptor;
  strcpy(opt_path_file_receptor, param->path_receptors); 
  strcat(opt_path_file_receptor, file_receptor_name);
  vina_args[4] = opt_path_file_receptor;
  //ligand
  set_ligand_file(file_ligand_name, docking->compound);
  vina_args[5] = opt_ligand;
  strcpy(opt_path_file_ligand, param->path_compounds); 
  strcat(opt_path_file_ligand, file_ligand_name);
  vina_args[6] = opt_path_file_ligand;
  //out
  set_out_file(file_out_name, docking->receptor, docking->compound);
  vina_args[7] = opt_out;
  strcpy(opt_path_file_out, param->path_out); 
  strcat(opt_path_file_out, file_out_name);
  vina_args[8] = opt_path_file_out;
  //log
  set_log_file(file_log_name, docking->receptor, docking->compound);
  vina_args[9] = opt_log;
  strcpy(opt_path_file_log, param->path_log); 
  strcat(opt_path_file_log, file_log_name);
  vina_args[10] = opt_path_file_log;

  vina_args[11] = NULL;

  if (!run_vina_program(program, vina_args)){
    fatal_error("Failed to run Vina program \n");
  }

}



void call_prepare_ligand(const input_parameters_t *param, const char *fmol2, const char *fpdbqt){
    char *ligand_args[12];
    char *script_ligand4=NULL;
    char *mol2_aux = NULL;
    char *pdbqt_aux = NULL;

    script_ligand4 = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
    mol2_aux = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
    pdbqt_aux = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);

    strcpy(program, param->pythonsh);
    ligand_args[0] = program;
    //script ligand4
    strcpy(script_ligand4, param->script_ligand4);
    ligand_args[1] = script_ligand4;
    ligand_args[2] = opt_l;
    //mol2 file
    strcpy(mol2_aux, fmol2);
    ligand_args[3] = mol2_aux;
    ligand_args[4] = opt_v;
    ligand_args[5] = opt_o;
    //pdbqt file
    strcpy(pdbqt_aux, fpdbqt);
    ligand_args[6] = pdbqt_aux;
    ligand_args[7] = opt_u;
    ligand_args[8] = opt_nphs_lps;
    ligand_args[9] = opt_a;
    ligand_args[10] = opt_hydrogens;

    ligand_args[11] = NULL;

    if (!run_vina_program(program, ligand_args) ){
      fatal_error("Failed to prepare ligand \n");     
    }

    free(mol2_aux);
    free(pdbqt_aux);
    free(script_ligand4);
}