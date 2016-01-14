/*
 * Routines to prepare ligands for Virtual Screening
 * These routines were developed by:
 * Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
 * Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com
*/
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "defines.h"
#include "messages.h"
#include "load_parameters.h"
#include "futil.h"
#include "string_owner.h"
#include "vina.h"


int main(int argc, char *argv[]) {
	//Preparing input parameters data structure
	input_parameters_t *param=NULL;
	//Stores mol2 extension
	char *ex_mol2 = NULL;
	//Stores pdbqt extension
	char *ex_pdbqt = NULL;
	//Array to all mol2 files
	owner_file_t *all_mol2_files = NULL;
	//Size of all_mol2_files
	int mol2_size;
	//Used at main for command
	int m;
	//Path file name of mol2
	char *path_file_mol2=NULL;
	//Path file name of pdbqt
	char *path_file_pdbqt=NULL;
	//Base name of File. It is file name without extension
	char *base_file_name=NULL;
	//File name used to prepare ligand
	char *file_pdbqt=NULL;

	param = (input_parameters_t*)malloc(sizeof(input_parameters_t));    
	//loading input parameters files  
	load_parameters_from_file(param,argv[1]);

	//Setting pdbqt extension
	ex_pdbqt = (char*)malloc(sizeof(char)*MAX_SUFIX);
	strcpy(ex_pdbqt, "pdbqt");

	//Loading all Mol2 files
	ex_mol2 = (char*)malloc(sizeof(char)*MAX_SUFIX);
	strcpy(ex_mol2, "mol2");
	mol2_size = how_many_files_directory_by_extension(param->path_mol2, ex_mol2);		
	all_mol2_files = allocate_owner_file(&mol2_size);
	insert_files_directory_by_extension(all_mol2_files, param->path_mol2, ex_mol2);

	initialize_vina_execution();
	//Calling prepare ligand 
    #pragma omp parallel shared(m,mol2_size) private(path_file_mol2, path_file_pdbqt, file_pdbqt, base_file_name)
    {           	
		path_file_mol2 = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
		path_file_pdbqt = (char*)malloc(sizeof(char)*MAX_PATH_FILE_NAME);
		file_pdbqt = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
		base_file_name  = (char*)malloc(sizeof(char)*MAX_FILE_NAME);     		    	
		#pragma omp for 		
		for (m = 0; m < mol2_size; m++){
			//Obtaing base file name. This is file name without its extension
			set_base_file_name(base_file_name, all_mol2_files[m].file_name, ex_mol2);
			strcpy(file_pdbqt, base_file_name);
			add_ext_file_name(file_pdbqt, ex_pdbqt);
			//mol2 file
			strcpy(path_file_mol2, param->path_mol2);
			strcat(path_file_mol2, all_mol2_files[m].file_name);
			//pdbqt files		
			strcpy(path_file_pdbqt, param->path_compounds);
			strcat(path_file_pdbqt, file_pdbqt);
			call_prepare_ligand(param, path_file_mol2, path_file_pdbqt);

		}
		free(base_file_name);
		free(file_pdbqt);
		free(path_file_pdbqt);		
		free(path_file_mol2);		    		
    }
	finish_vina_execution();

	free(ex_pdbqt);
	desallocate_owner_file(all_mol2_files, &mol2_size);
	free(ex_mol2);
	deAllocateload_parameters(param);
	return 0;
}