/*
 * Routines to performe virtual screening using mpi
 * These routines were developed by:
 * Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
 * Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include "defines.h"
#include "docking.h"
#include "messages.h"
#include "load_parameters.h"
#include "futil.h"
#include "string_owner.h"
#include "mpi_parameters_type.h"
#include "mpi_docking_type.h"
#include "vina.h"

// Calculates the number of docking for each process
void get_number_docking(int *n_dock, int *n_dock_root, 
  const int *world_size, const int *full_dock){  
  if ( (*full_dock  % *world_size) == 0){
    *n_dock = (int) *full_dock  / *world_size;
    *n_dock_root = *n_dock;
  }else{
    *n_dock = (int) *full_dock  / *world_size;
    *n_dock_root = *n_dock + 1;
  }

}

int main(int argc, char const *argv[]) {

 
  MPI_Init(NULL, NULL);	

  //Creating mpi struct types
/*************  mpi_input_parameters_t ***************************/
  const int nitems=7;
  int blocklengths[nitems] = {MAX_PATH, MAX_PATH, MAX_PATH, MAX_PATH, MAX_PATH, MAX_PATH_FILE_NAME, MAX_PATH_FILE_NAME};
  MPI_Datatype types[nitems] = {MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR};  
  MPI_Aint offsets[nitems];

  offsets[0] = offsetof(input_parameters_t, local_execute);
  offsets[1] = offsetof(input_parameters_t, path_receptors);
  offsets[2] = offsetof(input_parameters_t, path_compounds);
  offsets[3] = offsetof(input_parameters_t, path_out);
  offsets[4] = offsetof(input_parameters_t, path_log);
  offsets[5] = offsetof(input_parameters_t, config_vina);
  offsets[6] = offsetof(input_parameters_t, vina_program);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_input_parameters_t);
  MPI_Type_commit(&mpi_input_parameters_t); 
/*************  mpi_input_parameters_t end ***************************/

/*************  mpi_docking_t ***************************/  
  const int nitems_dock=2;
  int blocklengths_dock[nitems_dock] = {MAX_PATH, MAX_PATH};
  MPI_Datatype types_dock[nitems_dock] = {MPI_CHAR, MPI_CHAR};  
  MPI_Aint offsets_dock[nitems_dock];

  offsets_dock[0] = offsetof(docking_t, receptor);
  offsets_dock[1] = offsetof(docking_t, compound);

  MPI_Type_create_struct(nitems_dock, blocklengths_dock, offsets_dock, types_dock, &mpi_docking_t);
  MPI_Type_commit(&mpi_docking_t); 
/*************  mpi_docking_t end ***************************/


  int root = 0;
  const int tag_docking = 1;
  int number_dock = -1;
  int number_dock_root = -1;      
  int world_size;
  int world_rank;  

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// It is assuming at least 2 processes for performing the Virtual Screening
  if (world_size < 2) {
    fatal_error("World size must be greater than 1 for performing the Virtual Screening\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  //Preparing input parameters data structure
  input_parameters_t *param=NULL;
  param = (input_parameters_t*)malloc(sizeof(input_parameters_t));    
  //loading input parameters files
  if (world_rank == root) {    
    load_parameters_from_file(param,argv[1]);        
  }
  //Broadcast input parameters
  MPI_Bcast(param, 1, mpi_input_parameters_t, root, MPI_COMM_WORLD);  

  //Preparing the docking number for each proecess
  if (world_rank == root) {    
    int full_dock;
    full_dock = get_number_docking(argv[2]);    
    get_number_docking(&number_dock, &number_dock_root, &world_size, &full_dock);     
  } 
  //Broadcast the docking number for each proecess
  MPI_Bcast(&number_dock, 1, MPI_INT, root, MPI_COMM_WORLD);

  number_dock = 2; //REMOVER
  docking_t *v_docking = NULL;
  v_docking = (docking_t*)malloc(number_dock*sizeof(docking_t));

  //Guarantee that all process received all data
  MPI_Barrier(MPI_COMM_WORLD);

  //Preparing buffer to be sent
  docking_t  *buff = NULL;
  int buffer_dock;
  buffer_dock = sizeof(docking_t)*number_dock*MPI_BSEND_OVERHEAD;
  buff = (docking_t *) malloc(buffer_dock);

  //Sending data for performing the Virtual Screening
  if (world_rank == root){    
    //docking_t *docking_root = NULL;
    FILE *f_dock=NULL;    
    char *line=NULL;
    int num_line_ref;
      
    //docking_root = (docking_t*)malloc(number_dock_root*sizeof(docking_t));    
    line = (char*)malloc(MAX_LINE_FILE);

    f_dock = open_file(argv[2], fREAD);
    //Ignoring first line of file
    fgets(line, MAX_LINE_FILE, f_dock);
    num_line_ref = -1;    
    while (fgets(line, MAX_LINE_FILE, f_dock) != NULL){
      num_line_ref = num_line_ref + 1;
      set_receptor_compound(v_docking[num_line_ref].receptor, v_docking[num_line_ref].compound, line);
    }
    free(line);
    fclose(f_dock);
    MPI_Buffer_attach(buff, buffer_dock);
    MPI_Bsend(v_docking, number_dock, mpi_docking_t, 1, tag_docking, MPI_COMM_WORLD);
  }else{
    MPI_Status status;
    int i;
    initialize_vina_execution();
    MPI_Recv(v_docking, number_dock, mpi_docking_t, root,tag_docking, MPI_COMM_WORLD, &status);
    for (i = 0; i < number_dock; i++){
      printf("%s %s \n", v_docking[i].receptor, v_docking[i].compound);
      call_vina(param, &v_docking[i]);
    }
    finish_vina_execution();
    printf("REMOVER O VALOR number_dock QUE ESTA 2. REPENSAR EM COMO ENVIAR \n");        
  }
  
  MPI_Buffer_detach(buff, &buffer_dock);

  deAllocate_docking(v_docking);
  deAllocateload_parameters(param);
  
  MPI_Type_free(&mpi_docking_t);  
  MPI_Type_free(&mpi_input_parameters_t);  

  MPI_Finalize();
	return 0;
}