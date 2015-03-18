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
#include "load_parameters.h"
#include "mpi_parameters_type.h"
#include "mpi_owner.h"

int main(int argc, char const *argv[]) {

 
  MPI_Init(NULL, NULL);	

  //Creating mpi struct types
/*************  mpi_input_parameters_t ***************************/
  const int nitems=3;
  int blocklengths[nitems] = {MAX_PATH,MAX_PATH, MAX_PATH};
  MPI_Datatype types[nitems] = {MPI_CHAR, MPI_CHAR, MPI_CHAR};  
  MPI_Aint offsets[nitems];

  offsets[0] = offsetof(input_parameters_t, local_execute);
  offsets[1] = offsetof(input_parameters_t, path_receptors);
  offsets[2] = offsetof(input_parameters_t, path_compounds);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_input_parameters_t);
  MPI_Type_commit(&mpi_input_parameters_t); 
/*************  mpi_input_parameters_t end ***************************/

  int root = 0;
  int world_size;
  int world_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// It is assuming at least 2 processes for performing the Virtual Screening
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for performing the Virtual Screening\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  //Preparing input parameters data structure
  input_parameters_t *param=NULL;
  param = (input_parameters_t*)malloc(sizeof(input_parameters_t));    
  //loading input parameters
  if (world_rank == 0) {    
    load_parameters_from_file(param,argv[1]);        
  }

  MPI_Bcast(param, 1, mpi_input_parameters_t, root, MPI_COMM_WORLD);

  show_parameters(param);
/*
  else{
    MPI_Status status;
    input_parameters_t *param_r=NULL;
    param_r = (input_parameters_t*)malloc(sizeof(input_parameters_t));        
    //initialize_parameters(param_r);
    //MPI_Irecv(param, 1, mpi_input_parameters_t, 0, 0, MPI_COMM_WORLD, &request);
    MPI_Recv(param_r, 1, mpi_input_parameters_t, root, 0, MPI_COMM_WORLD, &status);
    //MPI_Wait(&request, &status); 
    show_parameters(param_r);
    deAllocateload_parameters(param_r);
  }
*/
  deAllocateload_parameters(param);
  MPI_Type_free(&mpi_input_parameters_t);
  MPI_Finalize();
	return 0;
}