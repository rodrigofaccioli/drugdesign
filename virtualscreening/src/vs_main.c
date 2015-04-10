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
#include "execution_information.h"

int main(int argc, char *argv[]) {

 
  MPI_Init(&argc, &argv);	

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

  int root = 0;
  const int tag_docking = 1;
  int number_dock = -1;
  int number_dock_root = -1;      
  int world_size;
  int world_rank;
  int nthreads;
  double started_time, finished_time;
  time_t started_date, finished_date;
  MPI_Request request_dock;  

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
    full_dock = get_number_docking_from_file(argv[2]);
    set_number_docking(&number_dock, &number_dock_root, &world_size, &full_dock);     
  } 
  //Broadcast the docking number for each proecess
  MPI_Bcast(&number_dock, 1, MPI_INT, root, MPI_COMM_WORLD);
  
  docking_t *v_docking = NULL;
  v_docking = (docking_t*)malloc(number_dock*sizeof(docking_t));

  //Guarantee that all process received all data
  MPI_Barrier(MPI_COMM_WORLD);

  //Decompose all docking against all process
  docking_t *docking_root = NULL;
  if (world_rank == root){
    FILE *f_dock=NULL;    
    char *line=NULL;
    int num_line_ref;
    int dest;    

    docking_root = (docking_t*)malloc(number_dock_root*sizeof(docking_t));
    //saving information of virtual screening execution
    nthreads = 1;//omp_get_num_threads();
    save_information(param->local_execute, &world_size,  &number_dock, &number_dock_root, &nthreads);  

    line = (char*)malloc(MAX_LINE_FILE);
    f_dock = open_file(argv[2], fREAD);
    //Ignoring first  line of file
    fgets(line, MAX_LINE_FILE, f_dock);
    //Obtaining docking that will be executed in root
    for (num_line_ref = 0; num_line_ref < number_dock_root; num_line_ref++){
      fgets(line, MAX_LINE_FILE, f_dock);
      set_receptor_compound(docking_root[num_line_ref].receptor, docking_root[num_line_ref].compound, line);
    }    
    dest = 0;
    save_file_docking_from_array(docking_root, &number_dock_root, param->local_execute, &dest);        
    //Creating for all process    
    for (dest = 1; dest < world_size; dest++){
      num_line_ref = 0;      
      do{
        fgets(line, MAX_LINE_FILE, f_dock);
        set_receptor_compound(v_docking[num_line_ref].receptor, v_docking[num_line_ref].compound, line);
        num_line_ref = num_line_ref + 1;
      }while ( num_line_ref < number_dock  );     
      save_file_docking_from_array(v_docking, &number_dock, param->local_execute, &dest);      
    } 
    fclose(f_dock);
    free(line);

  }
  //Guarantee that all process received all data
  MPI_Barrier(MPI_COMM_WORLD);

  //Call Docking from all process
  if (world_rank == root){
    started_time = MPI_Wtime();
    started_date = time(NULL);

    //printf("Rank %d num_dock %d \n", world_rank, number_dock_root);
    load_docking_from_file(docking_root, &number_dock_root, param->local_execute, &world_rank);
    int i;
    initialize_vina_execution();    
    //printf("rank %d nthreads %d\n", world_rank, nthreads);
    for (i = 0; i < number_dock_root; i++){
      call_vina(param, &docking_root[i]);
    }    
    finish_vina_execution();    
    deAllocate_docking(docking_root);
  }else{
    //printf("Rank %d num_dock %d \n", world_rank, number_dock);
    load_docking_from_file(v_docking, &number_dock, param->local_execute, &world_rank);
    int i;
    initialize_vina_execution();    
    //printf("rank %d nthreads %d\n", world_rank, nthreads);      
    for (i = 0; i < number_dock; i++){
      call_vina(param, &v_docking[i]);
    }
    finish_vina_execution(); 
  }

  MPI_Barrier(MPI_COMM_WORLD);
  finished_time = MPI_Wtime();
  finished_date = time(NULL);

  if (world_rank == root){
    saving_time_execution(param->local_execute, &finished_time, &started_time, &finished_date, &started_date);    
  }

  MPI_Type_free(&mpi_input_parameters_t);  

  deAllocate_docking(v_docking);
  deAllocateload_parameters(param);

  MPI_Finalize();
	return 0;
}