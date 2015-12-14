/*
 * Routines to performe virtual screening using mpi with 1 processes
 * These routines were developed by:
 * Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
 * Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com
*/
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <time.h>

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

  int root = 0;  
  int number_dock = -1;
  int number_dock_root = -1;      
  int world_size;
  int world_rank;
  int nthreads;
  int num_line_ref;
  double started_time, finished_time;
  time_t started_date, finished_date;
  FILE *f_dock=NULL;    
  char *line=NULL;
  docking_t *v_docking = NULL;  
  char *aux_fgets = NULL; //used for getting return of fgets

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // It is assuming the number if processes for performing the Virtual Screening must be 1
  if (world_size > 1) {
    fatal_error("World size must be 1 for performing the Virtual Screening\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

	//Preparing input parameters data structure
  input_parameters_t *param=NULL;
  param = (input_parameters_t*)malloc(sizeof(input_parameters_t));    
  load_parameters_from_file(param,argv[1]);        
  
  //saving information of virtual screening execution
  nthreads = omp_get_num_threads();
  number_dock_root = 0;
  number_dock = get_number_docking_from_file(argv[2]);  
  save_information(param->local_execute, &world_size,  &number_dock, &nthreads);  

  started_time = MPI_Wtime();
  started_date = time(NULL);

  v_docking = (docking_t*)malloc(number_dock*sizeof(docking_t));

  line = (char*)malloc(MAX_LINE_FILE);
  f_dock = open_file(argv[2], fREAD);
  //Ignoring first  line of file
  aux_fgets = fgets(line, MAX_LINE_FILE, f_dock);
  //Obtaining docking that will be executed in root  
  for (num_line_ref = 0; num_line_ref < number_dock; num_line_ref++){
      aux_fgets = fgets(line, MAX_LINE_FILE, f_dock);
      set_receptor_compound(v_docking[num_line_ref].receptor, 
        v_docking[num_line_ref].compound,
        &v_docking[num_line_ref].num_torsion_angle,
        &v_docking[num_line_ref].num_atom,
        line);
  }    
  fclose(f_dock);
  free(line);

  int i;
  initialize_vina_execution();      
  for (i = 0; i < number_dock; i++){
    call_vina(param, &v_docking[i]);
  }
  finish_vina_execution(); 

  finished_time = MPI_Wtime();
  finished_date = time(NULL);

  saving_time_execution(param->local_execute, &finished_time, &started_time, &finished_date, &started_date);

  deAllocate_docking(v_docking);
  deAllocateload_parameters(param);

	return 0;
}