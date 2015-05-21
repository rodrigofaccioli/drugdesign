/*
 * Routines to performe virtual screening using mpi by
 * apdaptative task distribution
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
#include "docking_log.h"

int main(int argc, char *argv[]) {

 
  MPI_Init(&argc, &argv);	

  //Creating mpi struct types
/*************  mpi_input_parameters_t ***************************/
  const int nitems=8;
  int blocklengths[nitems] = {MAX_PATH, MAX_PATH, MAX_PATH, MAX_PATH, MAX_PATH, MAX_PATH_FILE_NAME, MAX_PATH_FILE_NAME, MAX_PATH_FILE_NAME};
  MPI_Datatype types[nitems] = {MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR};  
  MPI_Aint offsets[nitems];

  offsets[0] = offsetof(input_parameters_t, local_execute);
  offsets[1] = offsetof(input_parameters_t, path_receptors);
  offsets[2] = offsetof(input_parameters_t, path_compounds);
  offsets[3] = offsetof(input_parameters_t, path_out);
  offsets[4] = offsetof(input_parameters_t, path_log);
  offsets[5] = offsetof(input_parameters_t, config_vina);
  offsets[6] = offsetof(input_parameters_t, vina_program);
  offsets[7] = offsetof(input_parameters_t, compound_database);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_input_parameters_t);
  MPI_Type_commit(&mpi_input_parameters_t); 
/*************  mpi_input_parameters_t end ***************************/

  int root = 0;
  const int tag_docking = 1;
  const int tag_finsihed_docking = 2;
  int *v_number_dock = NULL;
  char *path_file_log = NULL;
  int number_dock;
  int world_size;
  int world_rank;
  int nthreads;
  MPI_Request request_dock;  
  MPI_Status status;

  char **all_docking = NULL; //store whole docking file
  int *dock_dist = NULL; //represents all docking number 
  int *doc_ref_root = NULL; //represents reference docking
  docking_t *docking_rank = NULL; //it will be performed in each rank
  char *task_from_root = NULL; //represents char that was received by root 
  char *task_executed_by_rank = NULL; //represents char that was executed by rank

  double started_time_global, finished_time_global;
  time_t started_date_global, finished_date_global;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  //Preparing input parameters data structure
  input_parameters_t *param=NULL;
  param = (input_parameters_t*)malloc(sizeof(input_parameters_t));    
  //loading input parameters files
  if (world_rank == root) {    
    load_parameters_from_file(param,argv[1]);        
  }
  //Broadcast input parameters
  MPI_Bcast(param, 1, mpi_input_parameters_t, root, MPI_COMM_WORLD);  

  //Guarantee that all process received all data
  MPI_Barrier(MPI_COMM_WORLD);

  //Preparing to send first docking from root to all ranking  
  if (world_rank == root){
    started_time_global = MPI_Wtime();
    started_date_global = time(NULL);

    FILE *f_dock=NULL;    
    char *line=NULL;
    int i;
    int num_line_ref;
    int rank_ref;

    dock_dist = (int*)malloc(sizeof(int));
    *dock_dist = get_number_docking_from_file(argv[2]);
    line = (char*)malloc(MAX_LINE_FILE);
    //allocating all_docking
    all_docking = (char**)malloc(*dock_dist*sizeof(char*) );
    for (i = 0; i < *dock_dist; i++){
      all_docking[i] = (char*)malloc(sizeof(char)*MAX_LINE_FILE);
    } 

    //Opening docking file and storing into all_docking
    f_dock = open_file(argv[2], fREAD);
    //Ignoring first line of file, because the number of docking we are know at dock_dist variable
    fgets(line, MAX_LINE_FILE, f_dock);
    //Storing whole docking file into all_docking
    num_line_ref=0;
    i = -1;
    while (num_line_ref < *dock_dist){
      i = i + 1;
      fgets(line, MAX_LINE_FILE, f_dock);
      strcpy(all_docking[i], line);
      num_line_ref = num_line_ref + 1;
    }
    fclose(f_dock);
    free(line);

    //sending first docking
    doc_ref_root = (int*)malloc(sizeof(int));
    *doc_ref_root = 0;
    for (rank_ref = 1; rank_ref < world_size; rank_ref++){
      MPI_Send(all_docking[*doc_ref_root], strlen(all_docking[*doc_ref_root]), 
        MPI_CHAR, rank_ref, tag_docking, MPI_COMM_WORLD);
      *doc_ref_root = *doc_ref_root + 1;
    }
    //saving information of virtual screening execution
    nthreads = 1;//omp_get_num_threads();
    save_information_adaptive(param->local_execute, &world_size, dock_dist, &nthreads);

  }else{
    //Allocating char that receives task from root
    task_from_root = (char*)malloc(MAX_LINE_FILE);
    //Allocating docking_t which will perform docking
    docking_rank = (docking_t*)malloc(sizeof(docking_t));
  }

  //Guarantee that all process received the first docking to start
  MPI_Barrier(MPI_COMM_WORLD);

  //starting all ranking to work
  if (world_rank != root){
    time_t f_time_docking, s_time_docking;
    char *docking_execution_info; //represents execution information of docking

    docking_execution_info = (char*)malloc(sizeof(char)*MAX_LINE_FILE);
    initialize_vina_execution();

    MPI_Recv(task_from_root, MAX_LINE_FILE, MPI_CHAR, root, tag_docking, 
      MPI_COMM_WORLD, &status);
    if (strcmp(task_from_root, "") != 0 ){
      set_receptor_compound(docking_rank->receptor, 
            docking_rank->compound,
            &docking_rank->num_torsion_angle,
            &docking_rank->num_atom,          
            task_from_root);          
      s_time_docking = time(NULL);
      call_vina(param, docking_rank);
      f_time_docking = time(NULL);            
      // Sending to root docking executing information
      set_log_file_line(docking_execution_info, docking_rank, &f_time_docking, &s_time_docking);
      MPI_Send(docking_execution_info, strlen(docking_execution_info), 
        MPI_CHAR, root, tag_finsihed_docking, MPI_COMM_WORLD);
      //Receive next either next docking or "" to finish
      MPI_Recv(task_from_root, MAX_LINE_FILE, MPI_CHAR, root, tag_docking, 
        MPI_COMM_WORLD, &status);
    }    
    free(docking_execution_info);
    finish_vina_execution();    
    deAllocate_docking(docking_rank);
    free(task_from_root);    
  }else{
    //root is watting messages from any rank    
    char *log_file_name = NULL;
    char *finished_docking = NULL;
    int received_docking;
    task_executed_by_rank = (char*)malloc(MAX_LINE_FILE);

    //Preparing log file
    log_file_name = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
    set_log_file_name_adaptive(log_file_name);
    path_file_log = path_join_file(param->local_execute, log_file_name);
    free(log_file_name);
    initialize_log(path_file_log);

    received_docking = 0;
    finished_docking = (char*)malloc(sizeof(char)*4);
    strcpy(finished_docking, " ");
    do{
      //Receiving docking from rank
      MPI_Recv(task_executed_by_rank, MAX_LINE_FILE, MPI_CHAR, MPI_ANY_SOURCE,
        tag_finsihed_docking, MPI_COMM_WORLD, &status);
      received_docking = received_docking + 1;
      save_log_by_line(path_file_log,task_executed_by_rank);

      if (*doc_ref_root == *dock_dist){
        MPI_Send(finished_docking, strlen(finished_docking), 
          MPI_CHAR, 1, tag_docking, MPI_COMM_WORLD);
      }else{
        printf("ENVIAR NOVO DOCKING PRO RANK\n");
      }
    }while (received_docking < *dock_dist);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  finished_time_global = MPI_Wtime();
  finished_date_global = time(NULL);

  if (world_rank == root){
    //Desalocating
    int i;
    for (i = 0; i < *dock_dist; i++){
      free(all_docking[i]);
    }
    free(all_docking);    
    free(doc_ref_root);
    free(dock_dist);  
    free(task_executed_by_rank);
    free(path_file_log);
    
    saving_time_execution(param->local_execute, &finished_time_global, &started_time_global, &finished_date_global, &started_date_global); 
  }

  MPI_Type_free(&mpi_input_parameters_t);  
  
  deAllocateload_parameters(param);

  MPI_Finalize();
	return 0;
}