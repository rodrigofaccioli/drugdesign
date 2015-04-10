#ifndef OLD_EXECUTION_INFORMATION_H
#define OLD_EXECUTION_INFORMATION_H

#include <time.h>

void save_information(const char *local_execute, const int *num_proc, 
	const int *num_dock, const int *num_dock_root, const int *nthreads);
void saving_time_execution(const char *local, const double *finished, const double *started,
	const time_t *f_time, const time_t *s_time);
#endif
