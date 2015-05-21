#ifndef OLD_DOCKING_LOG_H
#define OLD_OLD_DOCKING_LOG_H

#include <time.h>
#include "docking_type.h"

void set_log_file_line(char *line, const docking_t *docking, 
	const time_t *f_time, const time_t *s_time);
void set_log_file_name_adaptive(char *file_name);
void build_log_file_name_from_rank(char *file_name, const int *rank);
void initialize_log(const char *path_file_name);
void save_log_by_line(const char *path_file_name, const char *line);
void save_log(const char *path_file_name, const docking_t *docking, 
	const time_t *f_time, const time_t *s_time);

#endif