#ifndef OLD_FUTIL_H
#define OLD_FUTIL_H

#include "enums.h"

typedef struct s_owner_file{	
	// Name of file
	char *file_name;
}owner_file_t;

owner_file_t *allocate_owner_file(const int *size);
void desallocate_owner_file(owner_file_t *aux, const int *size);
void show_owner_file(const owner_file_t *aux, const int *size);
char *path_join_file(const char *path, const char *f);
FILE *open_file(__const char *filename, mode_files_t mode);
boolean_t file_is_empty(FILE *file_aux);
boolean_t check_exists_file(const char *path_file_name);
long get_file_size(FILE *fp);
char *get_last_line(const char *filepath);
const char *get_filename_ext(const char *filename);
int how_many_files_directory_by_extension(const char *path, const char *ext);
void insert_files_directory_by_extension(owner_file_t *file_names, const char *path, const char *ext);
void set_base_file_name(char *filename_no_ext, const char *filename_source, const char *ext);
void add_ext_file_name(char *filename, const char *ext);

#endif
