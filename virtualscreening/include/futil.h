#ifndef OLD_FUTIL_H
#define OLD_FUTIL_H

#include "enums.h"

char *path_join_file(const char *path, const char *f);

FILE *open_file(__const char *filename, mode_files_t mode);
boolean_t file_is_empty(FILE *file_aux);
boolean_t check_exists_file(const char *path_file_name);
long get_file_size(FILE *fp);
char *get_last_line(const char *filepath);

#endif
