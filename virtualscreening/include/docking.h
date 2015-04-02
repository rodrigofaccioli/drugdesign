#ifndef OLD_DOCKING_H
#define OLD_OLD_DOCKING_H

#include "docking_type.h"

void deAllocate_docking(docking_t *d);
int get_number_docking_from_file(const char *file_name);
void set_receptor_compound(char *receptor, char *compound, char *source);
void save_file_docking_from_array(const docking_t *v_doc, const int *num_doc, 
  const char *local_execute, const int *suf);
void load_docking_from_file(docking_t *v_doc, const int *num_dock, 
  const char *local, const int *suf);
#endif

