#ifndef OLD_DOCKING_H
#define OLD_OLD_DOCKING_H

#include "docking_type.h"

void deAllocate_docking(docking_t *d);
int get_number_docking(const char *file_name);
void set_receptor_compound(char *receptor, char *compound, char *source);
#endif

