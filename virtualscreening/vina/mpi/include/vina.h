#ifndef OLD_VINA_H
#define OLD_VINA_H

#include "parameters_type.h"
#include "docking_type.h"

void initialize_vina_execution();
void finish_vina_execution();
void call_vina(const input_parameters_t *param, const docking_t *docking);
void call_prepare_ligand(const input_parameters_t *param, const char *fmol2, const char *fpdbqt);

#endif