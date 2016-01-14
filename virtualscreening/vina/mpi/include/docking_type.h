#ifndef OLD_DOCKING_TYPE_H
#define OLD_DOCKING_TYPE_H

#include "defines.h"

typedef struct s_docking{
	//Indicates receptor name
	char receptor[MAX_FILE_NAME];
	//Indicates compound (ligand) name
	char compound[MAX_FILE_NAME];
	//Number of torsion angle of ligand
	int num_torsion_angle;
	//Number of atoms of ligand
	int num_atom;
}docking_t;

#endif