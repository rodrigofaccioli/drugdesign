#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_PATH_FILE_NAME 500

int main(int argc, char *argv[]){
	
	typedef struct{
		float x,y,z;
		char res_name[10];
		char atm_name[10];
		char atm_type[10];
		int res_number;
		int prev;
		int h;
		int connectivity;
	} atom;
	
	// variable declaration
	float distance, distance_cutoff, angle_cutoff;
	float dist_prev_donor, dist_donor_h, dist_h_acceptor, dist_prev_h, dist_prev_acceptor, dist_donor_acceptor;
	float angle_prev_donor_h, angle_prev_donor_acceptor, angle_h_donor_acceptor;
	//char lig_filename[100], rec_filename[100];
	int total_atm_lig, total_atm_rec, atm_rec, atm_lig, atm;
	int total_lig_h_donor, total_lig_h_acceptor, total_rec_h_donor, total_rec_h_acceptor;
	int bonded_h, h_index ;

	char *lig_filename, *rec_filename, *output_filename, *f_path_temporary_rec_no,	
	     *f_path_temporary_lig_no, *f_path_temporary_rec_h, *f_path_temporary_lig_h;

	FILE *input_lig, *input_rec, *f_output_filename;

	atom *lig, *rec;
	atom *lig_h_donor, *lig_h_acceptor, *rec_h_donor, *rec_h_acceptor;

	//Allocating file names
	lig_filename = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );
	rec_filename = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );
	output_filename = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );
	f_path_temporary_rec_no = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );
	f_path_temporary_lig_no = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );
	f_path_temporary_rec_h = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );
	f_path_temporary_lig_h = (char*) malloc( MAX_PATH_FILE_NAME * sizeof(char) );

		
	// reading parameters from input command line	
	strcpy(rec_filename, argv[1]);	
	total_atm_rec = atoi(argv[2]);	
	strcpy(lig_filename, argv[3]);
	total_atm_lig =  atoi(argv[4]);	
	distance_cutoff = atof(argv[5]);
	angle_cutoff = atof(argv[6]);
	strcpy(output_filename, argv[7]);
	strcpy(f_path_temporary_rec_no, argv[8] );
	strcpy(f_path_temporary_lig_no, argv[9] );
	strcpy(f_path_temporary_rec_h, argv[10] );
	strcpy(f_path_temporary_lig_h, argv[11] );	
	
	// Setting output filename
	f_output_filename = fopen(output_filename,"w");

	// dynamic allocation
	lig = (atom*) malloc( (total_atm_lig+1) * sizeof(atom) );
	
	lig_h_acceptor = (atom*) malloc( (total_atm_lig+1) * sizeof(atom) );
	lig_h_donor = (atom*) malloc( (total_atm_lig+1) * sizeof(atom) );
	
	rec = (atom*) malloc( (total_atm_rec+1) * sizeof(atom) );
	
	rec_h_acceptor = (atom*) malloc( (total_atm_rec+1) * sizeof(atom) );
	rec_h_donor = (atom*) malloc( (total_atm_rec+1) * sizeof(atom) );
	
	
	// reading input files
	input_rec = fopen( f_path_temporary_rec_no , "r" );
	atm_rec = 1;
	while( fscanf( input_rec, "%s %d %s %f %f %f %s",
		&rec[atm_rec].res_name[0], 
		&rec[atm_rec].res_number, 
		&rec[atm_rec].atm_name[0],
		&rec[atm_rec].x, 
		&rec[atm_rec].y, 
		&rec[atm_rec].z, 
		&rec[atm_rec].atm_type[0]
	)!=EOF ) atm_rec++ ;
	fclose(input_rec);
	
	input_lig = fopen( f_path_temporary_lig_no , "r" );
	atm_lig = 1;
	while( fscanf( input_lig , "%s %f %f %f %s", 
		&lig[atm_lig].atm_name[0], 
		&lig[atm_lig].x, 
		&lig[atm_lig].y, 
		&lig[atm_lig].z, 
		&lig[atm_lig].atm_type[0]
	)!=EOF) atm_lig++ ;
	fclose(input_lig);



/*

AutoDock Atom Types for H, N and O:

H  Non H-bonding Hydrogen
HD Donor 1 H-bond Hydrogen
HS Donor S Spherical Hydrogen

N  Non H-bonding Nitrogen
NA Acceptor 1 H-bond Nitrogen
NS Acceptor S Spherical Nitrogen

OA Acceptor 2 H-bonds Oxygen
OS Acceptor S Spherical Oxygen

*/
	
	// Classify all ligand O and N as H donors and acceptors
	total_lig_h_acceptor = 0;
	total_lig_h_donor = 0;
	
	for( atm_lig=1 ; atm_lig<=total_atm_lig ; atm_lig++ ){
		
		// If the atom type is any N or O type, we must check if there is a H bonded to it.
		// If there is, this atom is a HB donor
		if(  (strcmp(lig[atm_lig].atm_type,"N")==0) ||
			(strcmp(lig[atm_lig].atm_type,"NA")==0) ||
			(strcmp(lig[atm_lig].atm_type,"NS")==0) ||
			(strcmp(lig[atm_lig].atm_type,"OA")==0) ||
			(strcmp(lig[atm_lig].atm_type,"OS")==0)
		){
			bonded_h = 0;	// Number of Hs covalently bonded to the current O or N atom
			
			// For each hydrogen ...
			for( atm=1 ; atm<=total_atm_lig; atm++ ){
				if( (strcmp(lig[atm].atm_type,"HD")==0) ||
				    (strcmp(lig[atm].atm_type,"HS")==0) 
				){
					distance=sqrt(
						pow( lig[atm_lig].x - lig[atm].x ,2) + 
						pow( lig[atm_lig].y - lig[atm].y ,2) + 
						pow( lig[atm_lig].z - lig[atm].z ,2)
					);
					// If the distance between the current O or N atom and this H is <= 1.25 A
					if( distance <= 1.25 ){
						bonded_h++;	// Increase bounded H counter
						h_index = atm ;	// Store the H atom number
					}
				}
			}
			// If there is 1 or more H bonded to the current N/O atom, Classify it as a HB donor.
			if( bonded_h >= 1 ){
				total_lig_h_donor++;
				lig_h_donor[total_lig_h_donor].x = lig[atm_lig].x;
				lig_h_donor[total_lig_h_donor].y = lig[atm_lig].y;
				lig_h_donor[total_lig_h_donor].z = lig[atm_lig].z;
				strcpy( lig_h_donor[total_lig_h_donor].atm_name , lig[atm_lig].atm_name );
				
				lig_h_donor[total_lig_h_donor].h = h_index ;
				
				// Now, it is necessary to find which heavy atom is bound to the donor.
				// This is stored in lig_h_donor[].prev as the index of such atom.
				// e.g, the heavy atom bound to lig_h_donor[i].atm_name is lig[ lig_h_donor[i].prev ].atm_name
				for( atm=1 ; atm<=total_atm_lig; atm++ ){
					
					if(  (strcmp(lig[atm].atm_type,"HD")!=0) &&
						(strcmp(lig[atm].atm_type,"HS")!=0) &&
						(strcmp(lig[atm].atm_type,"H")!=0)
					){
						distance=sqrt(
							pow( lig[atm_lig].x - lig[atm].x ,2) + 
							pow( lig[atm_lig].y - lig[atm].y ,2) + 
							pow( lig[atm_lig].z - lig[atm].z ,2)
						);
						if( (distance<=1.90) && (distance>=0.3) ){
							lig_h_donor[total_lig_h_donor].prev = atm;	// >=0.3 to avoid bonding with itself
							lig_h_donor[total_lig_h_donor].connectivity++;
						}
					}
					
				}
				
			}
			
			
			
			// Types NA, NS, OA and OS are always acceptors
			if( strcmp( lig[atm_lig].atm_type, "N" )!=0 ){
				total_lig_h_acceptor++;
				lig_h_acceptor[total_lig_h_acceptor].x = lig[atm_lig].x;
				lig_h_acceptor[total_lig_h_acceptor].y = lig[atm_lig].y;
				lig_h_acceptor[total_lig_h_acceptor].z = lig[atm_lig].z;
				strcpy( lig_h_acceptor[total_lig_h_acceptor].atm_name , lig[atm_lig].atm_name );
			}	
		}
	}
		
		
	
	
	
	
	// Classify all receptor O and N as H donors and acceptors
	total_rec_h_acceptor = 0;
	total_rec_h_donor = 0;
	
	for( atm_rec=1 ; atm_rec<=total_atm_rec ; atm_rec++ ){
		
		// If the atom type is any N or O type, we must check if there is a H bonded to it.
		// If there is, this atom is a HB donor
		if(  (strcmp(rec[atm_rec].atm_type,"N")==0) ||
			(strcmp(rec[atm_rec].atm_type,"NA")==0) ||
			(strcmp(rec[atm_rec].atm_type,"NS")==0) ||
			(strcmp(rec[atm_rec].atm_type,"OA")==0) ||
			(strcmp(rec[atm_rec].atm_type,"OS")==0)
		){
			bonded_h = 0;	// Number of Hs covalently bonded to the current O or N atom
			for( atm=1 ; atm<=total_atm_rec; atm++ ){
				if( (strcmp(rec[atm].atm_type,"HD")==0) ||
				    (strcmp(rec[atm].atm_type,"HS")==0) 
				){
					distance=sqrt(
						pow( rec[atm_rec].x - rec[atm].x ,2) + 
						pow( rec[atm_rec].y - rec[atm].y ,2) + 
						pow( rec[atm_rec].z - rec[atm].z ,2)
					);
					if( distance <= 1.25 ){
						bonded_h++;
						h_index = atm ;
					}
				}
			}
			if( bonded_h >= 1 ){
				total_rec_h_donor++;
				rec_h_donor[total_rec_h_donor].x = rec[atm_rec].x;
				rec_h_donor[total_rec_h_donor].y = rec[atm_rec].y;
				rec_h_donor[total_rec_h_donor].z = rec[atm_rec].z;
				rec_h_donor[total_rec_h_donor].res_number = rec[atm_rec].res_number;
				strcpy( rec_h_donor[total_rec_h_donor].atm_name , rec[atm_rec].atm_name );
				strcpy( rec_h_donor[total_rec_h_donor].res_name , rec[atm_rec].res_name );
				
				rec_h_donor[total_rec_h_donor].h = h_index ;
				
				// Now, it is necessary to find which heavy atom is bound to the donor.
				// This is stored in rec_h_donor[].prev as the index of such atom.
				// e.g, the heavy atom bound to rec_h_donor[i].atm_name is rec[ rec_h_donor[i].prev ].atm_name
				for( atm=1 ; atm<=total_atm_rec; atm++ ){
					
					if(  (strcmp(rec[atm].atm_type,"HD")!=0) &&
						(strcmp(rec[atm].atm_type,"HS")!=0) &&
						(strcmp(rec[atm].atm_type,"H")!=0)
					){
						distance=sqrt(
							pow( rec[atm_rec].x - rec[atm].x ,2) + 
							pow( rec[atm_rec].y - rec[atm].y ,2) + 
							pow( rec[atm_rec].z - rec[atm].z ,2)
						);
						if( (distance<=1.90) && (distance>=0.3) ){
							rec_h_donor[total_rec_h_donor].prev = atm;	// >=0.3 to avoid bonding with itself
							rec_h_donor[total_rec_h_donor].connectivity++;
						}
					}
					
				}
				
			}
			
			
			// Types NA, NS, OA and OS are always acceptors
			if( strcmp( rec[atm_rec].atm_type, "N" )!=0 ){
				total_rec_h_acceptor++;
				rec_h_acceptor[total_rec_h_acceptor].x = rec[atm_rec].x;
				rec_h_acceptor[total_rec_h_acceptor].y = rec[atm_rec].y;
				rec_h_acceptor[total_rec_h_acceptor].z = rec[atm_rec].z;
				rec_h_acceptor[total_rec_h_acceptor].res_number = rec[atm_rec].res_number;
				strcpy( rec_h_acceptor[total_rec_h_acceptor].atm_name , rec[atm_rec].atm_name );
				strcpy( rec_h_acceptor[total_rec_h_acceptor].res_name , rec[atm_rec].res_name );
			}
			
			
		}
	}
	
	

	

	
	
	// Find the ligand-receptor h-bonds in which the ligand is donating H
	for( atm_lig=1 ; atm_lig<=total_lig_h_donor ; atm_lig++ ){
		for( atm_rec=1 ; atm_rec<=total_rec_h_acceptor; atm_rec++ ){
			
			dist_prev_donor = sqrt(
				pow( lig[ lig_h_donor[atm_lig].prev ].x - lig_h_donor[atm_lig].x , 2) +
				pow( lig[ lig_h_donor[atm_lig].prev ].y - lig_h_donor[atm_lig].y , 2) +
				pow( lig[ lig_h_donor[atm_lig].prev ].z - lig_h_donor[atm_lig].z , 2)
			);
			
			dist_donor_h = sqrt(
				pow( lig_h_donor[atm_lig].x - lig[ lig_h_donor[atm_lig].h ].x , 2) +
				pow( lig_h_donor[atm_lig].y - lig[ lig_h_donor[atm_lig].h ].y , 2) +
				pow( lig_h_donor[atm_lig].z - lig[ lig_h_donor[atm_lig].h ].z , 2)
			);
			
			dist_prev_h = sqrt(
				pow( lig[ lig_h_donor[atm_lig].prev ].x - lig[ lig_h_donor[atm_lig].h ].x , 2) +
				pow( lig[ lig_h_donor[atm_lig].prev ].y - lig[ lig_h_donor[atm_lig].h ].y , 2) +
				pow( lig[ lig_h_donor[atm_lig].prev ].z - lig[ lig_h_donor[atm_lig].h ].z , 2)

			);
			
			dist_h_acceptor = sqrt(
				pow( lig[ lig_h_donor[atm_lig].h ].x - rec_h_acceptor[atm_rec].x , 2) +
				pow( lig[ lig_h_donor[atm_lig].h ].y - rec_h_acceptor[atm_rec].y , 2) +
				pow( lig[ lig_h_donor[atm_lig].h ].z - rec_h_acceptor[atm_rec].z , 2)

			);
			
			dist_prev_acceptor = sqrt(
				pow( lig[ lig_h_donor[atm_lig].prev ].x - rec_h_acceptor[atm_rec].x , 2) +
				pow( lig[ lig_h_donor[atm_lig].prev ].y - rec_h_acceptor[atm_rec].y , 2) +
				pow( lig[ lig_h_donor[atm_lig].prev ].z - rec_h_acceptor[atm_rec].z , 2)

			);
			
			dist_donor_acceptor = sqrt(
				pow( lig_h_donor[atm_lig].x - rec_h_acceptor[atm_rec].x , 2) +
				pow( lig_h_donor[atm_lig].y - rec_h_acceptor[atm_rec].y , 2) +
				pow( lig_h_donor[atm_lig].z - rec_h_acceptor[atm_rec].z , 2) 
			);
			
			
			angle_prev_donor_h = acos(
				( pow(dist_prev_donor,2) + pow(dist_donor_h,2) - pow(dist_prev_h,2) )
				/ ( 2 * dist_prev_donor * dist_donor_h )
			);
			
			
			angle_prev_donor_acceptor = acos(
				( pow(dist_prev_donor,2) + pow(dist_donor_acceptor ,2) - pow(dist_prev_acceptor,2) )
				/ ( 2 * dist_prev_donor * dist_donor_acceptor )
			);
			
			
			if( lig_h_donor[atm_lig].connectivity >= 2 ){
				angle_h_donor_acceptor = acos(
					( pow(dist_donor_h,2) + pow(dist_donor_acceptor,2) - pow(dist_h_acceptor,2) )
					/ ( 2 * dist_donor_h * dist_donor_acceptor )
				);
			}else{
				angle_h_donor_acceptor = fabs( angle_prev_donor_h - angle_prev_donor_acceptor );
			} 
			
			angle_h_donor_acceptor*=(180.0/3.14155);
			
			
			
			
			if( (dist_donor_acceptor <= distance_cutoff) && (angle_h_donor_acceptor <= angle_cutoff) ){
				fprintf(f_output_filename,"LIG-%s\tdonates_to\t%s-%d %s\t%.1f\t%.1f\t%s %s\n", 
					lig_h_donor[atm_lig].atm_name, 
					rec_h_acceptor[atm_rec].res_name, 
					rec_h_acceptor[atm_rec].res_number,
					rec_h_acceptor[atm_rec].atm_name,
					dist_donor_acceptor,
					angle_h_donor_acceptor,
					rec_filename, 
					lig_filename
				);
			
			}	
		}		
	}
	
	
	
	
	
	// Find the ligand-receptor h-bonds in which the ligand is accepting H
	for( atm_lig=1 ; atm_lig<=total_lig_h_acceptor ; atm_lig++ ){
		for( atm_rec=1 ; atm_rec<=total_rec_h_donor; atm_rec++ ){
			
			dist_prev_donor = sqrt(
				pow( rec[ rec_h_donor[atm_rec].prev ].x - rec_h_donor[atm_rec].x , 2) +
				pow( rec[ rec_h_donor[atm_rec].prev ].y - rec_h_donor[atm_rec].y , 2) +
				pow( rec[ rec_h_donor[atm_rec].prev ].z - rec_h_donor[atm_rec].z , 2)
			);
			
			dist_donor_h = sqrt(
				pow( rec_h_donor[atm_rec].x - rec[ rec_h_donor[atm_rec].h ].x , 2) +
				pow( rec_h_donor[atm_rec].y - rec[ rec_h_donor[atm_rec].h ].y , 2) +
				pow( rec_h_donor[atm_rec].z - rec[ rec_h_donor[atm_rec].h ].z , 2)
			);
			
			dist_prev_h = sqrt(
				pow( rec[ rec_h_donor[atm_rec].prev ].x - rec[ rec_h_donor[atm_rec].h ].x , 2) +
				pow( rec[ rec_h_donor[atm_rec].prev ].y - rec[ rec_h_donor[atm_rec].h ].y , 2) +
				pow( rec[ rec_h_donor[atm_rec].prev ].z - rec[ rec_h_donor[atm_rec].h ].z , 2)

			);
			
			dist_h_acceptor = sqrt(
				pow( rec[ rec_h_donor[atm_rec].h ].x - lig_h_acceptor[atm_lig].x , 2) +
				pow( rec[ rec_h_donor[atm_rec].h ].y - lig_h_acceptor[atm_lig].y , 2) +
				pow( rec[ rec_h_donor[atm_rec].h ].z - lig_h_acceptor[atm_lig].z , 2)

			);
			
			dist_prev_acceptor = sqrt(
				pow( rec[ rec_h_donor[atm_rec].prev ].x - lig_h_acceptor[atm_lig].x , 2) +
				pow( rec[ rec_h_donor[atm_rec].prev ].y - lig_h_acceptor[atm_lig].y , 2) +
				pow( rec[ rec_h_donor[atm_rec].prev ].z - lig_h_acceptor[atm_lig].z , 2)

			);
			
			dist_donor_acceptor = sqrt(
				pow( rec_h_donor[atm_rec].x - lig_h_acceptor[atm_lig].x , 2) +
				pow( rec_h_donor[atm_rec].y - lig_h_acceptor[atm_lig].y , 2) +
				pow( rec_h_donor[atm_rec].z - lig_h_acceptor[atm_lig].z , 2) 
			);
			
			
			angle_prev_donor_h = acos(
				( pow(dist_prev_donor,2) + pow(dist_donor_h,2) - pow(dist_prev_h,2) )
				/ ( 2 * dist_prev_donor * dist_donor_h )
			);
			
			
			angle_prev_donor_acceptor = acos(
				( pow(dist_prev_donor,2) + pow(dist_donor_acceptor ,2) - pow(dist_prev_acceptor,2) )
				/ ( 2 * dist_prev_donor * dist_donor_acceptor )
			);
			
			
			if( rec_h_donor[atm_lig].connectivity >= 2 ){
				angle_h_donor_acceptor = acos(
					( pow(dist_donor_h,2) + pow(dist_donor_acceptor,2) - pow(dist_h_acceptor,2) )
					/ ( 2 * dist_donor_h * dist_donor_acceptor )
				);
			}else{
				angle_h_donor_acceptor = fabs( angle_prev_donor_h - angle_prev_donor_acceptor );
			} 
			
			angle_h_donor_acceptor*=(180.0/3.14155);
			
			
			
			
			if( (dist_donor_acceptor <= distance_cutoff) && (angle_h_donor_acceptor <= angle_cutoff) ){
				fprintf(f_output_filename, "LIG-%s\taccepts_from\t%s-%d %s\t%.1f\t%.1f\t%s %s\n", 
					lig_h_acceptor[atm_lig].atm_name, 
					rec_h_donor[atm_rec].res_name, 
					rec_h_donor[atm_rec].res_number,
					rec_h_donor[atm_rec].atm_name,
					dist_donor_acceptor,
					angle_h_donor_acceptor,
					rec_filename, 
					lig_filename
				);
			}	
		}		
	}
	
	fclose(f_output_filename);
	
	free(lig_filename);
	free(rec_filename);
	free(output_filename);
	free(f_path_temporary_rec_no);
	free(f_path_temporary_lig_no);
	free(f_path_temporary_rec_h);
	free(f_path_temporary_lig_h);

	return 0;
}