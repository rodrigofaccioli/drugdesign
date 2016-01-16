#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[]){
	
	typedef struct{
		float x,y,z;
		char res_name[10];
		char atm_name[10];
		char atm_type[10];
		int res_number;
	} atom;
	
	// variable declaration
	int total_lig_no, total_lig_h, total_rec_no, total_rec_h ;
	int total_lig_h_acceptor, total_lig_h_donor, total_rec_h_acceptor, total_rec_h_donor ;
	int atm_lig_no, atm_lig_h, atm_rec_no, atm_rec_h;
	int bonded_h;
	float distance, distance_cutoff ;
	char lig_filename[200], rec_filename[200];
	
	FILE *input_lig_no, *input_lig_h, *input_rec_no, *input_rec_h;
	
	atom *lig_no, *lig_h, *lig_h_acceptor, *lig_h_donor;
	atom *rec_no, *rec_h, *rec_h_acceptor, *rec_h_donor;

		
	// reading parameters from input command line
	sscanf( argv[1] , "%s", &rec_filename[0] );
	sscanf( argv[2] , "%d", &total_rec_no );
	sscanf( argv[3] , "%d", &total_rec_h );
	sscanf( argv[4] , "%s", &lig_filename[0] );
	sscanf( argv[5] , "%d", &total_lig_no );
	sscanf( argv[6] , "%d", &total_lig_h );
	sscanf( argv[7] , "%f", &distance_cutoff );
	
	
	// dynamic allocation
	lig_no = (atom*) malloc( (total_lig_no+1) * sizeof(atom) );
	lig_h = (atom*) malloc( (total_lig_h+1) * sizeof(atom) );
	
	lig_h_acceptor = (atom*) malloc( (total_lig_no+1) * sizeof(atom) );
	lig_h_donor = (atom*) malloc( (total_lig_no+1) * sizeof(atom) );
	
	
	rec_no = (atom*) malloc( (total_rec_no+1) * sizeof(atom) );
	rec_h = (atom*) malloc( (total_rec_h+1) * sizeof(atom) );
	
	rec_h_acceptor = (atom*) malloc( (total_rec_no+1) * sizeof(atom) );
	rec_h_donor = (atom*) malloc( (total_rec_h+1) * sizeof(atom) );
	
	
	// reading input files
	input_rec_no = fopen( "temporary_rec_no" , "r" );
	atm_rec_no = 1;
	while( fscanf( input_rec_no , "%s %d %s %f %f %f %s",
		&rec_no[atm_rec_no].res_name[0], 
		&rec_no[atm_rec_no].res_number, 
		&rec_no[atm_rec_no].atm_name[0],
		&rec_no[atm_rec_no].x, 
		&rec_no[atm_rec_no].y, 
		&rec_no[atm_rec_no].z, 
		&rec_no[atm_rec_no].atm_type[0]
	)!=EOF ) atm_rec_no++;
	fclose(input_rec_no);
	
	input_lig_no = fopen( "temporary_lig_no" , "r" );
	atm_lig_no=1;
	while( fscanf( input_lig_no , "%s %f %f %f %s", 
		&lig_no[atm_lig_no].atm_name[0], 
		&lig_no[atm_lig_no].x, 
		&lig_no[atm_lig_no].y, 
		&lig_no[atm_lig_no].z, 
		&lig_no[atm_lig_no].atm_type[0]
	)!=EOF) atm_lig_no++;
	fclose(input_lig_no);

	input_rec_h = fopen( "temporary_rec_h" , "r" );
	atm_rec_h = 1;
	while( fscanf( input_rec_no , "%s %d %s %f %f %f %s", 
		&rec_h[atm_rec_h].res_name[0], 
		&rec_h[atm_rec_h].res_number, 
		&rec_h[atm_rec_h].atm_name[0], 
		&rec_h[atm_rec_h].x, 
		&rec_h[atm_rec_h].y, 
		&rec_h[atm_rec_h].z, 
		&rec_h[atm_rec_h].atm_type[0]
	)!=EOF ) atm_rec_h++;
	fclose(input_rec_h);
	
	input_lig_h = fopen( "temporary_lig_h" , "r" );
	atm_lig_h=1;
	while( fscanf( input_lig_h , "%s %f %f %f %s", 
		&lig_h[atm_lig_h].atm_name[0], 
		&lig_h[atm_lig_h].x, 
		&lig_h[atm_lig_h].y, 
		&lig_h[atm_lig_h].z, 
		&lig_h[atm_lig_h].atm_type[0]
	)!=EOF ) atm_lig_h++;
	fclose(input_lig_h);
	
	
	// Classify all ligand O and N as H donors and acceptors
	total_lig_h_acceptor = 0;
	total_lig_h_donor = 0;
	for( atm_lig_no=1 ; atm_lig_no<=total_lig_no ; atm_lig_no++ ){
		
		bonded_h = 0;
		for( atm_lig_h=1 ; atm_lig_h<=total_lig_h; atm_lig_h++ ){
			
			distance=sqrt(
				pow( lig_no[atm_lig_no].x - lig_h[atm_lig_h].x ,2) + 
				pow( lig_no[atm_lig_no].y - lig_h[atm_lig_h].y ,2) + 
				pow( lig_no[atm_lig_no].z - lig_h[atm_lig_h].z ,2)
			);
			// If the distance of a O or N atom to a H that can be donated in a h-bond is < 1.5 A, the O/N atom is a h-bond donor.
			if( distance <= 1.25 ) bonded_h++;
				
		}
		if( bonded_h >= 1 ){
			total_lig_h_donor++;
			lig_h_donor[total_lig_h_donor].x = lig_no[atm_lig_no].x;
			lig_h_donor[total_lig_h_donor].y = lig_no[atm_lig_no].y;
			lig_h_donor[total_lig_h_donor].z = lig_no[atm_lig_no].z;
			strcpy( lig_h_donor[total_lig_h_donor].atm_name , lig_no[atm_lig_no].atm_name );
		}
		// Additionally, if the atom type is not N, it is a h-bond acceptor, i.e. types NA, NS, OA and OS are always acceptors
		if( strcmp( lig_no[atm_lig_no].atm_type, "N" )!=0 ){
			total_lig_h_acceptor++;
			lig_h_acceptor[total_lig_h_acceptor].x = lig_no[atm_lig_no].x;
			lig_h_acceptor[total_lig_h_acceptor].y = lig_no[atm_lig_no].y;
			lig_h_acceptor[total_lig_h_acceptor].z = lig_no[atm_lig_no].z;
			strcpy( lig_h_acceptor[total_lig_h_acceptor].atm_name , lig_no[atm_lig_no].atm_name );
		}	
	}
	
	// Classify all receptor O and N as H donors and acceptors
	total_rec_h_acceptor = 0;
	total_rec_h_donor = 0;
	for( atm_rec_no=1 ; atm_rec_no<=total_rec_no ; atm_rec_no++ ){
		
		bonded_h = 0;
		for( atm_rec_h=1 ; atm_rec_h<=total_rec_h; atm_rec_h++ ){
			
			distance=sqrt( 
				pow( rec_no[atm_rec_no].x - rec_h[atm_rec_h].x ,2) + 
				pow( rec_no[atm_rec_no].y - rec_h[atm_rec_h].y ,2) + 
				pow( rec_no[atm_rec_no].z - rec_h[atm_rec_h].z ,2) );
			// If the distance of a O or N atom to a H that can be donated in a h-bond is < 1.5 A, the O/N atom is a h-bond donor.
			if( distance <= 1.25 ) bonded_h++;
		}
		if( bonded_h >= 1 ){
			total_rec_h_donor++;
			rec_h_donor[total_rec_h_donor].x = rec_no[atm_rec_no].x;
			rec_h_donor[total_rec_h_donor].y = rec_no[atm_rec_no].y;
			rec_h_donor[total_rec_h_donor].z = rec_no[atm_rec_no].z;
			rec_h_donor[total_rec_h_donor].res_number = rec_no[atm_rec_no].res_number;
			strcpy( rec_h_donor[total_rec_h_donor].res_name , rec_no[atm_rec_no].res_name );
			strcpy( rec_h_donor[total_rec_h_donor].atm_name , rec_no[atm_rec_no].atm_name );
		}
		// Additionally, if the atom type is not N, it is a h-bond acceptor, i.e. types NA, NS, OA and OS are always acceptors
		if( strcmp( rec_no[atm_rec_no].atm_type, "N" )!=0 ){
			total_rec_h_acceptor++;
			rec_h_acceptor[total_rec_h_acceptor].x = rec_no[atm_rec_no].x;
			rec_h_acceptor[total_rec_h_acceptor].y = rec_no[atm_rec_no].y;
			rec_h_acceptor[total_rec_h_acceptor].z = rec_no[atm_rec_no].z;
			rec_h_acceptor[total_rec_h_acceptor].res_number = rec_no[atm_rec_no].res_number;
			strcpy( rec_h_acceptor[total_rec_h_acceptor].res_name , rec_no[atm_rec_no].res_name );
			strcpy( rec_h_acceptor[total_rec_h_acceptor].atm_name , rec_no[atm_rec_no].atm_name );
		}	
	}
	
	
	// Find the ligand-receptor h-bonds in which the ligand is donating H
	for( atm_lig_no=1 ; atm_lig_no<=total_lig_h_donor ; atm_lig_no++ ){
		for( atm_rec_no=1 ; atm_rec_no<=total_rec_h_acceptor; atm_rec_no++ ){
			
			distance=sqrt( 
				pow( lig_h_donor[atm_lig_no].x - rec_h_acceptor[atm_rec_no].x ,2) +
				pow( lig_h_donor[atm_lig_no].y - rec_h_acceptor[atm_rec_no].y ,2) +
				pow( lig_h_donor[atm_lig_no].z - rec_h_acceptor[atm_rec_no].z ,2)
			);
			
			if( distance <= distance_cutoff ){
				printf("LIG-%s\tdonates_to\t%s-%d %s\t%.2f\t%s %s\n", 
					lig_h_donor[atm_lig_no].atm_name, 
					rec_h_acceptor[atm_rec_no].res_name, 
					rec_h_acceptor[atm_rec_no].res_number,
					rec_h_acceptor[atm_rec_no].atm_name,
					distance,
					rec_filename, 
					lig_filename
				);
			}	
		}		
	}
	
	
	
	// Find the ligand-receptor h-bonds in which the ligand is accepting H
	for( atm_lig_no=1 ; atm_lig_no<=total_lig_h_acceptor ; atm_lig_no++ ){
		for( atm_rec_no=1 ; atm_rec_no<=total_rec_h_donor; atm_rec_no++ ){
			
			distance=sqrt( 
				pow( lig_h_acceptor[atm_lig_no].x - rec_h_donor[atm_rec_no].x ,2) +
				pow( lig_h_acceptor[atm_lig_no].y - rec_h_donor[atm_rec_no].y ,2) +
				pow( lig_h_acceptor[atm_lig_no].z - rec_h_donor[atm_rec_no].z ,2)
			);
			
			if( distance <= distance_cutoff ){
				printf("LIG-%s\taccepts_from\t%s-%d\t%s\t%.2f\t%s %s\n", 
					lig_h_acceptor[atm_lig_no].atm_name, 
					rec_h_donor[atm_rec_no].res_name, 
					rec_h_donor[atm_rec_no].res_number, 
					rec_h_donor[atm_rec_no].atm_name, 
					distance,
					rec_filename, 
					lig_filename
				);
			}	
		}		
	}
	
	
	
return 0;
}