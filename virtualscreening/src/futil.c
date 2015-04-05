#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>

#include "messages.h"
#include "defines.h"
#include "futil.h"
#include "string_owner.h"

owner_file_t *allocate_owner_file(const int *size){
    int i;
    owner_file_t * aux = NULL;
    aux = (owner_file_t*)malloc(sizeof(owner_file_t)* *size);
    for (i=0; i < *size; i++){
        aux[i].file_name = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
    }
    return aux;
}

void desallocate_owner_file(owner_file_t *aux, const int *size){
    int i;
    for (i=0; i < *size; i++){
        free(aux[i].file_name);
    }
    free(aux);
}

void show_owner_file(const owner_file_t *aux, const int *size){
    int i;
    for (i=0; i < *size; i++){
        printf("%s\n", aux[i].file_name);        
    }    
}


char * path_join_file(const char *path, const char *f){
	char *path_file;
	int len;
	len = (strlen(path)+ strlen(f)+1);
	if (len > MAX_PATH_FILE_NAME){
		fatal_error("In path_join_file function, len variable is more than MAX_PATH_FILE_NAME. Please check it. \n");
	}
	path_file = (char*)malloc(len);
    strcpy(path_file,path);
    strcat(path_file,f);
    return path_file;
}

FILE * open_file(__const char * filename, mode_files_t mode){
	FILE *aux_file=NULL;
	if (mode == fREAD){
		aux_file = fopen(filename, "r+");
		if (!aux_file){
			printf("ERROR: Can not open File %s\n ",filename);
			perror("Error when trying open file \n");
		}
	}else if ((mode == fWRITE)){
		aux_file = fopen(filename, "w+");
	}else if ((mode == fAPPEND)){
		aux_file = fopen(filename, "a+");
	}
	return aux_file;
}

boolean_t file_is_empty(FILE *file_aux){
	/* Checks if file is empty.
	 * Returns true when file is empty. Otherwise, false.
	 */
	//char line_aux[MAX_LINE_FILE];
	if ( feof(file_aux)){		
		return btrue;
	}else{
		return bfalse;
	}
}

/** Returns true if exists a file. Otherwise, returns false.
*/
boolean_t check_exists_file(const char *path_file_name){	
    struct stat buffer;
    int exist = stat(path_file_name,&buffer);
    if(exist == 0){
    	return btrue;
    }        
	return bfalse;
}

/** Get the size of FILE
*/
long get_file_size(FILE *fp){
    long fsize = 0;

    fseek(fp,0,SEEK_END);
    fsize = ftell(fp); 
    fseek(fp,0,SEEK_SET);

    return fsize;
}

/** Get the last line of file
*/
char *get_last_line(const char *filepath){
    FILE *fp;
    char buff[4096+1];
    int size,i;
    long fsize;
    if(NULL==(fp=fopen(filepath, "r"))){
        perror("file cannot open at lastline");
        return NULL;
    }
    fsize= -1L*get_file_size(fp);
    if(size=fseek(fp, MAX(fsize, -4096L), SEEK_END)){
        perror("cannot seek");
        exit(1);
    }
    size=fread(buff, sizeof(char), 4096, fp);
    fclose(fp);
    buff[size] = '\0';
    i=size-1;
    if(buff[i]=='\n'){
        buff[i] = '\0';
    }
    while(i >=0 && buff[i] != '\n')
        --i;
    ++i;
    return strdup(&buff[i]);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

/* How Many files are at directory by extension
*/
int how_many_files_directory_by_extension(const char *path, const char *ext){   
    const char *ext_aux = NULL; 
    struct dirent *ent = NULL;
    int r = 0;

    DIR *dir = opendir(path);
    while (ent = readdir(dir)) {
        ext_aux = get_filename_ext(ent->d_name);
        if (strcmp (ext_aux, ext) == 0){
            r++;
        }   
    }
    closedir(dir);
    return r;
}

void insert_files_directory_by_extension(owner_file_t *file_names, const char *path, const char *ext){
    const char *ext_aux = NULL;
    int index = -1;

    DIR *dir = opendir(path);
    struct dirent *ent;
    while (ent = readdir(dir)) {
        ext_aux = get_filename_ext(ent->d_name);
        if (strcmp (ext_aux, ext) == 0){
            index = index +1;
            strcpy(file_names[index].file_name, ent->d_name);
        }   
    }
    closedir(dir);
}

//Set file name without extension
void set_base_file_name(char *filename_no_ext, const char *filename_source, const char *ext){
    char aux_ext[MAX_SUFIX];
    char *aux=NULL;
    char *ch=NULL;
    aux = (char*)malloc(sizeof(char)*MAX_FILE_NAME);
    strcpy(aux_ext, ".");
    strcat(aux_ext, ext);
    strcpy(aux, filename_source);
    ch=strrchr(aux,'.');
    substring(filename_no_ext, filename_source, 0, strlen(filename_source)-strlen(ch));
    free(aux);
}

//Add extension in file name
void add_ext_file_name(char *filename, const char *ext){
    strcat(filename, ".");
    strcat(filename, ext);
}