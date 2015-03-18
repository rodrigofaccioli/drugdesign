#ifndef OLD_MPI_OWNER_H
#define OLD_MPI_OWNER_H

void send_broadcast(void* data, int count, MPI_Datatype datatype, int root,
              MPI_Comm communicator);

#endif