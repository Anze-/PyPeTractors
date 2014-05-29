#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "@function@"

#include <fstream>
using namespace std;


int main(int argc, char *argv[])
{
    int numprocessors, rank, namelen, i,limit,slaves;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
i=0;
  p_drop();
    if ( rank == 0 )
    {
        std::cout << "Master core: " << processor_name << " : " << rank <<"\n";
// std::cout << "master (" << rank << "/" << numprocessors << ")\n";
    } else {
// std::cout << "slave (" << rank << "/" << numprocessors << ")\n";


slaves=(numprocessors-1);
while(i<(img_s/slaves)){
int slice=rank+i*slaves;
doTheMath(slice);
//char c[10];
//sprintf(c, "%d", slice);
//ofstream myfile;
//myfile.open(c);
//riga=my_ext_function(slice)
//myfile << riga;
//myfile.close();
i++;
   }
   }
   if(rank <= (img_s%slaves))
   {
    int slice=rank+i*slaves;
    doTheMath(slice);
   }
   MPI_Finalize();
    if ( rank == 0)
    {
        std::cout << "Master core: " << processor_name << " : " << rank << "@function@ : END" <<"\n";
// printf("%d",i);
// ofstream myfile;
// myfile.open("200");
// myfile << 15;
// myfile.close();
//char output[200000] = system("ls");
    }

   return 0;
}
