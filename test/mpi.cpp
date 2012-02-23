#include <iostream>
#include <vector>
#include "mpi.hpp"

int main(int argc, char *argv[]) {
   
    MPIi::Init(argc, argv);
    unsigned int id = MPI::COMM_WORLD.Get_rank();
    unsigned int size = MPI::COMM_WORLD.Get_size();

    if(id == 0)
    {
        std::cout << id << ": We have " << size << " processors\n";
        for(i = 1;i < size; ++i)
        {
            std::cout << "Calling processor " << i << "\n";
            std::string s = "My message";
            MPI::COMM_WORLD.Send( &s, s.size(), MPI::CHAR, i, 1);
        }
        for(i = 1;i < size; ++i)
        {
            std::string s;
            MPI::COMM_WORLD.Recv(&s, 100, MPI::CHAR, MPI::ANY_SOURCE, MPI::ANY_TAG);
            std::cout << id << ": Message back from " << i << " : " << s << "\n";
        }
    }
    else
    {
        /* receive from rank 0: */
        MPI::COMM_WORLD.Recv(&s, 100, MPI::CHAR, 0, 1);
        std::string s = "Processor "; s += id; s+= " reporting for duty\n";
        /* send to rank 0: */
        MPI_Send(s, s.size(), MPI::CHAR, 0, 1);
    }

    MPI::Finalize(); 

    return 0;
}
