#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <iostream>
#include <omp.h>

using namespace std;


int test()
{
    

}

int main()
{
    
#pragma omp parallel
    {
#pragma omp critical
        {
            printf("Hello OpenMP from thread %i\n", omp_get_thread_num());
        }
    }
    return 0;
}

