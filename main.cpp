#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include "network.h"


#define PSEC 0
#define	PCO 1

using namespace std;


void psecOrpcoVslambda(Simulator s, Parameters pm, int myid, int numprocs, int loopNum, int flag)
{
   int cnt = 21;
   double myres[21];
   double delta_lambdaP = 0.00025;
   int myloopNum = loopNum/numprocs + 1;
   s.setLoopNum(myloopNum);

   for(int n = 0; n < cnt; n++)
   {
      pm.den_pj = n * delta_lambdaP;
      s.initParameters(pm);
      (flag == PSEC)? myres[n] = s.secrecyProbabiltiy(): myres[n] = s.connectionProbability();
   }

    if(0 != myid)
    {
     	MPI_Send(myres, cnt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
     	double resfromOthers[numprocs - 1][cnt];
        (flag == PSEC)?cout<<"lambdaP, psec(sim)"<<endl: cout<<"lambdaP, pco(sim)"<<endl;
        for(int i = 1; i < numprocs; i++)
        {
            MPI_Recv(resfromOthers[i-1], cnt, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int n = 0; n < cnt; n++)
            {
               myres[n] += resfromOthers[i-1][n];
               if (i == numprocs - 1)
                   cout<<n * delta_lambdaP<< ", "<< myres[n]/numprocs<<endl;
            }
	}
    }

}

int main(int argc, char* argv[])
{

    double rt = 8, re = 4;
    double lambdaT = 0.00001, lambdaP = 0.002, lambdaE = 0.001;
    double q = 0.5;
    Parameters pmts = {
        /*network-related parameters*/
        500,                           ///network radius
        lambdaT,                       ///density of D2D pair
        lambdaP,                       ///density of potential jammers    
        lambdaE,                       ///density of eavesdroppers
        /*transmission-related parameters*/
        1.0,                           ///D2D transmit power
        pow(2,re)-1,                        ///eavesdropper decoding threshold
        pow(2,rt)-1,                         ///legitimate node decoding threshold
        /*LDCJ parameters*/
        q,                             ///jammer selection probability inside LoS balls
        /*antenna-related parameters*/
        PI/6,                          ///tx's beamwidth of main lobe
        PI/6,                          ///rx's beamwidth of main lobe
        2 * PI,                        ///eve's beamwidth of main lobe
        0.1,                           ///tx's back lobe gain
        0.1,                           ///rx's back lobe gain
        10,                            ///eve's back lobe gain
        10,                            ///tx's main lobe gain
        10,                            ///rx's main lobe gain
        10,                            ///eve's main lobe gain
    };
    Simulator s;
    int loopNum = 10000;
    /*s.setLoopNum(loopNum);
    s.initParameters(pmts);
    cout<<s.connectionProbability()<<endl;*/
    
    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if(0 == myid) cout<<"lambdaT: "<<lambdaT<<" lambdaE: "<<lambdaE<<" re: "<<re<<" rt: "<<rt<<" q: "<<q<<endl;
 /* double myres, myloopNum;  
    myloopNum = loopNum/numprocs + 1;
    s.initParameters(pmts);
    s.setLoopNum(myloopNum);
    myres = s.connectionProbability();
    myres = s.connectionProbability();
    myres = s.secrecyProbabiltiy();
    if(0 != myid)
    {
        MPI_Send(&myres, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        double resfromOthers[numprocs - 1];
        for(int i = 1; i < numprocs; i++)
        {
            MPI_Recv(&resfromOthers[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            myres += resfromOthers[i];
        }
        cout<<"lambdaT: "<<lambdaT<<" lambdaJ: "<<lambdaJ<<" lambdaE: "<<lambdaE<<" q: "<<q<<" re: "<<re<<endl;
        cout<<"sum secrecy probability: "<< myres/numprocs<<endl;
    }
*/

    psecOrpcoVslambda(s, pmts, myid, numprocs, loopNum, PCO);
    MPI_Finalize();

    return 0;
}
