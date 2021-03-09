#include "itensor/all.h"
#include "stdlib.h"

#include "omp.h"
#include "thread"

#include "cmath"
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#define MAX_DIGITS 10


using namespace itensor;
using std::vector;
using namespace std::chrono;

int main(int argc, char* argv[])
{

    printfln("\nModel parameters of Hubbardteste.cc\n");

    //Dimensionality
    int L = 17;
    printfln("L = ", L);

    //Filling
    int Npart = L;
    printfln("Npart = ", Npart);

    //int L;
    //printf("L = ");
    //scanf("%d",&L);
    //Filling
    //int Npart;
    //printf("Npart = ");
    //scanf("%d",&Npart);

    //Parameters

    int t1 = 1; //hopping

    auto U0 = 0.0; //onsite inicial
    printfln("U0 = ", U0);

    auto  V0 = 0.0; //intersite inicial
    printfln("V0 = ",V0);

    auto Uf = 1.; //onsite final
    printfln("Uf = ", Uf);

    auto Vf = 2.0; //intersite final
    printfln("Vf = ",Vf);

    //Details of time evolution

    Real tstep = 0.2;
    Real ttotal = 5.0;
    Real cutoff = 1E-12;

    //Variável loop, determina sobre qual instante de tempo t
    //a evolução psi(t) irá começar.
    //É necessário saber onde o cálculo anterior parou.
    //Por exemplo, se parou no instante t=2, então set a variável abaixo
    //para tevol = 2.0.
    double tevol = 0.0;


    //Type model

    auto sites = Electron(L);

    auto state = InitState(sites);

    int p = Npart; //Filling

    for(int i = L; i >= 1; --i)
        {
        if(p > i)
            {
            println("Doubly occupying site ",i);
            state.set(i,"UpDn");
            p -= 2;
            }
        else
        if(p > 0)
            {
            println("Singly occupying site ",i);
            state.set(i,(i%2==1 ? "Up" : "Dn"));
            p -= 1;
            }
        else
            {
            state.set(i,"Emp");
            }
        }

    auto psi0 = MPS(state);

    ////////////////////////////////////////////
    //                                        //
    //Contructing Hamiltonian of Ground State //
    //                                        //
    ////////////////////////////////////////////

    auto ampo = AutoMPO(sites);

    for(int i = 1; i <= L; ++i)
        {
        ampo += U0,"Nupdn",i;
        }
    for(int b = 1; b < L; ++b)
        {
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;
        ampo += V0,"Ntot",b,"Ntot",b+1;
        }
    auto H = toMPO(ampo);

    ///////////////////////////////////////////
    //
    // Constructing Hamiltonian of final state
    //
    ///////////////////////////////////////////

    auto ampof = AutoMPO(sites);

    for(int i = 1; i <= L; ++i)
        {
        ampof += Uf,"Nupdn",i;
        }
    for(int b = 1; b < L; ++b)
        {
        ampof += -t1,"Cdagup",b,"Cup",b+1;
        ampof += -t1,"Cdagup",b+1,"Cup",b;
        ampof += -t1,"Cdagdn",b,"Cdn",b+1;
        ampof += -t1,"Cdagdn",b+1,"Cdn",b;
        ampof += Vf,"Ntot",b,"Ntot",b+1;
        }
    auto Hf = toMPO(ampof);

    /////////////////
    //             //
    //   DMRG  GS  //
    //             //
    /////////////////

    auto sweeps = Sweeps(30);
         sweeps.maxdim() = 50,100,200,400,600,600,800,800;
         sweeps.cutoff() = 1E-8;
    auto [energyGS,psi] = dmrg(H,psi0,sweeps,{"Quiet",true}); //output psi updated

    auto psiGS = psi; //Salve the Ground State in psiGS


    /////////////////////////////////
    //                             //
    // Save ground state and sites //
    //                             //
    /////////////////////////////////

    writeToFile("sites_t_0.000000",sites);
    writeToFile<MPS>("psi_t_0.000000",psi); //Writes ground-state into "psi_t_0.000000"


    ///////////////////////////
    //                       //
    // DMRG final state (FS) //
    //                       //
    ///////////////////////////

    auto [energyFS,psiFS] = dmrg(Hf,psi0,sweeps,{"Quiet",true}); //output psi updated

   //Parametro m_cdw do GS e do FS

    Vector upd0(L),dnd0(L); //Para GS representado por psiGS
    for(int j = 1; j <= L; ++j)
        {
        psiGS.position(j);
        upd0(j-1) = eltC(dag(prime(psiGS(j),"Site"))*op(sites,"Nup",j)*psiGS(j)).real();
        dnd0(j-1) = eltC(dag(prime(psiGS(j),"Site"))*op(sites,"Ndn",j)*psiGS(j)).real();
        }

   Vector updf(L),dndf(L); //Para FS representado por psif

    for(int j = 1; j <= L; ++j)
       {
       psiFS.position(j);
       updf(j-1) = eltC(dag(prime(psiFS(j),"Site"))*op(sites,"Nup",j)*psiFS(j)).real();
       dndf(j-1) = eltC(dag(prime(psiFS(j),"Site"))*op(sites,"Ndn",j)*psiFS(j)).real();
       }

    auto mcdwGS = 0.; //Para contar o GS

    for(int j = 0; j < L; ++j)
    {
        mcdwGS += pow(-1,j+1)*((upd0(j)+dnd0(j))-1);
    }

    auto mcdwFS = 0.; //Para contar o FS

    for(int j = 0; j < L; ++j)
    {
        mcdwFS += pow(-1,j+1)*((updf(j)+dndf(j))-1);
    }


    auto mocGS = sqrt(pow((mcdwGS/L),2)); //modulo de mcdw de GS
    auto mocFS = sqrt(pow((mcdwFS/L),2)); //module de mcdw de FS


    printfln("\n|mcdw_GS| = ", mocGS);
    printfln("|mcdw_FS| = ", mocFS);

    printfln("\nEnergy GS = ", energyGS);
    printfln("Energy FS = ", energyFS);

    printfln("|EnergyFS| - |EnergyGS| = ", sqrt(pow(energyFS,2)) - sqrt(pow(energyGS,2)));

    
    ////////////////////////////////////////////////
    //                                            //
    // Rotina que salvar psi à medida que o tempo //
    // evolui e avalia o paramentro de ordem CDW  //
    //                                            //
    ////////////////////////////////////////////////

    char num_char[MAX_DIGITS + sizeof(char)];

    Real ttstep = 1.0; //Intervalo de tempo em que os psi(t) serão salvos. Não confudir tstep.

    //A variável anterior menor o ttotal, pois não queremos ultrapassar
    //o limite de tempo total desejado.
    double ttime = ttotal - tevol;

    //Variável do for, para determinar quantos loops serão executados
    auto tloop = (ttime/ttstep);

    for(auto time = 0.; time <= tloop; time += ttstep)
    {
        if(tevol < (ttotal))
        {
            printfln("\ntevol = %f",tevol);

            readFromFile(format("sites_t_%f",tevol),sites);
            auto psiEvol = readFromFile<MPS>(format("psi_t_%f",tevol));

            //Parâmetro de ordem
            auto mc = 0.;

            for(int j = 1; j <= L ; j++)
            {
                psiEvol.position(j);

                auto ket = psiEvol(j);
                auto bra = dag(prime(ket,"Site"));

                auto N_op = op(sites,"Ntot",j);

                auto Ntott = eltC(bra*N_op*ket).real();

                mc += pow(-1,j)*(Ntott - 1.0);

            }

            printfln("mc = ", mc);


            auto mcc = sqrt(pow((mc/L),2)); //module CDW parameter

            printfln("## %d %.12f",(tevol), mcc); //prints module CDW_parameter in time
            printfln("#1 %d %.12f",(tevol), mocGS); //prints module CDW_parameterGS in time
            printfln("#* %d %.12f",(tevol), mocFS); //prints module CDW_parameterFS in time
            printfln("#0 %d %.12f",(tevol), mc/L); //prints CDW_parameter in time
            printfln("#2 %d %.12f",(tevol), mcdwFS/L); //prints CDW_parameter in time
            printfln("#3 %d %.12f",(tevol), mcdwGS/L); //prints CDW_parameter in time


            // construindo os gates da evolução
            auto gates = vector<BondGate>();

            for(int b = 1; b < L; ++b)
            {
                auto hterm = -t1*sites.op("Adagup*F",b)*sites.op("Aup",b+1);
                     hterm += -t1*sites.op("Adagdn",b)*sites.op("F*Adn",b+1);
                     hterm +=  t1*sites.op("Aup*F",b)*sites.op("Adagup",b+1);
                     hterm +=  t1*sites.op("Adn",b)*sites.op("F*Adagdn",b+1);

                 hterm += Vf*sites.op("Ntot",b)*sites.op("Ntot", b+1);

                 hterm += Uf*sites.op("Nupdn",b)*sites.op("Id",b+1);

                auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                gates.push_back(g);
            }

            for(int b = L-1; b <= L-1; ++b)
            {
                 auto hterm = Uf*sites.op("Id",b)*sites.op("Nupdn", b+1);

                auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                     gates.push_back(g);
            }

            for(int b = L-1 ; b >= 1; --b)
            {
                auto hterm = -t1*sites.op("Adagup*F",b)*sites.op("Aup",b+1);
                     hterm += -t1*sites.op("Adagdn",b)*sites.op("F*Adn",b+1);
                     hterm +=  t1*sites.op("Aup*F",b)*sites.op("Adagup",b+1);
                     hterm +=  t1*sites.op("Adn",b)*sites.op("F*Adagdn",b+1);

                     hterm += Vf*sites.op("Ntot",b)*sites.op("Ntot", b+1);

                     hterm += Uf*sites.op("Nupdn",b)*sites.op("Id",b+1);

                auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                     gates.push_back(g);
            }

            for(int b = L-1; b >= L-1; --b)
            {
                 auto hterm = Uf*sites.op("Id",b)*sites.op("Nupdn", b+1);

                 auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                     gates.push_back(g);
            }
            //Finish the bulding of gates

            gateTEvol(gates,ttstep,tstep,psiEvol,{"Cutoff=",cutoff,"Verbose=",true});

            tevol += ttstep;

            printfln("Final do lopp, tevol = ", tevol);

            writeToFile(format("sites_t_%f",tevol),sites);
            writeToFile<MPS>(format("psi_t_%f",tevol),psiEvol); //Writes psi evol by gateTEvol

        }
        else
        {
            printfln("\ntevol = %f",tevol);

            readFromFile(format("sites_t_%f",tevol),sites);
            auto psiEvol = readFromFile<MPS>(format("psi_t_%f",tevol));

            //Parâmetro de ordem
            auto mc = 0.;

            for(int j = 1; j <= L ; j++)
            {
                psiEvol.position(j);

                auto ket = psiEvol(j);
                auto bra = dag(prime(ket,"Site"));

                auto N_op = op(sites,"Ntot",j);

                auto Ntott = eltC(bra*N_op*ket).real();

                mc += pow(-1,j)*(Ntott - 1.0);

            }

            printfln("mc = ", mc);


            auto mcc = sqrt(pow((mc/L),2)); //module CDW parameter

            printfln("## %d %.12f",(tevol), mcc); //prints module CDW_parameter in time
            printfln("#1 %d %.12f",(tevol), mocGS); //prints module CDW_parameterGS in time
            printfln("#* %d %.12f",(tevol), mocFS); //prints module CDW_parameterFS in time
            printfln("#0 %d %.12f",(tevol), mc/L); //prints CDW_parameter in time
            printfln("#2 %d %.12f",(tevol), mcdwFS/L); //prints CDW_parameter in time
            printfln("#3 %d %.12f",(tevol), mcdwGS/L); //prints CDW_parameter in time

        }
    }


return 0;
}

