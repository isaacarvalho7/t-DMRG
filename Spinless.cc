#include "itensor/all.h"

using namespace itensor;
using std::vector;

int main()
    {
    
    //Model parameter

    int N = 50;
    printfln("N = ",N);
               
    Real tstep = 0.1;           
    Real ttotal = 1.0;           
    Real cutoff = 1E-10;           

    //Quench parameter
    
    auto th = 5.0;
    auto V0=0.5;
    auto Vf=5.;


    auto sites = Fermion(N);

    auto state = InitState(sites);

    for(auto j=1; j <=N ; j++)
        {
        if(j%2==1)
            {
            state.set(j,"Occ");
            }
        else
            {
            state.set(j,"Emp");
            }
        }

    auto psi0 = MPS(state);

    auto ampo = AutoMPO(sites);

    //Contructed Hamiltonian

    for(int j = 1; j < N; ++j)
        {
        ampo += -th,"Adag",j,"A",j+1;
        ampo += -th,"Adag",j+1,"A",j;

        ampo += V0,"N",j,"N",j+1;

        }

    auto H = toMPO(ampo);

    //Constructed DMRG
    auto sweeps = Sweeps(5); //number of sweeps is 5
         sweeps.maxdim() = 10,20,100,100,200; //gradually increase states kept
         sweeps.cutoff() = 1E-10; //desired truncation error

    //auto [energy,psi] = dmrg(H,psi,sweeps);
    auto [energy,psi] = dmrg(H,psi0,sweeps); //DMRG

    println("Ground State Energy = ",energy); 


    //occupation density of GroundState

    println("\nj N(j) = ");
    for( auto j : range1(N) ) 
        {
        //re-gauge psi to get ready to measure at position j
        psi.position(j);
        auto ket = psi(j);
        auto bra = dag(prime(ket,"Site"));
        auto Njop = op(sites,"N",j);
        //take an inner product 
        auto nj = elt(bra*Njop*ket);
        printfln("%d %.12f",j,nj);
        }

    //Create a std::vector (dynamically sizeable array)
    //to hold the Trotter gates

    auto gates = vector<BondGate>();

    //Create the gates exp(-i*tstep/2*hterm)
    //and add them to gates

    for(int b = 1; b <= N-1; ++b)
    {
    auto hterm = -th*op(sites,"Adag",b)*op(sites,"A",b+1);
    hterm += -th*op(sites,"Adag",b+1)*op(sites,"A",b);

    hterm += Vf*op(sites,"N",b)*op(sites,"N",b+1);

    auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
    gates.push_back(g);
    }

    //Create the gates exp(-i*tstep/2*hterm) in reverse
    //order (to get a second order Trotter breakup which
    //does a time step of "tstep") and add them to gates

    for(int b = N-1; b >= 1; --b)
    {
    auto hterm = -th*op(sites,"Adag",b)*op(sites,"A",b+1);
    hterm += -th*op(sites,"Adag",b+1)*op(sites,"A",b);

    hterm += Vf*op(sites,"N",b)*op(sites,"N",b+1);

    auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
    gates.push_back(g);
    }

    //Save initial state GroundState (GS);
    auto psiGS = psi;

    //Time evolve, overwriting psi when done
    gateTEvol(gates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});

    //Occupation density GS evol in time
    
    println("\nj N(j) = ");
    for( auto j : range1(N) )
        {
        //re-gauge psi to get ready to measure at position j
        psi.position(j);
        auto ket = psi(j);
        auto bra = dag(prime(ket,"Site"));
        auto Njop = op(sites,"N",j);
        //take an inner product 
        //
        auto nj = eltC(bra*Njop*ket);
        //printfln("%d %.12f",j,nj); //output is a real number
        
        auto nreal = nj.real();
        printfln("%d %.12f",j,nreal);

       
        }
        


return 0;

}
