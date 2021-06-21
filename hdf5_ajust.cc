#include "itensor/all.h"
#include "stdlib.h"
#include "omp.h"
#include "thread"
#include "cmath"
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#include "H5Cpp.h"

#define MAX_DIGITS 10


using namespace itensor;
using std::vector;
using namespace std::chrono;
using namespace std;

int main(int argc, char* argv[])
{

    //Args::global().add("UseSVD=",true);
    //Args::global().add("SVDMethod=", "gesdd");

    printfln("\nModel parameters of Hubbardteste.cc\n");


    int L =5;//Dimensionality
    printfln("L = ", L);

    int Npart = L;//Filling
    printfln("Npart = ", Npart);

    auto t1 = 1.0; //hopping
    printfln("t_hopp = ", t1);

    auto U0 = 0.; //onsite inicial
    printfln("U_0 = ", U0);

    auto  V0 = 0.0; //intersite inicial
    printfln("V_0 = ",V0);

    auto Uf = 1.; //onsite final
    printfln("U_f = ", Uf);

    auto Vf = 2.0; //intersite final
    printfln("V_f = ",Vf);

    //pinning fields (Placed at the ends of the chain )
    auto h_z = -0.01;
    printfln("pinning field = ", h_z);

    Real tstep = 0.05; //accuracy of time evolution
    Real cutoff = 1E-12; // error dicarded
    int sweep_n = 1; //number'sweep
    printfln("Sweep_number = ", sweep_n);

    ///Time evolution after finite time quench

    auto time_evolution = 0.5;
    printfln("time evolution after quench = ", time_evolution);


    auto tau_U = 0.01;

    printfln("tau__U ",tau_U);

    auto tau_V = (tau_U*(sqrt(pow(Uf-U0,2))))/(sqrt(pow(Vf-V0,2)));

    printfln("tau__V ",tau_V);






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



    //defining sgn function
    auto sgn_U = Uf-U0;
    auto sgn_V = Vf-V0;


    if(sgn_U < 0.){sgn_U = -1.;}
    if(sgn_U == 0.){sgn_U = 0.;}
    if(sgn_U > 0){sgn_U = 1.;}


    if(sgn_V < 0.){sgn_V = -1.;}
    if(sgn_V == 0.){sgn_V = 0.;}
    if(sgn_V > 0){sgn_V = 1.;}

    double time_f = (sqrt(pow(Uf - U0,2))*tau_U);


    //Condition so that we have at least 20 (p_t = 20 below) points on the graph
    //Increases the tstep precision so that this occurs when tau_U is too small.

    int p_t = 1;
    if(time_f <= 0.09){tstep = time_f/p_t;}
    if(time_f > 0.09 && time_f < 0.2){tstep = time_f/p_t;}


    int n_loop = int(time_f/tstep); //number of steps of evolution


    Real del_t = tstep;

    Vector U_t(n_loop + 1) , V_t(n_loop + 1) ;

    for(int j = 0; j <= n_loop; j += 1)
    {
        double d_t = j*del_t;

        U_t(j)= U0 + (sgn_U*d_t)/tau_U;
        cout << "U_t(" << j << ") = " << U_t(j) << "\n";
    }
    for(int j = 0; j <= n_loop; j += 1)
    {
        double d_t = j*del_t;

        V_t(j)= V0 + (sgn_V*d_t)/tau_V;
        cout << "V_t(" << j << ") = " << V_t(j) << "\n";
    }


    //Ground State (GS) Hamiltonian
    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= L; ++i)
        {
        ampo += U_t(0),"Nupdn",i;
        }
    for(int b = 1; b < L; ++b)
        {
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;
        ampo += V_t(0),"Ntot",b,"Ntot",b+1;
        }
        ampo += h_z, "Sz", 1;
        ampo += h_z, "Sz", L;

    auto H0 = toMPO(ampo);
    //DRMG routine
    auto sweeps = Sweeps(sweep_n);
         sweeps.maxdim() = 10,10,20,20,40,40,80,80,160,160,320,320,400,400;
         sweeps.cutoff() = 1E-10;
    auto [energyGS,psi] = dmrg(H0,psi0,sweeps,{"Quiet",true}); //output psi updated

    auto psi_ini = psi; //Salve the Ground State in psiGS
    auto psi_Evol = psi; // Defines psi_Evol that I will use to realize the finite time quench


    //////////////////////
    ///                 //
    /// Very attention! //
    ///                 //
    //////////////////////

    //I set this variable to
    //if a_measure=1 measures light cones during scanning
    //if a_measure=0 does not measure light cones during scanning.
    //Very attention!

    int a_measure_cones_varredura = 0; //todos os cones =1 mede ; =0 não mede


    int a_measure_emaranhamento= 0;
    int a_measure_den_carga= 0;
    int a_measure_den_sz= 0;
    int a_measure_corr_carga= 0;
    int a_measure_corr_sz= 0;
    int a_measure_ener_quench= 0;



	//////////////////////////////////
	//								//
	//	Finite time quantum quench  //
	//								//
	///////////////////////////////////

    for(int lt = 0; lt <= n_loop; lt ++)
    {

        //Ground State (GS) Hamiltonian
        auto ampo = AutoMPO(sites);
        for(int i = 1; i <= L; ++i)
            {
            ampo += U_t(lt),"Nupdn",i;
            }
        for(int b = 1; b < L; ++b)
            {
            ampo += -t1,"Cdagup",b,"Cup",b+1;
            ampo += -t1,"Cdagup",b+1,"Cup",b;
            ampo += -t1,"Cdagdn",b,"Cdn",b+1;
            ampo += -t1,"Cdagdn",b+1,"Cdn",b;
            ampo += V_t(lt),"Ntot",b,"Ntot",b+1;
            }
            ampo += h_z, "Sz", 1;
            ampo += h_z, "Sz", L;

        auto H = toMPO(ampo);
        //DRMG routine
        auto sweeps = Sweeps(sweep_n);
             sweeps.maxdim() = 10,10,20,20,40,40,80,80,160,160,320,320,400,400;
             sweeps.cutoff() = 1E-10;

        steady_clock::time_point t0 = steady_clock::now(); //time from here

        auto [energyGS,psi_GS] = dmrg(H,psi0,sweeps,{"Quiet",true}); //output psi updated

        steady_clock::time_point t2 = steady_clock::now(); //final calculation time
        duration<double> time_span = duration_cast<duration<double>>(t2 - t0);
        cout << "#time_drmg  " << sqrt(pow(U_t(lt)-U0,2)) << " " << time_span.count() <<"\n";
        cout << "#time_drmg_t  " << lt*del_t << " " << time_span.count() <<"\n";


        ///Measures the variance or energy error GS
        auto C = nmultMPO(prime(H),H,{"MaxDim",5000,"Cutoff",1E-14});
        auto h_0 = inner(psi_GS,H,psi_GS);
        auto h_2 = inner(psi_GS,C,psi_GS);
        auto variance = h_2 - pow(h_0,2);

        printfln("<psi_GS,H,psi_GS>");
        cout << "#h_0 " << sqrt(pow(U_t(lt)-U0,2)) << " " <<  h_0 <<"\n";
        cout << "#h_0_t " << lt*del_t << " " << h_0 <<"\n";

        printfln("<psi_GS,H,H,psi_GS>");
        cout << "#h_2 " << sqrt(pow(U_t(lt)-U0,2)) << " " <<  h_2 <<"\n";
        cout << "#h_2_t " << lt*del_t << " " << h_2 <<"\n";


        printfln("Variance");
        cout << "#variance " << sqrt(pow(U_t(lt)-U0,2)) << " " <<  variance <<"\n";
        cout << "#variance_t " << lt*del_t << " " << variance <<"\n";

        printfln("Ground state energy");
        cout << "%energy_GS " << sqrt(pow(U_t(lt)-U0,2)) << " " <<  energyGS <<"\n";
        cout << "%energy_GS_t " << lt*del_t << " " << energyGS <<"\n";



        //Building gates of evolution
        auto gates = vector<BondGate>();
        for(int b = 1; b < L; ++b)
        {
            auto hterm = -t1*sites.op("Adagup*F",b)*sites.op("Aup",b+1);
                 hterm += -t1*sites.op("Adagdn",b)*sites.op("F*Adn",b+1);
                 hterm +=  t1*sites.op("Aup*F",b)*sites.op("Adagup",b+1);
                 hterm +=  t1*sites.op("Adn",b)*sites.op("F*Adagdn",b+1);
             hterm += V_t(lt)*sites.op("Ntot",b)*sites.op("Ntot", b+1);
             hterm += U_t(lt)*sites.op("Nupdn",b)*sites.op("Id",b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
            gates.push_back(g);
        }
        for(int b = L-1; b <= L-1; ++b)
        {
             auto hterm = U_t(lt)*sites.op("Id",b)*sites.op("Nupdn", b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }
        ///modifiquei daqui * aplicação campo h_z nas extremidades
        for(int b = L-1; b <= L-1; ++b)
        {
             auto hterm = h_z*sites.op("Id",b)*sites.op("Sz", b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }

        for(int b = 1; b <= 1; ++b)
        {
             auto hterm = h_z*sites.op("Sz",b)*sites.op("Id", b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }
        ///modifiquei até aqui * aplicação campo h_z nas extremidades
        for(int b = L-1 ; b >= 1; --b)
        {
            auto hterm = -t1*sites.op("Adagup*F",b)*sites.op("Aup",b+1);
                 hterm += -t1*sites.op("Adagdn",b)*sites.op("F*Adn",b+1);
                 hterm +=  t1*sites.op("Aup*F",b)*sites.op("Adagup",b+1);
                 hterm +=  t1*sites.op("Adn",b)*sites.op("F*Adagdn",b+1);
                 hterm += V_t(lt)*sites.op("Ntot",b)*sites.op("Ntot", b+1);
                 hterm += U_t(lt)*sites.op("Nupdn",b)*sites.op("Id",b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }
        for(int b = L-1; b >= L-1; --b)
        {
             auto hterm = U_t(lt)*sites.op("Id",b)*sites.op("Nupdn", b+1);
             auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }
        ///modifiquei daqui * aplicação de campo externo nas extremidades
        for(int b = L-1; b >= L-1; --b)
        {
            auto hterm = h_z*sites.op("Id",b)*sites.op("Sz", b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }

        for(int b = 1; b >= 1; --b)
        {
            auto hterm = h_z*sites.op("Sz",b)*sites.op("Id", b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
                 gates.push_back(g);
        }
        ///modifiquei até aqui * aplicação de campo externo nas extremidades
        ///Finish the bulding of gates


        ////////////////////////////////////////
        //                                    //
        //                                    //
        // GATE DE EVOLUÇÃO ESTÁ LOGO ABAIXO  //
        // GATE DE EVOLUÇÃO ESTÁ LOGO ABAIXO  //
        // GATE DE EVOLUÇÃO ESTÁ LOGO ABAIXO  //
        // GATE DE EVOLUÇÃO ESTÁ LOGO ABAIXO  //
        //                                    //
        //                                    //
        ////////////////////////////////////////


        gateTEvol(gates,tstep,tstep,psi_Evol,{"Quiet",true,"Cutoff=",cutoff,"Verbose=",true,"UseSVD=",true,"SVDMethod=","gesdd"});

		//Print o instante de tempo
		printfln("\n\nt = ", lt*del_t);

        printfln("\n\ntime do loop = ", lt*del_t);

        //Para entender como U e V estão variando no tempo
        cout << "#!U " << lt*del_t << " " << U_t(lt) <<"\n";
        cout << "#!V " << lt*del_t << " " << V_t(lt) <<"\n";


        ///to measure or not to measure light cones?
        ///If a_measure=1, it is measuring
        ///if a_measre=0, it is not measuring
        if(a_measure_cones_varredura==1) //open  if(a_measure==1)
        {

        if(a_measure_emaranhamento==1) //open a_measure_emaranhamento==1
        {
        //Para começar a medir no centro da cadeia
        printfln("Cone-de-luz entropia");

        //Reconstrui as quantidade abaixo para medir
        //as quantidade ao longo de toda a cadeia

        int i_med = 1; //contador p/ inicio da medida abaixo

        //if(L%2==1){i_med = (L-1)/2;} //cadeia par
        //else{i_med = L/2;} //cadeia ímpar

        for(int s_A = i_med; s_A < L ; s_A++) // < L, pois precisamos subsistema A e B
        {


            //Entropia GS
            auto b = s_A;
            printfln("Subsistema A de tamanho = ", b);

            psi_GS.position(b);
            auto l = leftLinkIndex(psi_GS,b); //subsistema tamanho b
            auto s = siteIndex(psi_GS,b);
            auto [U,S,V] = svd(psi_GS(b),{l,s});
            auto u = commonIndex(U,S);
            Real SvN = 0.;
            for(auto n : range1(dim(u)))
                {
                auto Sn = elt(S,n,n);
                auto p = sqr(Sn);
                if(p > 1E-12) SvN += -p*log(p);
                }
            printfln("Entropia psi_GS");

            cout << "#svn0_cone " << b << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << SvN <<"\n";
            cout << "#svn0_conet " << b << " " << lt*del_t << " " << SvN <<"\n";

            //Entropia psi_Evol
            auto b1 = s_A;

            psi_Evol.position(b1);
            auto l1 = leftLinkIndex(psi_Evol,b1);
            auto s1 = siteIndex(psi_Evol,b1);
            auto [U1,S1,V1] = svd(psi_Evol(b1),{l1,s1});
            auto u1 = commonIndex(U1,S1);
            Real SvN1 = 0.;
            for(auto n : range1(dim(u1)))
                {
                auto Sn = elt(S1,n,n);
                auto p = sqr(Sn);
                if(p > 1E-12) SvN1 += -p*log(p);
                }
            printfln("Entropia psi_Evol");
            cout << "#svn1 " << b1 << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << SvN1 <<"\n";
            cout << "#svntime_1 "<< b1 << " " << lt*del_t << " " << SvN1 <<"\n";
            //Diferença módulo entropia de psi_Evol e psi_GS
            printfln("Diferença entropia psi_Evol e psi_GS");

            cout << "#dife_svn_cone " << b1 << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(SvN1-SvN,2)) <<"\n";
            cout << "#dife_svn_conet " << b1 << " " << lt*del_t << " " << sqrt(pow(SvN1-SvN,2)) <<"\n";


            //Emaranhamento estado inicial psi_ini
            auto b2 = s_A;

            psi_ini.position(b2);
            auto l2 = leftLinkIndex(psi_ini,b2);
            auto s2 = siteIndex(psi_ini,b2);
            auto [U2,S2,V2] = svd(psi_ini(b2),{l2,s2});
            auto u2 = commonIndex(U2,S2);
            Real SvN2 = 0.;
            for(auto n : range1(dim(u2)))
                {
                auto Sn = elt(S2,n,n);
                auto p = sqr(Sn);
                if(p > 1E-12) SvN2 += -p*log(p);
                }
            printfln("Entropia psi_inicial");
            cout << "#svn2 " << b2 << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << SvN2 <<"\n";
            cout << "#svntime_2 " << b2 << " " << lt*del_t << " " << SvN2 <<"\n";
            //Diferença módulo entropia de psi_Evol e psi_inicial
            printfln("Diferença entropia psi_Evol e psi_inicial");

            cout << "#dife_svn2_cone " << b2 << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(SvN1-SvN2,2)) <<"\n";
            cout << "#dife_svn2_conet " << b2 << " " << lt*del_t << " " << sqrt(pow(SvN1-SvN2,2)) <<"\n";

        } //fim loop da entropia
        }//close if(a_measure_emaranhamento==1)

        if(a_measure_den_carga==1) //open a_measure_emaranhamento==1
        {
        printfln("Cone-de-luz densidades de carga");

        int l_med = 1; //contador p/ inicio da medida abaixo

        //if(L%2==1){l_med = (L-1)/2;} //cadeia par
        //else{l_med = L/2;} //cadeia ímpar

        for(int s_A = l_med; s_A <= L ; s_A++)
        {
                //Estado evoluído

                int j = s_A;

                psi_Evol.position(j);
                auto ket = psi_Evol(j);
                auto bra = dag(prime(ket,"Site"));
                auto N_op = op(sites,"Ntot",j);
                auto Ntot = eltC(bra*N_op*ket).real();
                auto mag_Evol = Ntot;

                psi_GS.position(j);
                auto ket1 = psi_GS(j);
                auto bra1 = dag(prime(ket1,"Site"));
                auto N_op1 = op(sites,"Ntot",j);
                auto Ntot1 = eltC(bra1*N_op1*ket1).real();
                auto mag_GS = Ntot1;

                psi_ini.position(j);
                auto ket2 = psi_ini(j);
                auto bra2 = dag(prime(ket2,"Site"));
                auto N_op2 = op(sites,"Ntot",j);
                auto Ntot2 = eltC(bra2*N_op2*ket2).real();
                auto mag_ini = Ntot2;

                printfln("DensiCarga_Evol");
                cout << "densevol_ev " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << mag_Evol << "\n";
                cout << "densevol_ev_t " << j << " " << lt*del_t << " " << mag_Evol <<"\n";

                printfln("DensiCarga_GS");
                cout << "densevol_GS " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << mag_GS << "\n";
                cout << "densevol_GS_t " << j << " " << lt*del_t << " " << mag_GS <<"\n";

                printfln("DensiCarga_ini");
                cout << "densevol_Ini " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << mag_ini << "\n";
                cout << "densevol_Ini_t " << j << " " << lt*del_t << " " << mag_ini <<"\n";

                printfln("Diferença densiCarga_Evol e mag_GS");
                cout << "densevol_gs " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(mag_Evol-mag_GS,2)) << "\n";
                cout << "densevol_gs_t " << j << " " << lt*del_t << " " << sqrt(pow(mag_Evol - mag_GS,2)) <<"\n";

                printfln("Diferença mag_Evol e mag_ini");
                cout << "densevol_ini " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(mag_Evol - mag_ini,2)) << "\n";
                cout << "densevol_ini_t " << j << " " << lt*del_t << " " << sqrt(pow(mag_Evol - mag_ini,2)) << "\n";


        } //fecha loop da densidade de carga local

        }//close if(a_measure_den_carga==1)


        if(a_measure_den_sz==1) //open if a_measure_den_sz==1
        {

        printfln("Cone-de-luz da magnetização");

        int s_med = 1; //contador p/ inicio da medida abaixo

        //if(L%2==1){s_med = (L-1)/2;} //cadeia par
        //else{s_med = L/2;} //cadeia ímpar

        for(int s_A = s_med; s_A <= L ; s_A++)
        {
            //Estado evoluído

            int j = s_A;

            psi_Evol.position(j);
            auto ket = psi_Evol(j);
            auto bra = dag(prime(ket,"Site"));
            auto N_op = op(sites,"Sz",j);
            auto Ntot = eltC(bra*N_op*ket).real();
            auto mag_Evol = Ntot;

            psi_GS.position(j);
            auto ket1 = psi_GS(j);
            auto bra1 = dag(prime(ket1,"Site"));
            auto N_op1 = op(sites,"Sz",j);
            auto Ntot1 = eltC(bra1*N_op1*ket1).real();
            auto mag_GS = Ntot1;

            psi_ini.position(j);
            auto ket2 = psi_ini(j);
            auto bra2 = dag(prime(ket2,"Site"));
            auto N_op2 = op(sites,"Sz",j);
            auto Ntot2 = eltC(bra2*N_op2*ket2).real();
            auto mag_ini = Ntot2;

            printfln("mag_Evol");
            cout << "magevol_Evolu " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << mag_Evol << "\n";
            cout << "magevol_Evolu_t " << j << " " << lt*del_t << " " << mag_Evol <<"\n";

            printfln("mag_GS");
            cout << "magevol_GSu " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << mag_GS << "\n";
            cout << "magevol_GSu_t " << j << " " << lt*del_t << " " << mag_GS <<"\n";

            printfln("mag_ini");
            cout << "magevol_iniu " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << mag_ini << "\n";
            cout << "magevol_iniu_t " << j << " " << lt*del_t << " " << mag_ini <<"\n";



            printfln("Diferença mag_Evol e mag_GS");
            cout << "magevol_gs " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(mag_Evol-mag_GS,2)) << "\n";
            cout << "magevol_gs_t " << j << " " << lt*del_t << " " << sqrt(pow(mag_Evol - mag_GS,2)) <<"\n";

            printfln("Diferença mag_Evol e mag_ini");
            cout << "magevol_ini " << j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(mag_Evol - mag_ini,2)) << "\n";
            cout << "magevol_ini_t " << j << " " << lt*del_t << " " << sqrt(pow(mag_Evol - mag_ini,2)) << "\n";

        } //fecha loop da densidade de carga local

        }//close if(a_measure_den_sz)


        if(a_measure_corr_carga==1) //if open a_measure_corr_carga==1
        {
        printfln("Cone-de-luz funções de correlação de carga");

        int c_med = 0; //contador p/ inicio da medida abaixo

        if(L%2==1){c_med = (L-1)/2;} //cadeia par
        else{c_med = L/2;} //cadeia ímpar


        for(int i_j = 0; i_j <= (L - c_med); i_j += 1) //abre loop das funções de correlação de carga
        {
            if(i_j == 0)
            {
                cout << "C(" << c_med << "," << c_med << ")" << "\n";

                int i1 = c_med;

                psi_Evol.position(i1);
                auto ket = psi_Evol(i1);
                auto bra = dag(prime(ket,"Site"));
                auto Njop = op(sites,"Ntot*Ntot",i1);
                auto njop = eltC(bra*Njop*ket).real();

                psi_Evol.position(i1);
                auto ket1 = psi_Evol(i1);
                auto bra1 = dag(prime(ket1,"Site"));
                auto Njop1 = op(sites,"Ntot",i1);
                auto njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr = njop - pow(njop1,2);

                printfln("Corre psi_Evol");
                cout << "Corr_cone_psiEvol " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr << "\n";
                cout << "Corr_cone_psiEvol_t " << i_j << " " << lt*del_t << " " << corr << "\n";

                psi_GS.position(i1);
                ket = psi_GS(i1);
                bra = dag(prime(ket,"Site"));
                Njop = op(sites,"Ntot*Ntot",i1);
                njop = eltC(bra*Njop*ket).real();

                psi_GS.position(i1);
                ket1 = psi_GS(i1);
                bra1 = dag(prime(ket1,"Site"));
                Njop1 = op(sites,"Ntot",i1);
                njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr1 = njop - pow(njop1,2);

                printfln("Corre psi_GS");
                cout << "Corr_cone_psiGS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr1 << "\n";
                cout << "Corr_cone_psiGS_t " << i_j << " " << lt*del_t << " " << corr1 << "\n";

                psi_ini.position(i1);
                ket = psi_ini(i1);
                bra = dag(prime(ket,"Site"));
                Njop = op(sites,"Ntot*Ntot",i1);
                njop = eltC(bra*Njop*ket).real();

                psi_ini.position(i1);
                ket1 = psi_ini(i1);
                bra1 = dag(prime(ket1,"Site"));
                Njop1 = op(sites,"Ntot",i1);
                njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr2 = njop - pow(njop1,2);

                printfln("Correlações psi_ini");
                cout << "Corr_cone_psiini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr2 << "\n";
                cout << "Corr_cone_psiini_t " << i_j << " " << lt*del_t << " " << corr2 << "\n";

                printfln("Correlações psi_Evol - psi_GS");
                cout << "cone_diff_Evol_GS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr1,2)) << "\n";
                cout << "cone_diff_Evol_GSt " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr1,2)) << "\n";

                printfln("Correlações psi_Evol - psi_ini");
                cout << "cone_diff_Evol_ini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr2,2)) << "\n";
                cout << "cone_diff_Evol_init " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j == 0)
            if(i_j > 0)
            {

                int i1 = c_med;//posição em que começo a medir a função de correlação
                int j1 = i1 + i_j;

                cout << "C(" << i1 << "," << j1 << ")" << "\n";

                //Para correlações psi_Evol

                auto op_i = op(sites,"Ntot",i1);
                auto op_j = op(sites,"Ntot",j1);

                psi_Evol.position(i1);
                auto psidag = dag(psi_Evol);
                psidag.prime("Link");
                auto li_1 = leftLinkIndex(psi_Evol,i1);
                auto C = prime(psi_Evol(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                    {
                    C *= psi_Evol(k);
                    C *= psidag(k);
                    }
                auto lj = rightLinkIndex(psi_Evol,j1);
                C *= prime(psi_Evol(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                auto result = eltC(C).real();

                psi_Evol.position(i1);
                auto keti = psi_Evol(i1);
                auto brai = dag(prime(keti,"Site"));
                auto Njopi = op(sites,"Ntot",i1);
                auto njopi = eltC(brai*Njopi*keti).real();

                psi_Evol.position(j1);
                auto ketj = psi_Evol(j1);
                auto braj = dag(prime(ketj,"Site"));
                auto Njopj = op(sites,"Ntot",j1);
                auto njopj= eltC(braj*Njopj*ketj).real();

                auto corr = result - njopi*njopj;

                printfln("Correlações psi_Evol");
                cout << "Corr_cone_psiEvol " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr << "\n";
                cout << "Corr_cone_psiEvol_t " << i_j << " " << lt*del_t << " " << corr << "\n";

                psi_GS.position(i1);
                psidag = dag(psi_GS);
                psidag.prime("Link");
                li_1 = leftLinkIndex(psi_GS,i1);
                C = prime(psi_GS(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                    {
                    C *= psi_GS(k);
                    C *= psidag(k);
                    }
                lj = rightLinkIndex(psi_GS,j1);
                C *= prime(psi_GS(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                result = eltC(C).real();

                psi_GS.position(i1);
                keti = psi_GS(i1);
                brai = dag(prime(keti,"Site"));
                Njopi = op(sites,"Ntot",i1);
                njopi = eltC(brai*Njopi*keti).real();

                psi_GS.position(j1);
                ketj = psi_GS(j1);
                braj = dag(prime(ketj,"Site"));
                Njopj = op(sites,"Ntot",j1);
                njopj= eltC(braj*Njopj*ketj).real();

                auto corr1 = result - njopi*njopj;

                printfln("Correlações psi_GS");
                cout << "Corr_cone_psiGS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr1 << "\n";
                cout << "Corr_cone_psiGS_t " << i_j << " " << lt*del_t << " " << corr1 << "\n";

                psi_ini.position(i1);
                psidag = dag(psi_ini);
                psidag.prime("Link");
                li_1 = leftLinkIndex(psi_ini,i1);
                C = prime(psi_ini(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                    {
                    C *= psi_ini(k);
                    C *= psidag(k);
                    }
                lj = rightLinkIndex(psi_ini,j1);
                C *= prime(psi_ini(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                result = eltC(C).real();

                psi_ini.position(i1);
                keti = psi_ini(i1);
                brai = dag(prime(keti,"Site"));
                Njopi = op(sites,"Ntot",i1);
                njopi = eltC(brai*Njopi*keti).real();

                psi_ini.position(j1);
                ketj = psi_ini(j1);
                braj = dag(prime(ketj,"Site"));
                Njopj = op(sites,"Ntot",j1);
                njopj= eltC(braj*Njopj*ketj).real();

                auto corr2 = result - njopi*njopj;

                printfln("Correlações psi_ini");
                cout << "Corr_cone_psiini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr2 << "\n";
                cout << "Corr_cone_psiini_t " << i_j << " " << lt*del_t << " " << corr2 << "\n";

                printfln("Correlações psi_Evol - psi_GS");
                cout << "cone_diff_Evol_GS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr1,2)) << "\n";
                cout << "cone_diff_Evol_GSt " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr1,2)) << "\n";

                printfln("Correlações psi_Evol - psi_ini");
                cout << "cone_diff_Evol_ini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr2,2)) << "\n";
                cout << "cone_diff_Evol_init " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j > 0)

         } //fecha loop das correlações

        } //close if(a_measure_corr_carga)

        if(a_measure_corr_sz==1) //open if a_measure_corr_sz == 1
        {

        printfln("Cone-de-luz funções de correlação de Sz");

        int c_med = 0; //contador p/ inicio da medida abaixo, já defini acima, então é necessário criar um novo

        if(L%2==1){c_med = (L-1)/2;} //cadeia par
        else{c_med = L/2;} //cadeia ímpar


        for(int i_j = 0; i_j <= (L - c_med); i_j += 1) //abre loop das funções de correlação de Sz
        {
            if(i_j == 0) //auto correlações, no mesmo sítio
            {
                cout << "C_sz(" << c_med << "," << c_med << ")" << "\n";

                int i1 = c_med;

                psi_Evol.position(i1);
                auto ket = psi_Evol(i1);
                auto bra = dag(prime(ket,"Site"));
                auto Njop = op(sites,"Sz*Sz",i1);
                auto njop = eltC(bra*Njop*ket).real();

                psi_Evol.position(i1);
                auto ket1 = psi_Evol(i1);
                auto bra1 = dag(prime(ket1,"Site"));
                auto Njop1 = op(sites,"Sz",i1);
                auto njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr = njop - pow(njop1,2);

                printfln("Corre psi_Evol Sz");
                cout << "CorrSz_cone_psiEvol " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr << "\n";
                cout << "CorrSz_cone_psiEvol_t " << i_j << " " << lt*del_t << " " << corr << "\n";

                psi_GS.position(i1);
                ket = psi_GS(i1);
                bra = dag(prime(ket,"Site"));
                Njop = op(sites,"Sz*Sz",i1);
                njop = eltC(bra*Njop*ket).real();

                psi_GS.position(i1);
                ket1 = psi_GS(i1);
                bra1 = dag(prime(ket1,"Site"));
                Njop1 = op(sites,"Sz",i1);
                njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr1 = njop - pow(njop1,2);

                printfln("Corre psi_GS Sz");
                cout << "CorrSz_cone_psiGS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr1 << "\n";
                cout << "CorrSz_cone_psiGS_t " << i_j << " " << lt*del_t << " " << corr1 << "\n";

                psi_ini.position(i1);
                ket = psi_ini(i1);
                bra = dag(prime(ket,"Site"));
                Njop = op(sites,"Sz*Sz",i1);
                njop = eltC(bra*Njop*ket).real();

                psi_ini.position(i1);
                ket1 = psi_ini(i1);
                bra1 = dag(prime(ket1,"Site"));
                Njop1 = op(sites,"Sz",i1);
                njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr2 = njop - pow(njop1,2);

                printfln("Correlações psi_ini Sz");
                cout << "CorrSz_cone_psiini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr2 << "\n";
                cout << "CorrSz_cone_psiini_t " << i_j << " " << lt*del_t << " " << corr2 << "\n";

                printfln("Correlações psi_Evol - psi_GS Sz");
                cout << "coneSz_diff_Evol_GS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr1,2)) << "\n";
                cout << "coneSz_diff_Evol_GSt " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr1,2)) << "\n";

                printfln("Correlações psi_Evol - psi_ini Sz");
                cout << "coneSz_diff_Evol_ini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr2,2)) << "\n";
                cout << "coneSz_diff_Evol_init " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j == 0)
            if(i_j > 0)
            {

                int i1 = c_med;//posição em que começo a medir a função de correlação
                int j1 = i1 + i_j;

                cout << "C_sz(" << i1 << "," << j1 << ")" << "\n";

                //Para correlações psi_Evol

                auto op_i = op(sites,"Sz",i1);
                auto op_j = op(sites,"Sz",j1);

                psi_Evol.position(i1);
                auto psidag = dag(psi_Evol);
                psidag.prime("Link");
                auto li_1 = leftLinkIndex(psi_Evol,i1);
                auto C = prime(psi_Evol(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                    {
                    C *= psi_Evol(k);
                    C *= psidag(k);
                    }
                auto lj = rightLinkIndex(psi_Evol,j1);
                C *= prime(psi_Evol(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                auto result = eltC(C).real();

                psi_Evol.position(i1);
                auto keti = psi_Evol(i1);
                auto brai = dag(prime(keti,"Site"));
                auto Njopi = op(sites,"Sz",i1);
                auto njopi = eltC(brai*Njopi*keti).real();

                psi_Evol.position(j1);
                auto ketj = psi_Evol(j1);
                auto braj = dag(prime(ketj,"Site"));
                auto Njopj = op(sites,"Sz",j1);
                auto njopj= eltC(braj*Njopj*ketj).real();

                auto corr = result - njopi*njopj;

                printfln("Correlações psi_Evol Sz");
                cout << "CorrSz_cone_psiEvol " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr << "\n";
                cout << "CorrSz_cone_psiEvol_t " << i_j << " " << lt*del_t << " " << corr << "\n";

                psi_GS.position(i1);
                psidag = dag(psi_GS);
                psidag.prime("Link");
                li_1 = leftLinkIndex(psi_GS,i1);
                C = prime(psi_GS(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                    {
                    C *= psi_GS(k);
                    C *= psidag(k);
                    }
                lj = rightLinkIndex(psi_GS,j1);
                C *= prime(psi_GS(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                result = eltC(C).real();

                psi_GS.position(i1);
                keti = psi_GS(i1);
                brai = dag(prime(keti,"Site"));
                Njopi = op(sites,"Sz",i1);
                njopi = eltC(brai*Njopi*keti).real();

                psi_GS.position(j1);
                ketj = psi_GS(j1);
                braj = dag(prime(ketj,"Site"));
                Njopj = op(sites,"Sz",j1);
                njopj= eltC(braj*Njopj*ketj).real();

                auto corr1 = result - njopi*njopj;

                printfln("Correlações psi_GS Sz");
                cout << "CorrSz_cone_psiGS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr1 << "\n";
                cout << "CorrSz_cone_psiGS_t " << i_j << " " << lt*del_t << " " << corr1 << "\n";

                psi_ini.position(i1);
                psidag = dag(psi_ini);
                psidag.prime("Link");
                li_1 = leftLinkIndex(psi_ini,i1);
                C = prime(psi_ini(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                    {
                    C *= psi_ini(k);
                    C *= psidag(k);
                    }
                lj = rightLinkIndex(psi_ini,j1);
                C *= prime(psi_ini(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                result = eltC(C).real();

                psi_ini.position(i1);
                keti = psi_ini(i1);
                brai = dag(prime(keti,"Site"));
                Njopi = op(sites,"Sz",i1);
                njopi = eltC(brai*Njopi*keti).real();

                psi_ini.position(j1);
                ketj = psi_ini(j1);
                braj = dag(prime(ketj,"Site"));
                Njopj = op(sites,"Sz",j1);
                njopj= eltC(braj*Njopj*ketj).real();

                auto corr2 = result - njopi*njopj;

                printfln("Correlações psi_ini Sz");
                cout << "CorrSz_cone_psiini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << corr2 << "\n";
                cout << "CorrSz_cone_psiini_t " << i_j << " " << lt*del_t << " " << corr2 << "\n";

                printfln("Correlações psi_Evol - psi_GS Sz");
                cout << "coneSz_diff_Evol_GS " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr1,2)) << "\n";
                cout << "coneSz_diff_Evol_GSt " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr1,2)) << "\n";

                printfln("Correlações psi_Evol - psi_ini Sz");
                cout << "coneSz_diff_Evol_ini " << i_j << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(corr - corr2,2)) << "\n";
                cout << "coneSz_diff_Evol_init " << i_j << " " << lt*del_t << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j > 0)

         } //fecha loop das correlações

        } //close if(a_measure_corr_sz==1)

        if(a_measure_ener_quench==1) //open if a_measure_ener_quench==1
        {
        printfln("Energia que quench colocada no sistema: ");

        if(lt == 0)
        {
            auto overlap0 = innerC(psi_ini,H,psi_ini).real(); //H neste caso é H(t)
            auto overlap1 = innerC(psi_ini,H,psi_ini).real(); //H neste caso é H(t-delta_t)

            cout << "ener_quench " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(overlap0-overlap1,2)) << "\n";
            cout << "ener_quench_t " << lt*del_t << " " << sqrt(pow(overlap0-overlap1,2)) << "\n";
        }
        if(lt > 0)
        {
            //Ground State (t - delta_t) Hamiltonian
            //Usarei o Hamiltoniano H1 para calcular a energia do quench abaixo
            auto ampo1 = AutoMPO(sites);
            for(int i = 1; i <= L; ++i)
                {
                ampo1 += U_t(lt-1),"Nupdn",i;
                }
            for(int b = 1; b < L; ++b)
                {
                ampo1 += -t1,"Cdagup",b,"Cup",b+1;
                ampo1 += -t1,"Cdagup",b+1,"Cup",b;
                ampo1 += -t1,"Cdagdn",b,"Cdn",b+1;
                ampo1 += -t1,"Cdagdn",b+1,"Cdn",b;
                ampo1 += V_t(lt-1),"Ntot",b,"Ntot",b+1;
                }
            auto H1 = toMPO(ampo1);

            auto overlap0 = innerC(psi_GS,H,psi_GS).real(); //H neste caso é H(t)
            auto overlap1 = innerC(psi_GS,H1,psi_GS).real(); //H neste caso é H(t-delta_t)

            cout << "ener_quench " << sqrt(pow(U_t(lt)-U0,2)) << " " << sqrt(pow(overlap0-overlap1,2)) << "\n";
            cout << "ener_quench_t " << lt*del_t << " " << sqrt(pow(overlap0-overlap1,2)) << "\n";

        }//fecha if(lt > 0)

        }//close if(a_measure_ener_quench==1)

        }//fecha if(a_measure_cones_varredura) //relativo a medida de todos os cones de luz



        if(lt == n_loop)
        {
            printfln("\nQuantidades finais\n");

            for(int j = 1; j <= L ; j++)
                {
                    psi_Evol.position(j);
                    auto ket = psi_Evol(j);
                    auto bra = dag(prime(ket,"Site"));
                    auto N_op = op(sites,"Sz",j);
                    auto Ntot = eltC(bra*N_op*ket).real();
                    auto mag_Evol = Ntot;

                    psi_Evol.position(j);
                    ket = psi_Evol(j);
                    bra = dag(prime(ket,"Site"));
                    N_op = op(sites,"Nup",j);
                    Ntot = eltC(bra*N_op*ket).real();
                    auto nup_Evol = Ntot;

                    psi_Evol.position(j);
                    ket = psi_Evol(j);
                    bra = dag(prime(ket,"Site"));
                    N_op = op(sites,"Ndn",j);
                    Ntot = eltC(bra*N_op*ket).real();
                    auto ndn_Evol = Ntot;

                    psi_GS.position(j);
                    ket = psi_GS(j);
                    bra = dag(prime(ket,"Site"));
                    N_op = op(sites,"Sz",j);
                    Ntot = eltC(bra*N_op*ket).real();
                    auto mag_GS = Ntot;

                    psi_GS.position(j);
                    ket = psi_GS(j);
                    bra = dag(prime(ket,"Site"));
                    N_op = op(sites,"Nup",j);
                    Ntot = eltC(bra*N_op*ket).real();
                    auto nup_GS = Ntot;

                    psi_GS.position(j);
                    ket = psi_GS(j);
                    bra = dag(prime(ket,"Site"));
                    N_op = op(sites,"Ndn",j);
                    Ntot = eltC(bra*N_op*ket).real();
                    auto ndn_GS = Ntot;

                    cout << "#Nevol_Ntot " << j << " " << nup_Evol+ndn_Evol  << "\n";
                    cout << "#Nevol_Nup " << j << " " << nup_Evol << "\n";
                    cout << "#Nevol_Ndn " << j << " " << ndn_Evol << "\n";
                    cout << "#Nevol_magz " << j << " " << mag_Evol << "\n";

                    cout << "#Ngs_Ntot " << j << " " << nup_GS+ndn_GS  << "\n";
                    cout << "#Ngs_Nup " << j << " " << nup_GS << "\n";
                    cout << "#Ngs_Ndn " << j << " " << ndn_GS << "\n";
                    cout << "#Ngs_magz " << j << " " << mag_GS << "\n";

                    auto d_dens = sqrt(pow(nup_Evol + ndn_Evol - nup_GS - ndn_GS,2));
                    auto d_mag = sqrt(pow(mag_GS-mag_Evol,2));

                    cout << "dif_evol_GS_n " << j << " " << d_dens << "\n";
                    cout << "dif_evol_GS_m " << j << " " << d_mag << "\n";

                }//fecha for(int j = 1; j <= L ; j++)

        }//fecha if(lt == n_loop)



    }//fecha loop total


    //Falta implementar a rotina que salva psi



    writeToFile("sites.h5",sites);
    writeToFile<MPS>("psi_Evol.h5",psi_Evol);

    //readFromFile("sites",sites);
    //auto psi_Evol = readFromFile<MPS>("psi");





return 0;
}

