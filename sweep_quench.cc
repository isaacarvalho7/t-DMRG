/*!
    Código para o modelo de Hubbard estendido 1D, construído utilizando a biblioteca ITensor.
    Configure o input_file para definir os parâmetros do cálculo e ativar as medidas desejadas.
    O código avalia as correlações de carga e spin, os parâmetros de ordem, CDW e SDW,
    e o emaranhamento durante e após o quench de varredura ou sweep quench.

    \author Isaac Martins Carvalho.
    \since 2/07/21
    \version 1.0
 */
#include "itensor/all.h"
#include "stdlib.h"
#include "omp.h"
#include "thread"
#include "cmath"
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

using namespace itensor;
using std::vector;
using namespace std::chrono;
using namespace std;

int main(int argc, char* argv[])
{
   //!Parâmetros para importar os dados do input_file

   /*!
      argc será definido pelo compilador C++
      igual ao número de entradas para o programa e argv é uma matriz dessas entradas.
      argv[0] é sempre o nome do programa.
      argv[1] é a primeira entrada não trivial na linha de comando.

      Então ao rodar o código usando "./program filename", teremos

      argv[0] == "program"

      argv[1] == "filename"

      A declaração if(argc != 2) garante que fornecemos o que o nome do arquivo.
      O InputGroup (argv [1], "input") constrói um InputGroup a partir do arquivo de entrada
      (cujo nome é armazenado em argv [1]). A string "input" nomeia a coleção de entradas
      que queremos ler.
    */

    //! if(argc != 2)
    /*!
       Garante que fornecemos o que o nome do arquivo.
     */
    if(argc != 2)
    {
        printfln("Usage: %s inputfile",argv[0]);
        return 0;
    }

    //! InputGroup(argv[1],"input)
    /*!
       constrói um InputGroup a partir do arquivo de entrada
       (cujo nome é armazenado em argv [1]). A string "input" nomeia a coleção de entradas
       que queremos ler
     */
    auto input = InputGroup(argv[1],"input");

    /*!
       A funções getInt ou getReal extrai o parâmetro "a" int ou Real do input_file: a configurado no input_file.
       b atribuído a variável caso "a" não seja definido.
     */


    auto L = input.getInt("L");
    auto Npart = input.getInt("Npart",L);
    auto sweeps = input.getInt("nsweeps");
    auto t1 = input.getReal("t1",1.);
    auto U0 = input.getReal("U0",0.);
    auto V0 = input.getReal("V0",0.);
    auto Uf = input.getReal("Uf",1.);
    auto Vf = input.getReal("Vf",2.0);

    auto h_z = input.getReal("h_z",0.0);
    auto tau_U = input.getReal("tau_U",0.01);

    auto tau_V = (tau_U*(sqrt(pow(Uf-U0,2))))/(sqrt(pow(Vf-V0,2))); ///< Definição de tau_U

    printfln("tau_V = ",tau_V);

    //Parâmetro para evolução após o quench de varredura.
    auto tstep = input.getReal("tstep",0.05);
    auto cutoff = input.getReal("cutoff",1E-12);
    auto ttstep = input.getReal("ttstep",tstep);

    auto time_evolution = input.getReal("time_evolution",0.0);

    auto quiet = input.getYesNo("quiet",false);
    println(sweeps);

    //////////////////////////////////////
    ///          ATENÇÃO
    //////////////////////////////////////

    /*!
       Configure no input igual 1 para ligar e
       igual 0 para desligar.
     */


    int performs_sweep_quench = input.getInt("performs_sweep_quench",0);

    int DRMG_varredura_measure = input.getInt("DRMG_varredura_measure",0);
    int m_cdw_measure = input.getInt("m_cdw_measure",1);
    int m_sdw_measure = input.getInt("m_sdw_measure",1);
    int emaranhamento_parcial_measure = input.getInt("emaranhamento_parcial_measure",1);

    int a_measure_emaranhamento = input.getInt("a_measure_emaranhamento",0);
    int a_measure_den_carga = input.getInt("a_measure_den_carga",0);
    int a_measure_den_sz = input.getInt("a_measure_den_sz",0);
    int a_measure_corr_carga = input.getInt("a_measure_corr_carga",0);
    int a_measure_corr_sz = input.getInt("a_measure_corr_sz",0);
    int a_measure_ener_quench = input.getInt("a_measure_ener_quench",0);

    int turn_on_totalEvol_after = input.getInt("turn_on_totalEvol_after",0);

    int turn_emaranhamento = input.getInt("turn_emaranhamento",0);
    int turn_dens_carga = input.getInt("turn_dens_carga",0);
    int turn_corr_carga = input.getInt("turn_corr_carga",0);
    int turn_dens_sz = input.getInt("turn_dens_sz",0);
    int turn_corr_sz = input.getInt("turn_corr_sz",0);

    int salve_psi_final_sweep_quench = input.getInt("salve_psi_final_sweep_quench",0);
    int salve_psi_evolution = input.getInt("salve_psi_evolution",0);

    ////////////////////////////
    ///          CODE
    ////////////////////////////

    auto sites = Electron(L); ///<Inicializa o sites para ser do tipo EletronSite
    auto state = InitState(sites); ///< Inicializa MPS fornecendo um tipo de SiteSet

    int p = Npart; //Filling

    //! Função para configurar o preenchimento da cadeia
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


    auto psi0 = MPS(state); ///< Salva o estado inicial em psi0, utilizado para o DMRG a seguir.

    //! Definição da função sgn
    auto sgn_U = Uf-U0;
    auto sgn_V = Vf-V0;

    if(sgn_U < 0.){sgn_U = -1.;}
    if(sgn_U == 0.){sgn_U = 0.;}
    if(sgn_U > 0){sgn_U = 1.;}

    if(sgn_V < 0.){sgn_V = -1.;}
    if(sgn_V == 0.){sgn_V = 0.;}
    if(sgn_V > 0){sgn_V = 1.;}

    double time_f = (sqrt(pow(Uf - U0,2))*tau_U); ///< time_f é o tempo total do sweep quench

    printfln("\ntime_f = ", time_f);

    //! Número de pontos mínimos para o sweep quench.
    /*!
       Condition so that we have at least 20 (p_t = 20 below) points on the graph
       Increases the tstep precision so that this occurs when tau_U is too small.
     */
    int p_t = 20; ///<Número mínimo de ponto no gráfico.
    if(time_f <= 0.09){tstep = time_f/p_t;}
    if(time_f > 0.09 && time_f < 0.2){tstep = time_f/p_t;}

    int n_loop = int(time_f/tstep); ///< Valor de n_loop, utilizado para configurar o loop da varredura.


    //! Definição de U(t) e V(t) utilizados no sweep quench.
    /*!
       O valor de U(t) e V(t) é atualizado em cada step da varredura.
       A função abaixo configura um valor para U(t) e V(t) em cada step de
       tempo.
     */

    Real del_t = tstep; ///< Valor de del_t utilizado para construir os vetores U(t) e V(t).

    Vector U_t(n_loop + 1) , V_t(n_loop + 1);

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


    //!Ground State (GS) Hamiltonian
    /*!
        Definição do Hamiltoniano do estado inicial.
        Esse Hamiltoniano será utilizado para realizar DRMG a fim
        de rastrear o estado fundamental com a maior precisão possível.
     */
    auto ampo = AutoMPO(sites); ///<Inicializa o AutoMPO dos sites
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

    //! Realiza DMRG
    /*!
        https://itensor.org/docs.cgi?vers=cppv3&page=classes/dmrg

        Recebe as entradas H0, psi0 e sweeps.
        Retorna a energia energyGS do ground state psi.
     */
    auto [energyGS,psi] = dmrg(H0,psi0,sweeps,{"Quiet",true}); ///< DMRG, retorna energy GS e psi

    auto psi_ini = psi; ///< Salva o ground state psi em psi_ini
    auto psi_Evol = psi; ///< Salva o ground state psi em psi_Evol que será usado no sweep quench.



	//////////////////////////////////
	///	Finite time quantum quench
	//////////////////////////////////


	/*!
	    Rotina para os sweep quench. Nessa rotina está incluso
	    uma série de condições, localizada no escopo de cada medida,
	    que possuem a finalidade de controlar o que será medido durante a varredura.
	    A configuração do que será medido se encontra no arquivo input_file.
	    Além disso, é avaliado o tempo das rotinas DRMG e t-DRMG em cada instante da varredura.
	 */

    if(performs_sweep_quench == 1)
    {
    for(int lt = 0; lt <= n_loop; lt ++)
    {
        //!Ground State (GS) Hamiltonian
        /*!
           O Hamiltoniano abaixo é construído mediante os valores de U(t) e V(t).
           Ele será utilizado no DRMG para avaliar os estados de equilíbrio em cada instante de tempo,
           ou seja, os grounds states relacionados aos valores de U(t) e V(t) em cada instante de tempo.
         */

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

        steady_clock::time_point t0 = steady_clock::now(); //time from here

        //!DMRG da varredura desligado?
        /*!
           No caso em que não se realiza DMRG da varredura, a variavel psi_GS,
           relacionada aos estados de equilíbrio, recebe psi0.
         */
        auto psi_GS=psi0;
        if(DRMG_varredura_measure==1)
        {
        //!Realiza DMRG, retorna energy_GS e psi_GS
        auto [energy_GS,psi_GS] = dmrg(H,psi0,sweeps,{"Quiet",true}); ///<DMRG da varredura


        steady_clock::time_point t2 = steady_clock::now(); //final calculation time
        duration<double> time_span = duration_cast<duration<double>>(t2 - t0);
        cout << "#time_drmg  " << sqrt(pow(U_t(lt)-U0,2)) << " " << time_span.count() <<"\n";
        cout << "#time_drmg_t  " << lt*del_t << " " << time_span.count() <<"\n";


        //! Energia do estado de equilíbrio ao final do DMRG
        auto h_0 = inner(psi_GS,H,psi_GS); ///< Energia de psi_GS é igual h_0

        printfln("<psi_GS,H,psi_GS>");
        cout << "#h_0 " << sqrt(pow(U_t(lt)-U0,2)) << " " <<  h_0 <<"\n";
        cout << "#h_0_t " << lt*del_t << " " << h_0 <<"\n";

        printfln("Ground state energy");
        cout << "%energy_GS " << sqrt(pow(U_t(lt)-U0,2)) << " " <<  energy_GS <<"\n";
        cout << "%energy_GS_t " << lt*del_t << " " << energy_GS <<"\n";

        }

        //!Contrução dos gates da evolução
        /*!
            Abaixo é construído os gates utilizados na decomposição Trotter.
            Esses gates irão compor a função que realiza as evoluções do sweep quench.
         */
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
        //Finish the bulding of gates

        //! Função gateTEvol que realiza a evolução de psi_Evol.
        /*!
            A função gateTEvol é definida gateTEvol(gates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true}).

            https://itensor.org/docs.cgi?vers=cppv3&page=formulas/tevol_trotter

            Em nosso caso, configuramos o tempo de evolução ttotal para o mesmo valor de tstep.
            O sweep quench corresponde a uma sucessão de evoluções de tempo tstep, onde em cada evolução
            os valores de U(t) e V(t) são atualizados até completar o tempo total da varredura time_f.
         */
        gateTEvol(gates,tstep,tstep,psi_Evol,{"Quiet",true,"Cutoff=",cutoff,"Verbose=",true,"UseSVD=",true,"SVDMethod=","gesdd"}); ///< Realiza t-DMRG sobre os estados evoluídos ou de não equilíbrio.

		//! Saídas do cálculo:
		//Print o instante de tempo
	    printfln("\n\ntime do loop = ", lt*del_t); ///< Instante de tempo da varredura.

        //Para entender como U e V estão variando no tempo
        cout << "#!U " << lt*del_t << " " << U_t(lt) <<"\n"; ///< Valor de U(t) no tempo
        cout << "#!V " << lt*del_t << " " << V_t(lt) <<"\n"; ///< Valor de V(t) no tempo

        //Compute overlap
        printfln("Overlap <psi_gs|H|psi_gs> - <psi_evol|H|psi_evol> ");
        auto overlap_evol = innerC(psi_Evol,H,psi_Evol).real(); ///< Overlap de psi_Evol com H(t)
        auto overlap_gs = innerC(psi_GS,H,psi_GS).real(); ///< Overlpa de psi_GS com H(t)

        cout << "#over_ut " << sqrt(pow(U_t(lt)-U0,2)) << " " << overlap_gs-overlap_evol <<"\n";
        cout << "#over_vt " << sqrt(pow(V_t(lt)-V0,2)) << " " << overlap_gs-overlap_evol <<"\n";
        cout << "#over_time " << lt*del_t << " " << overlap_gs-overlap_evol <<"\n";


        //! Medição do parâmetro m_cdw para psi_Evol e psi_GS
        if(m_cdw_measure == 1)
        {
        //Parâmetro m_cdw
        printfln("Parâmetro m_cdw psi_Evol");
        auto mc = 0.;
        for(int j = 1; j <= L ; j++)
        {
            psi_Evol.position(j);
            auto ket = psi_Evol(j);
            auto bra = dag(prime(ket,"Site"));
            auto N_op = op(sites,"Ntot",j);
            auto Ntott = eltC(bra*N_op*ket).real();
            mc += pow(-1,j)*(Ntott - 1.0);
        }
        cout << "m_cdw_evol " << sqrt(pow(U_t(lt)-U0,2)) << " " << mc/L <<"\n";
        cout << "m_cdw_evolt " << lt*del_t << " " << mc/L <<"\n";

        printfln("Parâmetro m_cdw psi_GS");
        mc = 0.;
        for(int j = 1; j <= L ; j++)
        {
            psi_GS.position(j);
            auto ket = psi_GS(j);
            auto bra = dag(prime(ket,"Site"));
            auto N_op = op(sites,"Ntot",j);
            auto Ntott = eltC(bra*N_op*ket).real();
            mc += pow(-1,j)*(Ntott - 1.0);
        }
        cout << "m_cdwGS " << sqrt(pow(U_t(lt)-U0,2)) << " " << mc/L <<"\n";
        cout << "m_cdwGSt " << lt*del_t << " " << mc/L <<"\n";
        }

        //! Medição do parâmetro m_sdw para psi_Evol e psi_GS
        if(m_sdw_measure == 1)
        {
        //Parâmetro m_sdw
        printfln("Parâmetro m_sdw psi_Evol");
        auto msc = 0.;
        for(int j = 1; j <= L ; j++)
        {
            psi_Evol.position(j);
            auto ket = psi_Evol(j);
            auto bra = dag(prime(ket,"Site"));
            auto N_op = op(sites,"Sz",j);
            auto Ntott = eltC(bra*N_op*ket).real();
            msc += pow(-1,j)*(Ntott);
        }
        cout << "m_sdw_evol " << sqrt(pow(U_t(lt)-U0,2)) << " " << msc/L <<"\n";
        cout << "m_sdw_evolt " << lt*del_t << " " << msc/L <<"\n";

        printfln("Parâmetro m_sdw psi_GS");
        msc = 0.;
        for(int j = 1; j <= L ; j++)
        {
            psi_GS.position(j);
            auto ket = psi_GS(j);
            auto bra = dag(prime(ket,"Site"));
            auto N_op = op(sites,"Sz",j);
            auto Ntott = eltC(bra*N_op*ket).real();
            msc += pow(-1,j)*(Ntott);
        }
        cout << "m_sdwGS " << sqrt(pow(U_t(lt)-U0,2)) << " " << msc/L <<"\n";
        cout << "m_sdwGSt " << lt*del_t << " " << msc/L <<"\n";
        }


        //! Medição entropia de emaranhamento de von Neumann para psi_Evol e psi_GS
        if(emaranhamento_parcial_measure == 1)
        {
        printfln("Emaranhamento GS L/2");
        int s_A = 1;
        if(L%2 == 0){s_A = L/2;}
        else{s_A = (L-1)/2;}

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

        cout << "#svngs_l/2 " << b << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << SvN <<"\n";
        cout << "#svngs_l/2t " << b << " " << lt*del_t << " " << SvN <<"\n";

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
        cout << "#svnevol/2 " << b1 << " " << sqrt(pow(U_t(lt)-U0,2)) << " " << SvN1 <<"\n";
        cout << "#svnevol/2t "<< b1 << " " << lt*del_t << " " << SvN1 <<"\n";
        }

        //! Medição density plot do emaranhamento para psi_Evol, psi_GS e psi_inicial.
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

        //! Medição density plot para densidade local de carga de psi_Evol, psi_GS e psi_inicial.
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

        //! Medição density plot para a magnetização Sz local de psi_Evol, psi_GS e psi_inicial.
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

        //! Medição density plot para as correlações de carga de psi_Evol, psi_GS e psi_inicial.
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

                cout << "c_carga_all_t " << i_j << " " << lt*del_t << " " << corr << "\n";

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
                cout << "c_carga_all_t " << i_j << " " << lt*del_t << " " << corr << "\n";


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

        //! Medição density plot para as correlações da magnetização Sz de psi_Evol, psi_GS e psi_inicial.
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

        //! Mediação da energia colocado pelo quench no equilíbrio.
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

        //!Quantidades medidas ao final da varredura.
        /*!
           Ao final do quench de varredura são feitas as seguintes as medições das densidades locais de
           spin up e down ao longo da cadeia.
         */
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

    }//fecha if(performs_sweep_quench == 1)

    printfln("\n\nFINAL DO QUENCH DE VARREDURA\n\n");


    //! A função abaixo salva e lê a função de onda ao final da varredura
    /*!
        No input_file o usuário deve escolher se deseja salvar ou não a função de onda ao final
        da varredura. É necessário salvar a função de onda e o sites para que seja possível reutilizá-los depois.
     */
    if(salve_psi_final_sweep_quench == 1) //Salvar ou não salvar função de onda
    {
    string site_fl = format("sitesFinalSweepQuench_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U); ///< string para salvar sites
    string psi_fl = format("psiFinalSweepQuench_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U); ///< string para salvar função de onda psi_Evol

    writeToFile(site_fl,sites); ///< Salva na pasta local os sites
    writeToFile<MPS>(psi_fl,psi_Evol); ///< Salva na pasta local a função de onda psi_Evol

    readFromFile(site_fl,sites); ///< Lê na pasta local os sites
    psi_Evol = readFromFile<MPS>(psi_fl); ///< Lê na pasta local a função de onda psi_Evol.
    }

    /*!
       Salva psi_Evol em psi_Evol_ini, pois vou utilizar
       psi_Evol para continuar a evolução após o sweep quench
     */
    auto psi_Evol_ini = psi_Evol; ///< Salva psi_Evol em psi_Evol_ini

    //!Configure o input_file para ativar a opção de evolução pós quench de varredura.
    if(turn_on_totalEvol_after ==1)
    {
    //////////////////////////
    ///	AFTER SWEEP QUENCH  //
    //////////////////////////

    /*!
        Abaixo definimos os gates da evolução fora do loop em que as medidas serão realizadas, pois
        não é mais necessário atualizar os valores dos parâmetros U(t) e V(t). Neste caso, os valores das interações
        serão mantidos fixos, sendo Uf e Vf os valores atribuídos ao Hamiltoniano da evolução.
     */
    int lt = n_loop;

    //!Printa U(t) = Uf e V(t)=Vf ao final da varredura.
    cout << "\nU_t(" << lt << ") = " << U_t(lt) << "\n"; ///< U(n_loop) = Uf
    cout << "V_t(" << lt << ") = " << V_t(lt) << "\n"; ///< V(n_loop) = Vf


    //!Contrução dos gates da evolução
    /*!
        Abaixo é construído os gates utilizados na decomposição Trotter.
        Esses gates irão compor a função que realiza a evolução do pós sweep quench.
     */
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


    double ttime = time_evolution; ///< ttime recebe o tempo total da evolução pós sweep quench.

    //!Variável do loop da evolução
    /*!
       Abaixo o tempo total da evolução é fracionado ttstep, que é o intervalo de tempo
       no qual as quantidade serão medidas ao longo da evolução.
     */
    //Variável do for, para determinar quantos loops serão executados
    auto tloop = int((ttime/ttstep)+0.5); ///< tloop definida para ajustar possíveis arredondamento.

    printfln("\ntloop = ", tloop);

    //!Rotina para evolução pós sweep quench.
    /*!
        A rotina de evolução é fracionada em steps para que possamos medir a quantidades ao longo da evolução.
        Em cada loop a função de onda é evoluída sob ação de Uf e Vf por um tempo total ttstep.
        A soma do tempo ttstep de todos loops é igual ao tempo total. A rotina está configurada para
        quantificar o tempo gasto para evoluir a função de onda e para concluir
        todas medidas.
     */
    for(int time = 0; time <= tloop; time += 1)
    {
        printfln("\ntime = ", time*ttstep); ///< printa o instante de tempo da evolução.


        /*!
            O gateTEvol recebe psi_Evol e o evolui por um tempo total igual a ttstep, a uma precisão ditada por tstep.
            A condição if(tstep < ttstep){tstep = ttstep;} abaixo serve para resolver o problema em que o tempo total da evolução é superior a precisão,
            que resulta em travamento do cálculo.
         */
        if(tstep < ttstep){tstep = ttstep;}
        steady_clock::time_point t0 = steady_clock::now(); //time from here
        if(time > 0)
        {
            gateTEvol(gates,ttstep,tstep,psi_Evol, {"Quiet",true,"Cutoff=",cutoff,"Verbose=",true,"UseSVD=",true,"SVDMethod=","gesdd"}); ///< Realiza t-DRMG sobre a função de onda.

        }
        steady_clock::time_point t2 = steady_clock::now(); //final calculation time
        duration<double> time_span = duration_cast<duration<double>>(t2 - t0);
        cout << "#gatevol_after_t  " << time*ttstep << " " << time_span.count() <<"\n";


        steady_clock::time_point t3 = steady_clock::now(); //time from here

        //!Rotina para salvar a função de onda da evolução
        /*!
            As rotinas abaixo são definidas para salvar a função de onda na metade do tempo total da evolução e ao final da evolução.
            Caso deseje utilizá-las, configure o input_file para isso.
         */
         //!Salva a função na metade do tempo total da evolução.
        if(salve_psi_evolution == 1)
        {
        if(tloop%2 == 0) //caso tloop par
        {
            if(time == tloop/2)
            {
                string site_ev = format("sites_Evol_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U,"_time_",time*ttstep);
                string psi_ev = format("psi_Evol_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U,"_time_",time*ttstep);
                writeToFile(site_ev,sites);
                writeToFile<MPS>(psi_ev,psi_Evol);
            }
        }
        else
        {
            if(time == (tloop-1)/2) //caso tloop ímpar
            {
                string site_ev = format("sites_Evol_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U,"_time_",time*ttstep);
                string psi_ev = format("psi_Evol_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U,"_time_",time*ttstep);
                writeToFile(site_ev,sites);
                writeToFile<MPS>(psi_ev,psi_Evol);
            }
        } //fecha else

        //! Salva a função de onda ao final da evolução.
        if(time==tloop) //Salvar ao final da evolução a função de onda e o sites
        {
            string site_ev = format("sites_Evol_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U,"_time_",time*ttstep);
            string psi_ev = format("psi_Evol_L_",L,"_Npart_",Npart,"_V0_",V0,"_U0_",U0,"_Uf_",Uf,"_Vf_",Vf,"_tau_U_",tau_U,"_time_",time*ttstep);
            writeToFile(site_ev,sites);
            writeToFile<MPS>(psi_ev,psi_Evol);
        }
        }//fecha if(salve_psi_evolution == 1)

        //! Avalia pontos para o density plot para a entropia de emaranhamento de psi_Evol
        if(turn_emaranhamento == 1)
        {
        //Para começar a medir no centro da cadeia
        printfln("Cone-de-luz entropia");

        int i_med = 1; //contador p/ inicio da medida abaixo

        //if(L%2==1){i_med = (L-1)/2;} //cadeia par
        //else{i_med = L/2;} //cadeia ímpar

        for(int s_A = i_med; s_A < L ; s_A++) // < L, pois precisamos subsistema A e B
        {
            //Entropia psi_Evol
            auto b1 = s_A;

            psi_Evol.position(b1);
            auto l1 = leftLinkIndex(psi_Evol,b1);
            auto s1 = siteIndex(psi_Evol,b1);
            auto [U1,S1,V1] = svd(psi_Evol(b1), {l1,s1});
            auto u1 = commonIndex(U1,S1);
            Real SvN1 = 0.;
            for(auto n : range1(dim(u1)))
            {
                auto Sn = elt(S1,n,n);
                auto p = sqr(Sn);
                if(p > 1E-12) SvN1 += -p*log(p);
            }
            printfln("Entropia psi_Evol");
            cout << "#af_svntime_1 "<< b1 << " " << time*ttstep << " " << SvN1 <<"\n";

            //Emaranhamento estado inicial psi_ini
            auto b2 = s_A;

            psi_Evol_ini.position(b2);
            auto l2 = leftLinkIndex(psi_Evol_ini,b2);
            auto s2 = siteIndex(psi_Evol_ini,b2);
            auto [U2,S2,V2] = svd(psi_Evol_ini(b2), {l2,s2});
            auto u2 = commonIndex(U2,S2);
            Real SvN2 = 0.;
            for(auto n : range1(dim(u2)))
            {
                auto Sn = elt(S2,n,n);
                auto p = sqr(Sn);
                if(p > 1E-12) SvN2 += -p*log(p);
            }
            printfln("Diference Entropia psi_Evol and psi_Evol_ini");
            cout << "#af_dife_svn2_conet " << b2 << " " << time*ttstep << " " << sqrt(pow(SvN1-SvN2,2)) <<"\n";

        } //fim loop da entropia
        }

        //!Avalia pontos para o density plot da densidade de carga de psi_Evol
        if(turn_dens_carga==1)
        {
        printfln("Cone-de-luz da densidade de carga local");

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
            auto N_op = op(sites,"Ntot",j);
            auto Ntot = eltC(bra*N_op*ket).real();
            auto mag_Evol = Ntot;

            printfln("Densidade local de carga");
            cout << "af_local_t " << j << " " << time*ttstep << " " << mag_Evol <<"\n";

            psi_Evol_ini.position(j);
            ket = psi_Evol_ini(j);
            bra = dag(prime(ket,"Site"));
            N_op = op(sites,"Ntot",j);
            Ntot = eltC(bra*N_op*ket).real();
            auto mag_Evol_ini = Ntot;


            printfln("Diferença psi_Evol e psi_Evol_ini local de carga");
            cout << "diff_af_local_t " << j << " " << time*ttstep << " " << sqrt(pow(mag_Evol - mag_Evol_ini,2)) <<"\n";


        } //fecha loop da densidade de carga local
        }

        //!Avalia pontos para o density plot das correlações de carga de psi_Evol
        if(turn_corr_carga == 1)
        {
        printfln("Cone-de-luz funções de correlação de carga");

        int c_med = 0; //contador p/ inicio da medida abaixo

        if(L%2==1)
        {
            c_med = (L-1)/2;   //cadeia par
        }
        else
        {
            c_med = L/2;   //cadeia ímpar
        }


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
                cout << "af_Corr_cone_psiEvol_t " << i_j << " " << time*ttstep << " " << corr << "\n";
                cout << "c_carga_all_t " << i_j << " " << time*ttstep+n_loop*tstep << " " << corr << "\n";

                psi_Evol_ini.position(i1);
                ket = psi_Evol_ini(i1);
                bra = dag(prime(ket,"Site"));
                Njop = op(sites,"Ntot*Ntot",i1);
                njop = eltC(bra*Njop*ket).real();

                psi_Evol_ini.position(i1);
                ket1 = psi_Evol_ini(i1);
                bra1 = dag(prime(ket1,"Site"));
                Njop1 = op(sites,"Ntot",i1);
                njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr1 = njop - pow(njop1,2);

                printfln("Correlações psi_Evol - psi_Evol_ini");
                cout << "af_cone_diff_Evol_init " << i_j << " " << time*ttstep << " " << sqrt(pow(corr - corr1,2)) << "\n";

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
                cout << "af_Corr_cone_psiEvol_t " << i_j << " " << time*ttstep << " " << corr << "\n";
                cout << "c_carga_all_t " << i_j << " " << time*ttstep+n_loop*tstep << " " << corr << "\n";

                psi_Evol_ini.position(i1);
                psidag = dag(psi_Evol_ini);
                psidag.prime("Link");
                li_1 = leftLinkIndex(psi_Evol_ini,i1);
                C = prime(psi_Evol_ini(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                {
                    C *= psi_Evol_ini(k);
                    C *= psidag(k);
                }
                lj = rightLinkIndex(psi_Evol_ini,j1);
                C *= prime(psi_Evol_ini(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                result = eltC(C).real();

                psi_Evol_ini.position(i1);
                keti = psi_Evol_ini(i1);
                brai = dag(prime(keti,"Site"));
                Njopi = op(sites,"Ntot",i1);
                njopi = eltC(brai*Njopi*keti).real();

                psi_Evol_ini.position(j1);
                ketj = psi_Evol_ini(j1);
                braj = dag(prime(ketj,"Site"));
                Njopj = op(sites,"Ntot",j1);
                njopj= eltC(braj*Njopj*ketj).real();

                auto corr2 = result - njopi*njopj;

                printfln("Correlações psi_Evol - psi_Evol_ini");
                cout << "af_cone_diff_Evol_init " << i_j << " " << time*ttstep << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j > 0)

        } //fecha loop das correlações
        }

        //!Avalia pontos para o density plot da magnetização Sz local de psi_Evol.
        if(turn_dens_sz==1)
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


            psi_Evol_ini.position(j);
            auto ket2 = psi_Evol_ini(j);
            auto bra2 = dag(prime(ket2,"Site"));
            auto N_op2 = op(sites,"Sz",j);
            auto Ntot2 = eltC(bra2*N_op2*ket2).real();
            auto mag_ini = Ntot2;

            printfln("mag_Evol");
            cout << "af_magevol_Evolu_t " << j << " " << time*ttstep << " " << mag_Evol <<"\n";

            printfln("Diferença mag_Evol e mag_ini");
            cout << "af_magevol_ini_t " << j << " " << time*ttstep << " " << sqrt(pow(mag_Evol - mag_ini,2)) << "\n";

        } //fecha loop da densidade de carga local
        }

        //!Avalia pontos para o density plot das correlações Sz ao longo da cadeia.
        if(turn_corr_sz ==1)
        {
        printfln("Cone-de-luz funções de correlação de Sz");

        int cs_med = 0; //contador p/ inicio da medida abaixo, já defini acima, então é necessário criar um novo

        if(L%2==1)
        {
            cs_med = (L-1)/2;   //cadeia par
        }
        else
        {
            cs_med = L/2;   //cadeia ímpar
        }


        for(int i_j = 0; i_j <= (L - cs_med); i_j += 1) //abre loop das funções de correlação de Sz
        {
            if(i_j == 0) //auto correlações, no mesmo sítio
            {
                cout << "C_sz(" << cs_med << "," << cs_med << ")" << "\n";

                int i1 = cs_med;

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
                cout << "af_CorrSz_cone_psiEvol_t " << i_j << " " << time*ttstep << " " << corr << "\n";

                psi_Evol_ini.position(i1);
                ket = psi_Evol_ini(i1);
                bra = dag(prime(ket,"Site"));
                Njop = op(sites,"Sz*Sz",i1);
                njop = eltC(bra*Njop*ket).real();

                psi_Evol_ini.position(i1);
                ket1 = psi_Evol_ini(i1);
                bra1 = dag(prime(ket1,"Site"));
                Njop1 = op(sites,"Sz",i1);
                njop1 = eltC(bra1*Njop1*ket1).real();

                auto corr2 = njop - pow(njop1,2);

                printfln("Correlações psi_Evol - psi_Evol_ini Sz");
                cout << "af_coneSz_diff_Evol_init " << i_j << " " << time*ttstep  << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j == 0)
            if(i_j > 0)
            {

                int i1 = cs_med;//posição em que começo a medir a função de correlação
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
                cout << "af_CorrSz_cone_psiEvol_t " << i_j << " " << time*ttstep << " " << corr << "\n";

                psi_Evol_ini.position(i1);
                psidag = dag(psi_Evol_ini);
                psidag.prime("Link");
                li_1 = leftLinkIndex(psi_Evol_ini,i1);
                C = prime(psi_Evol_ini(i1),li_1)*op_i;
                C *= prime(psidag(i1),"Site");
                for(int k = i1+1; k < j1; ++k)
                {
                    C *= psi_Evol_ini(k);
                    C *= psidag(k);
                }
                lj = rightLinkIndex(psi_Evol_ini,j1);
                C *= prime(psi_Evol_ini(j1),lj)*op_j;
                C *= prime(psidag(j1),"Site");
                result = eltC(C).real();

                psi_Evol_ini.position(i1);
                keti = psi_Evol_ini(i1);
                brai = dag(prime(keti,"Site"));
                Njopi = op(sites,"Sz",i1);
                njopi = eltC(brai*Njopi*keti).real();

                psi_Evol_ini.position(j1);
                ketj = psi_Evol_ini(j1);
                braj = dag(prime(ketj,"Site"));
                Njopj = op(sites,"Sz",j1);
                njopj= eltC(braj*Njopj*ketj).real();

                auto corr2 = result - njopi*njopj;

                printfln("Correlações psi_Evol - psi_Evol_ini Sz");
                cout << "af_coneSz_diff_Evol_init " << i_j << " " << time*ttstep << " " << sqrt(pow(corr - corr2,2)) << "\n";

            }//fecha if(i_j > 0)

        } //fecha loop das correlações de magnetização
        }

        steady_clock::time_point t4 = steady_clock::now(); //final calculation time of cone measure
        duration<double> time_span_0 = duration_cast<duration<double>>(t4 - t3);
        cout << "#cone_after_t  " << time*ttstep << " " << time_span_0.count() <<"\n";

    } //fecha loop de evolução pós quench de varredura no tempo
    }//fecha if(turn_on_totalEvol_after == 1) que realizada toda medida pós varredura



return 0;
}

