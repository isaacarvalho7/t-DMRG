#include "itensor/all.h"
#include "stdlib.h"
#include "omp.h"
#include "thread"
#include "cmath"
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <list>
#include <bitset>
#include <sstream>
#include <string>
#include <armadillo>



#define MAX_DIGITS 10

using namespace itensor;
using std::vector;
using namespace std::chrono;
using namespace std;


string IntToString(int a)
{
    ostringstream temp;
    temp << a;
    return temp.str();
}



int main(int argc, char* argv[])
{


    ////////////
    ///CODE
    ////////////

    const int L = 9; //configure o tamanho da cadeia L = 8 limite máximo
    printfln("L = ", L);

    int N =L;
    printfln("N = ", N);

    int dim = pow(4,L);

    auto t1 = 1.0; //hopping
    printfln("t_hopp = ", t1);

    auto U0 = 2.0; //onsite inicial
    printfln("U_0 = ", U0);

    auto  V0 = 4.0; //intersite inicial
    printfln("V_0 = ",V0);

    //pinning fields (Placed at the ends of the chain )
    auto h_z = 0.0;
    printfln("pinning field = ", h_z);


    /// para fazer 111 se L=3, ou 11 se L=2,
    /// depois faço 111 binario para decimal para
    /// dar entrada no loop do soma string

    int var = 0;
    for(int i = 0 ; i <= L-1 ; i++)
    {
        var += pow(10,i);
    }

    int bin = var; //converter binario para decimal
    int dec = 0;
    for(int i = 0; bin > 0; i++)
    {
        dec = dec + pow(2, i) * (bin % 10);
        bin = bin / 10;
    }

    cout << "dec = " << dec << endl;


    string tal[dim]; ///criando um vetor de string
    int cot = 0;

    for(int i = 0 ; i <= dec ; i++ ) //for one
    {
        bitset<L> bn(i); //convert decimal i em binario escrito em L digitos
        //cout << "bn" << i << ") = " << bn <<endl;
        string soma_string1;
        for(int k = L-1; k >= 0; k--)
        {
            string bs = IntToString(bn[k]); //já convert bn[0] em string
            //cout << "bs " << bs << endl;
            soma_string1 += bs;
        }
        //cout << "soma_string1 = " << soma_string1 << endl;

        for(int j = 0 ; j <= dec ; j++)
        {
            bitset<L> bn(j); //O termo dentro de < > deve ser configurado para <N>
            string soma_string2;

            for(int k = L-1; k >= 0; k--)
            {
                string bs = IntToString(bn[k]); //ja convert bn[0] em string
                //cout << "bs " << bs << endl;
                soma_string2 += bs;
            }

            //cout << "soma_string2 = " << soma_string2 << endl;
            //cout << "juntando = " << soma_string1 + soma_string2 << endl;
            tal[cot] += soma_string1 + soma_string2;
            cot += 1;
        }

    } //fecha loop for one


    //for(int i = 0; i <= dim-1  ; i++) //para conferir se armazenou em tal
    //{
    //    cout <<"soma_string1 + soma_string2 = " << tal[i] << endl;
    //}


    ///////////////////////////////////////////////
    /// Para contar quantos vetor ha na base L=N
    ////////////////////////////////////////////////

    int total_vet = 0;
    for(int i = 0; i <= dim-1  ; i++) //para contar num vetores da base de L=N
    {
        int st = 0;
        for(int j = 0 ; j < 2*L ; j++ )
        {

            string rec = tal[i].substr(j,1); //s.substr(posição,comprimento); Para extracao
            int cl = stoi(rec); //stoi convert string em int
//          //cout << "resultado = "<< cl << endl;
            st += cl;
        }
        //cout << "st = " << st << endl;
        if(st == N) // condicao para L=N
        {
            //cout << "tal[i] = " << tal[i] << endl;
            total_vet += 1;
        }
    }

    cout << "vet_base " << total_vet << endl;
    cout << "Num vetores da base = " << total_vet << "\n\n";






    ///////////////////////////////////////////////////////////
    /// Configuracao preenchimento utilizando condicoes acima
    ///////////////////////////////////////////////////////////

    auto sites = Electron(L);

    // state = InitState(sites); // coloquei fora do loop

    auto psi_base = std::vector<MPS>(total_vet);
    //auto sites_base = std::vector<MPS>(total_vet);

    int cpi = 0; //usado no loop abaixo para para rotular os vetores da base


    for(int i = 0; i <= dim-1  ; i++) //loop para configurar o preenchimento
    {
        int st = 0;
        for(int j = 0 ; j < 2*L ; j++ )
        {

            string rec = tal[i].substr(j,1); //s.substr(posição,comprimento); Para extracao
            int cl = stoi(rec); //stoi convert string em int
//          //cout << "resultadopoHA = "<< cl << endl;
            st += cl;
        }
        //cout << "st = " << st << endl;
        if(st == N) //configura o preenchimento L=N
        {
            auto state = InitState(sites);

            //cout << "tal[" << i << "] = " << tal[i] << "\n\n";
            for(int j = 0 ; j < L ; j++ ) // loop que cobre o digitos do setor Spin down
            {
                string s_1 = tal[i].substr(j,1); //s.substr(posição,comprimento); Para extracao
                int cl_dn = stoi(s_1); //stoi convert string em int

                string s_2 = tal[i].substr(j+L,1);
                int cl_up = stoi(s_2);

                //cout << "cl_1 = "<< cl_dn << endl;
                //cout << "cl_2 = "<< cl_up << endl;

                //Inicia condicoes par preenchimento
                if(cl_dn == cl_up && cl_dn != 0)
                {
                    //cout << "UpDn " << j+1 << endl;
                    state.set(j+1,"UpDn");
                }
                if(cl_dn == 1 && cl_up == 0)
                {
                    //cout << "Dn   " << j+1 << endl;
                    state.set(j+1,"Dn");
                }
                if(cl_dn == 0 && cl_up ==1)
                {
                    //cout << "Up   " << j+1 << endl;
                    state.set(j+1,"Up");
                }
                if(cl_dn == 0 && cl_up ==0)
                {
                    //cout << "Emp  " << j+1 << endl;
                    state.set(j+1,"Emp");
                }

            }// loop que cobre o digitos do setor Spin down



            auto psi0 = MPS(state);
            //sites_base[cpi] = state;
            //cout << "cpi = " << cpi << endl;
            psi_base[cpi] = psi0; //quantos vetores da base nos temos?
            cpi += 1;

        } // fecha loop que configura o preenchimento L=N

    } //// fecha loop para configurar o preenchimento

    /////////////////////////////////////////////////
    ///Conferindo se os preenchimento estão corretos
    /////////////////////////////////////////////////

    /*
    for(int i = 0 ; i < total_vet ; i++)
    {

        for(int j = 1; j <= L ; j++)
        {
            psi_base[i].position(j);
            auto ket = psi_base[i](j);
            auto bra = dag(prime(ket,"Site"));
            auto N_op = op(sites,"Ndn",j);
            auto Ntot = eltC(bra*N_op*ket).real();
            auto ndn_Evoldn = Ntot;

            psi_base[i].position(j);
            ket = psi_base[i](j);
            bra = dag(prime(ket,"Site"));
            N_op = op(sites,"Nup",j);
            Ntot = eltC(bra*N_op*ket).real();
            auto ndn_Evolup = Ntot;

            //cout <<"psi["<<i<<"]_" << "site i " << j << " dn_ " << ndn_Evoldn  << " up_ " << ndn_Evolup  << endl;
            //cout <<"psi["<<i<<"]_" << "site i " << j << " dn_ " << ndn_Evoldn  << " up_ " << ndn_Evolup  << endl;

        } //fecha for densidade up e dn
    }

    */



    ////////////////////////////////////
    ///Utilizando os vetores da base
    ///para ensaduichar o hamiltoniano
    ///e escrever o hamiltoniano numerico
    ////////////////////////////////////

    ///Ground State (GS) Hamiltonian
    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= L; ++i)
        {
        ampo += U0,"Nupdn",i;
        //ampo += 0.001,"Nupdn",i;
        }
    for(int b = 1; b < L; ++b)
        {
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;
        ampo += V0,"Ntot",b,"Ntot",b+1;
        }
        //ampo += h_z, "Sz", 1;
        //ampo += h_z, "Sz", L;

    auto H0 = toMPO(ampo);


    //printando a primeira linha

    for(int i = 0 ; i < total_vet ; i++) //elementos linha abre
    {
         auto line = psi_base[i];

         for(int j = 0 ; j < total_vet ; j++) //elementos coluna abre
         {
             auto elem = inner(psi_base[i], H0, psi_base[j]);

             //cout << "i_(" << i << ")_j_(" << j << ")___" << elem << endl;
             //cout << elem << " ";

         } //elementos coluna fecha
         //cout << "\n";

    } //elementos linha fecha

    /*!
       Print de string dos elementos da matriz]
       vou utilizar a saida dos elementos para diagonalizar
       no python.
    */

    ///Armadillo para configurar elementos matriz
    arma::mat A; //matriz hamiltoniana
    A.set_size(total_vet,total_vet);

    for(int i = 0 ; i < total_vet ; i++) //elementos linha abre
    {
         auto line = psi_base[i];

         for(int j = 0 ; j < total_vet ; j++) //elementos coluna abre
         {
             auto elem = inner(psi_base[i], H0, psi_base[j]);

             //cout << "#* " << i << j << " " << elem << endl; //saída utilizada para comparar com python
             A(i,j)=elem;


         } //elementos coluna fecha
         //cout << "\n";

    } //elementos linha fecha


    cout << A << endl;

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);

    cout << "autovalores =  " << eigval << endl;






































return 0;
}

