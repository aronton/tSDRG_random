#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <uni10.hpp>
#include "../MPO/mpo.h"
#include "../tSDRG_tools/tSDRG_tools.h"
#include "../tSDRG_tools/measure.h"

using namespace std;

void tSDRG_XXZ(int L, int chi, int Pdis, double Jdis, string BC, double S, double Jz, double Dim, double h, int Jseed, string checkOrNot)
{
    random_device rd;             // non-deterministic generator.
    mt19937 genRandom(rd() );    // use mersenne twister and seed is rd.
    mt19937 genFixed(Jseed);     // use mersenne twister and seed is fixed!

    uniform_real_distribution<double> Dist_J(nextafter(0.0, 1.0), 1.0); // probability distribution of J rand(0^+ ~ 1)
     
    /// create coupling list and MPO chain for OBC or PBC
    vector<double> J_list;
    if (BC == "PBC")
    {
        for(int i=0; i<L; i++)
        {
            double jvar = Dist_J(genRandom);
            jvar = (1 + Dim * pow(-1,i+1)) * Distribution_Random_Variable(Pdis, jvar, Jdis);
            J_list.push_back(jvar);
        }
    }    
    else if (BC == "OBC")
    {
        for(int i=0; i<L-1; i++)
        {
            double jvar = Dist_J(genRandom);
            jvar = (1 + Dim * pow(-1,i+1)) * Distribution_Random_Variable(Pdis, jvar, Jdis);
            J_list.push_back(jvar);
        }

    }
    // for (int i=0; i<J_list.size(); i++)
    //     cout << "J_list : " << J_list[i] << endl;
    // cout << "J_list : " << J_list << endl;

    /// STEP.1: Decompose the Hamiltonian into MPO blocks
    vector<MPO> MPO_chain;
    MPO_chain = generate_MPO_chain(L, "XXZ_" + BC, S, J_list, Jz, h);

    /// create folder in order to save data
    string file, dis, dim, folder, file_name, file_name1, file_name2, file_nameS;

    // folder = /home/aronton/tSDRG_project/tSDRG/old

    if (Jdis < 1.0 && Jdis < 0.1)
        dis = "00" + to_string( (int)(round(Jdis*100)) );
    else if (Jdis < 1.0 && Jdis >= 0.1)
        dis = "0" +  to_string( (int)(round(Jdis*100)) ); 
    else if (Jdis >= 1.0)
        dis = to_string( (int)(round(Jdis*100)) );
    
    if ( Dim >= 0 )//for positive dimmeriaztion
    {
        if (Dim < 1.0 && Dim < 0.1)
            dim = "00" + to_string( (int)(round(Dim*100)) );
        else if (Dim < 1.0 && Dim >= 0.1)
            dim = "0" + to_string( (int)(round(Dim*100)) );
        else if (Dim >= 1.0)
            dim = to_string( (int)(round(Dim*100)) );
    }
    else //for negative dimmeriaztion
    {
        Dim = abs(Dim);
        if (Dim < 1.0 && Dim < 0.1)
            dim = "N00" + to_string( (int)(round(Dim*100)) );
        else if (Dim < 1.0 && Dim >= 0.1)
            dim = "N0" + to_string( (int)(round(Dim*100)) );
        else if (Dim >= 1.0)
            dim = "N"+to_string( (int)(round(Dim*100)) );
    }


    /// return: TTN(w_up and w_loc) and energy spectrum of top tensor
    vector<uni10::UniTensor<double> > w_up;      // w_up is isometry tensor = VT
    vector<int> w_loc;                           // location of each w
    vector<double> En = J_list;                  // J_list will earse to one, and return ground energy.
    bool info = 1;                               // True; if tSDRG can not find non-zero gap, info return 0, and stop this random seed.
    bool save_RG_info = 0;                       // save gaps at RG stage 

    // 儲存時間用的變數
    time_t start, end;   

    // string output = "rank:" + to_string(rank) + "\n";

    if( checkOrNot == "Y" )
    {
        if(BC == "PBC")
        {
            string CheckPathZL = "data_random/" + BC + "/Jdis" + dis + "/Dim" + dim + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "/ZL.csv";
            
            ifstream CheckZL(CheckPathZL);
            if (!CheckZL.is_open())
            {
                // cerr << "Could not open the file - ";
                cout << "\nNo_ZL_Path:" << CheckPathZL << "\n";
                // cout << "No ZL\n";
            }
            else
            {
                string line = "";
                cout << "\nZL_Path:" << CheckPathZL << "\n";
                while (getline(CheckZL, line))
                {
                    cout << line << "\n";
                }
                return ;
            }
        }
        else
        {
            string CheckPathCorr1_etoe = "data_random/" + BC + "/Jdis" + dis + "/Dim" + dim + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "/corr1_etoe.csv";
            
            ifstream CheckZL(CheckPathCorr1_etoe);
            if (!CheckZL.is_open())
            {
                // cerr << "Could not open the file - ";
                cout << "\nNo_Corr1_Path:" << CheckPathCorr1_etoe << "\n";
                // cout << "No ZL\n";
            }
            else
            {
                string line = "";
                cout << "\nCorr1_Exist:" << CheckPathCorr1_etoe << "\n";
                // while (getline(CheckPathCorr1_etoe, line))
                // {
                //     cout << line << "\n";
                // }
                return ;
            }            
        }
    }


    start = time(NULL);
    // std::cout << std::asctime(std::localtime(&start)) << start << "Tree start seconds\n";
    std::cout << "Tree start\n" << "seconds:" << start << "\ndate:" << std::asctime(std::localtime(&start)) << "\n";


    // printf("start = %f\n", start);

    tSDRG(MPO_chain, En, w_up, w_loc, chi, dis, Pdis, Jseed, save_RG_info, info);

    end = time(NULL);
    std::cout << "Tree end\n" << "seconds:" << end << "\ndate:" << std::asctime(std::localtime(&end)) << "\n\n";

    double diff1 = difftime(end, start);
    
    printf("Tree time = %f\n\n", diff1);

    // for (int i = 1; i < w_up.size(); i++)
    // {

        // string datasj("./data/jj/");
        // datasj.append(to_string(i));
        // cout << datasj << endl;
    //     w_up[i].Load(datasj);
    // }




    /// create folder
    //folder = "data/" + BC + "/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    folder = "data_random/" + BC + "/Jdis" + dis + "/Dim" + dim + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);

    string str = "mkdir -p " + folder;

    const char *mkdir = str.c_str();
    const int dir_err = system(mkdir);

    if (dir_err == -1)
    {
        cout << "Error creating directory!" << endl;
        exit(1);
    }

    /// check info if can not RG_J
    if (info == 0)
    {
        cout << "random seed " + to_string(Jseed)  + " died (can not find non-zero gap) " << endl;
        return;
    }
    else
    {
        cout << "finish inxyz " << folder << endl;
    }

    for (int i=0; i<10; i++)
        cout << En[i] << endl;
    cout << "\n";
    //string top1 = Decision_tree(w_loc, true);

    /// create isometry of other part
    vector<uni10::UniTensor<double> > w_down;    // w_down
    uni10::UniTensor<double> kara;
    w_down.assign(L-1, kara);
    for (int i = 0; i < w_up.size(); i++)
    {
        w_down[i] = w_up[i];
        uni10::Permute(w_down[i], {-3, -1, 1}, 2, uni10::INPLACE);
    }

    /// save TTN; Don't do this. Your HDD will be fill with TTN
    /*file = folder;
    for (int i = 0; i < w_up.size(); i++)
    {
        file = folder + "/w_up" + to_string(i);
        w_up[i].Save(file);
    }*/

    /// save disorder coupling J list 
    vector<double> corr12;            // corr <S1><S2>
    file = folder + "/J_list.csv";
    ofstream fout(file);              // == fout.open(file);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file << ")";
        throw runtime_error(err.str());
    }
    for (int i = 0; i < J_list.size(); i++)
    {
        fout << setprecision(16) << J_list[i] << endl;
        corr12.push_back(Correlation_St(i, w_up, w_down, w_loc) );
    }
    fout.flush();
    fout.close();


    /// save merge order
    file = folder + "/w_loc.csv";
    fout.open(file);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file << ")";
        throw runtime_error(err.str());
    }
    for (int i = 0; i < w_loc.size(); i++)
    {
        fout << w_loc[i] << endl;
    }
    fout.flush();
    fout.close();


    start = time(NULL);
    // std::cout << std::asctime(std::localtime(&start)) << start << "Tree start seconds\n";
    std::cout << "Obsevable start\n" << "seconds:" << start << "\ndate:" << std::asctime(std::localtime(&start)) << "\n";

    /// save end to end correlation (OBC only)
    if (BC == "OBC")
    {
        double corr;
        file = folder + "/corr1_etoe.csv";
        ofstream fout(file);

        if (!fout)
        {
            ostringstream err;
            err << "Error: Fail to save (maybe need mkdir " << file << ")";
            throw runtime_error(err.str());
        }
        int site1 = 0;
        int site2 = L-1;

        corr = Correlation_StSt(site1, site2, w_up, w_down, w_loc);

        fout << "x2-x1,corr"  << endl;
        fout << setprecision(16) << site2-site1 << "," << corr << endl;

        fout.flush();
        fout.close();
        cout << "OBC_corr_etoe" << endl;
        cout << corr << endl;
    }
    else
    {
        // start = time(NULL);
        // std::cout << std::asctime(std::localtime(&start)) << start << " start seconds since the Epoch\n";

        /// bulk correlation and string order parameter print
        double corr, corr1 ,corr2, sop;
        file_name1 = folder + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_corr1.csv";
        file_name2 = folder + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_corr2.csv";
        file_nameS = folder + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_string.csv";
        ofstream fout1(file_name1);
        ofstream fout2(file_name2);
        ofstream foutS(file_nameS);
        fout1 << "x1,x2,corr" << endl;
        fout2 << "x1,x2,corr" << endl;
        foutS << "x1,x2,corr" << endl;
        int site1;
        int site2;
        int r;
        /// TODO: for OBC. Now is good for PBC.
        for (site1 = 0; site1 < L/2; site1 += 1)
        {
            site2 = site1 + L/2;

            corr = Correlation_StSt(site1, site2, w_up, w_down, w_loc);
            corr1 = corr12[site1];
            corr2 = corr12[site2];  

            fout1 << setprecision(16) << site1 << "," << site2 << "," << corr << endl;
            fout2 << setprecision(16) << site1 << "," << site2 << "," << corr - corr1*corr2 << endl;

            sop = Correlation_String(site1, site2, w_up, w_down, w_loc);
            foutS << setprecision(16) << site1 << "," << site2 << "," << sop << endl;


            // if (site1 == 0)
            // {
            //     for (site2 = site1 + 1; site2 <= site1 + L/2; site2 += 1)
            //     {
            //         corr = Correlation_StSt(site1, site2, w_up, w_down, w_loc);
            //         corr1 = corr12[site1];
            //         corr2 = corr12[site2];  

            //         fout1 << setprecision(16) << site1 << "," << site2 << "," << corr << endl;
            //         fout2 << setprecision(16) << site1 << "," << site2 << "," << corr - corr1*corr2 << endl;                      
            //     }
            // }
            // else
            // {
            //     for (site2 = site1 + L/2 -10; site2 <= site1 + L/2; site2 += 1)
            //     {
            //         corr = Correlation_StSt(site1, site2, w_up, w_down, w_loc);
            //         corr1 = corr12[site1];
            //         corr2 = corr12[site2];  

            //         fout1 << setprecision(16) << site1 << "," << site2 << "," << corr << endl;
            //         fout2 << setprecision(16) << site1 << "," << site2 << "," << corr - corr1*corr2 << endl;
                    
            //         if (site2 == site1 + L/2)
            //         {
            //             // cout << "string" << time(NULL) << endl;
            //             sop = Correlation_String(site1, site2, w_up, w_down, w_loc);
            //             foutS << setprecision(16) << site1 << "," << site2 << "," << sop << endl;
            //             // cout << "string" << time(NULL) << endl;
            //         }                      
            //     }
            // }

            // for (site2 = site1 + L/2-10; site2 < site1 + L/2; site2 += 1)
            // {
            //     r = site2 - site1;
            //     //cout << r << endl;
            //     if (r <= L/2)
            //     {
            //         //cout << "TEST: " << site1 << " & " << site2 << endl;
            //         corr = Correlation_StSt(site1, site2, w_up, w_down, w_loc);
            //         corr1 = corr12[site1];
            //         corr2 = corr12[site2];

            //         fout1 << setprecision(16) << site1 << "," << site2 << "," << corr << endl;
            //         fout2 << setprecision(16) << site1 << "," << site2 << "," << corr - corr1*corr2 << endl;

            //         if (r == L/2)
            //         {
            //             // cout << "string" << time(NULL) << endl;
            //             sop = Correlation_String(site1, site2, w_up, w_down, w_loc);
            //             foutS << setprecision(16) << site1 << "," << site2 << "," << sop << endl;
            //             // cout << "string" << time(NULL) << endl;
            //         }
            //     }
            // }
        }
        fout1.flush();
        fout2.flush();
        foutS.flush();
        fout1.close();
        fout2.close();
        foutS.close();

        // end = time(NULL);
        // std::cout << std::asctime(std::localtime(&end)) << end << " end seconds since the Epoch\n";

        // double diff2 = difftime(end, start);
        // printf("Time2 = %f\n", diff2);
        
        /// save twist order parameter
        // if (BC == "PBC")
        // {

        string file_ZLC = folder + "/ZLC.csv";
        string file_ZLR = folder + "/ZL.csv";
        string file_ZLI = folder + "/ZLI.csv";

        uni10_complex128 ZLC;
        double ZL;
        double ZLI;

        ZLC = Correlation_ZL(w_up, w_down, w_loc);
        ZL = ZLC.real();
        ZLI = ZLC.imag();

        fout.open(file_ZLC);
        if (!fout)
        {
            ostringstream err;
            err << "Error: Fail to save (maybe need mkdir " << file << ")";
            throw runtime_error(err.str());
        }

        fout << "ZLC:" << endl;
        fout << setprecision(16) << ZLC << endl;

        cout << "ZLC:" << endl;
        cout << setprecision(16) << ZLC << endl;

        fout.flush();
        fout.close();

        fout.open(file_ZLR);
        if (!fout)
        {
            ostringstream err;
            err << "Error: Fail to save (maybe need mkdir " << file << ")";
            throw runtime_error(err.str());
        }

        fout << "ZL:" << endl;
        fout << setprecision(16) << ZL << endl;

        cout << "ZL:" << endl;
        cout << setprecision(16) << ZL << endl;

        fout.flush();
        fout.close();

        fout.open(file_ZLI);
        if (!fout)
        {
            ostringstream err;
            err << "Error: Fail to save (maybe need mkdir " << file << ")";
            throw runtime_error(err.str());
        }

        fout << "ZLI:" << endl;
        fout << setprecision(16) << ZLI << endl;

        cout << "ZLI:" << endl;
        cout << setprecision(16) << ZLI << endl;

        fout.flush();
        fout.close();



        // file = folder + "/ZL.csv";
        // fout.open(file);
        // if (!fout)
        // {
        //     ostringstream err;
        //     err << "Error: Fail to save (maybe need mkdir " << file << ")";
        //     throw runtime_error(err.str());
        // }
        // uni10_complex128 ZL;
        // ZL = Correlation_ZL(w_up, w_down, w_loc);
        // fout << "ZL" << endl;
        // fout << setprecision(16) << ZL << endl;

        // cout << "ZL" << endl;
        // cout << setprecision(16) << ZL << endl;

        // // cout << "w_loc" << endl;
        // // for(int i=0; i<w_loc.size(); i++ )
        // //     cout << w_loc[i] << endl;

        // fout.flush();
        // fout.close();
        // }
    }

    end = time(NULL);
    std::cout << "Observable end\n" << "seconds:" << end << "\ndate:" << std::asctime(std::localtime(&end)) << "\n";

    double diff2 = difftime(end, start);
    printf("Observable Time = %f\n\n", diff2);

    /// save ground energy
    file = folder + "/energy.csv";
    fout.open(file);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file << ")";
        throw runtime_error(err.str());
    }
    fout << "energy" << endl;
    for (int i=0; i<En.size(); i++)
    {
        fout << setprecision(16) << En[i] << endl;
    }
    fout.flush();
    fout.close();
}

void errMsg(char *arg) 
{
    cerr << "Usage: " << arg << " [options]" << endl;
    cerr << "Need 10-parameter:" << endl;
    cerr << "./job.exe <system size> <keep state of RG procedure> <Prob distribution> <disorder> <dimerization> <algo> <seed1> <seed2> <check Or Not (Y/N)>\n" << endl;
    cerr << "Example:" << endl;
    cerr << "./job.exe 32 8 10 0.1 0.1 PBC 1 1 Y(N)\n" << endl;
}

int main(int argc, char *argv[])
{
    int L;                      // system size
    int chi;                    // keep state of isometry
    string BC;                  // boundary condition
    int Pdis;                   // model of random variable disturbution
    double Jdis;                // J-coupling disorder strength
    double Dim;			        // Dimerization constant
    int seed1;                  // random seed number in order to repeat data
    int seed2;                  // random seed number in order to repeat data
    double S      = 1.5;        // spin dimension
    double Jz     = 1.0;        // XXZ model
    double h      = 0.0;        // XXZ model
    string checkOrNot;

    if (argc == 10)
    {
        stringstream(argv[1]) >> L;
        stringstream(argv[2]) >> chi;
        stringstream(argv[3]) >> Pdis;
        stringstream(argv[4]) >> Jdis;
	    stringstream(argv[5]) >> Dim;
        BC = argv[6];
        stringstream(argv[7]) >> seed1;
        stringstream(argv[8]) >> seed2;
        checkOrNot = argv[9];

        cout << "exe : " << argv[0] << endl; 
        cout << "S : " << S << endl; 
        cout << "L : " << L << endl; 
        cout << "chi : " << chi << endl; 
        cout << "Pdis : " << Pdis << endl; 
        cout << "Jdis : " << Jdis << endl;
        cout << "Dim : " << Dim << endl; 
        cout << "BC : " << BC << endl; 
        cout << "seed1 : " << seed1 << endl; 
        cout << "seed2 : " << seed2 << endl; 
        cout << "checkOrNot : " << checkOrNot << endl; 
    }
    else
    {
        errMsg(argv[0]);
        return 1;
    }

    for (int i=0;i<10;i++)
        cout << argv[i] << endl;

    // int rank, size; 
    // MPI_Init (&argc, &argv);  
    // //initialize MPI library 
    // MPI_Comm_size(MPI_COMM_WORLD, &size); 
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // int Nseed = seed2 - seed1+1;
    // int part = Nseed/size;
    // printf("part:%d", part);


    // for (int Jseed=rank*part+seed1+1; Jseed<=(rank+1)*part+seed1; Jseed++)
    // {
    //     tSDRG_XXZ(L, chi, Pdis, Jdis, BC, S, Jz, Dim, h, Jseed, checkOrNot, rank);
    // }

    // MPI_Finalize(); 


    for (int Jseed=seed1; Jseed<=seed2; Jseed++)
    {
        tSDRG_XXZ(L, chi, Pdis, Jdis, BC, S, Jz, Dim, h, Jseed, checkOrNot);
    }

    return 0;
}
