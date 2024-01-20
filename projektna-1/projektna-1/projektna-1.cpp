
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <regex>
#include <sstream>
#include <omp.h>
#include <chrono>

std::vector<double> X; // vektor X koordinat
std::vector<double> Y; // vektor Y koordinat
std::vector<std::vector<double>> celice; // vektor celic

std::vector<int> vrsteRP;
std::vector<int> vrednostiRP;
std::vector<std::vector<int>> tockeRP;

void PreberiTocke(std::string filename) {

    std::ifstream indata;
    indata.open(filename);

    std::vector<int> velikosti;

    int zaporedna;
    std::string text;
    double a, b, c, d;
    char semicolon, colon;

    //TOČKE

    //preberemo št. točk


    std::string vrstica;
    std::getline(indata, vrstica);
    std::stringstream tt(vrstica);

    tt >> text >> a;
    velikosti.push_back(a);

    tt.clear();

    //std::cout << "Stevilo tock: " << velikosti[0] << "\n";

    //gremo čez vse točke, shranmo v vektor
    for (int i = 0; i < velikosti[0]; i++) {

        std::string line;
        std::getline(indata, line);
        std::stringstream ss(line);

        ss >> zaporedna >> semicolon >> a >> colon >> b;

        X.push_back(a);
        Y.push_back(b);
    }

    //CELICE

    //preberemo št. celic
    std::getline(indata, vrstica);
    std::getline(indata, vrstica);

    std::stringstream cc(vrstica);

    cc >> text >> a;
    velikosti.push_back(a);

    cc.clear();

    //std::cout << "Stevilo celic: " << velikosti[1] << std::endl;

    //gremo čez vse celice
    for (int i = 0; i < velikosti[1]; i++) {

        std::string line;
        std::getline(indata, line);
        std::stringstream ss(line);

        ss >> zaporedna >> semicolon >> a >> colon >> b >> colon >> c >> colon >> d;

        celice.push_back({ a, b ,c ,d });
    }

    //ROBNI POGOJI

    //preberemo št. robnih pogojev

    std::getline(indata, vrstica);
    std::getline(indata, vrstica);

    std::stringstream rp(vrstica);

    rp >> text >> text >> a;
    velikosti.push_back(a);

    cc.clear();

    //std::cout << "Stevilo robnih pogojev: " << velikosti[2] << std::endl;

    //Gremo čez vse robne pogoje

    std::string uporabntext;

    for (int i = 0; i < velikosti[2]; i++) {

        int stTock;
        std::string stevilka = "";
        std::vector<int> trenutneTocke;

        //vrsta RP

        std::getline(indata, vrstica);
        for (int ii = 0; ii < vrstica.length(); ii++) {
            stevilka += vrstica[ii];

            if (vrstica[ii] == ' ') {
                stevilka = "";
            }
            else if (ii == vrstica.length() - 1) {

                if (stevilka == "temperatura") {
                    vrsteRP.push_back(1);
                }
                else if (stevilka == "tok") {
                    vrsteRP.push_back(2);
                }

                stevilka = "";
            }
        }

        //vrednsot RP
        std::getline(indata, vrstica);

        for (int ii = 0; ii < vrstica.length(); ii++) {
            stevilka += vrstica[ii];

            if (vrstica[ii] == ' ') {
                stevilka = "";
            }
            else if (ii == vrstica.length() - 1) {
                int st;
                st = std::stoi(stevilka);
                vrednostiRP.push_back(st);
                stevilka = "";
            }
        }

        //gremo čez vse točke RP

        //št. točk

        std::getline(indata, vrstica);
        std::stringstream st(vrstica);

        st >> stTock;

        for (int i = 0; i < stTock; i++) {

            std::string line;
            std::getline(indata, line);
            std::stringstream ss(line);

            ss >> a;

            trenutneTocke.push_back(a);
        }
        tockeRP.push_back(trenutneTocke);

        std::getline(indata, vrstica);

    }
}

int preveri(double x_trenutno, double y_trenutno, double x_sosed, double y_sosed) {
    //preveri kje se vozlisce nahaja
    int zap_sos;

    if (abs(x_trenutno - x_sosed) < 1e-9) {
        if ((y_trenutno - y_sosed) > 0) {
            zap_sos = 1;
        }
        else if ((y_trenutno - y_sosed) < 0) {
            zap_sos = 3;
        }
    }

    else if (abs(y_trenutno - y_sosed) < 1e-9) {
        if ((x_trenutno - x_sosed) > 0) {
            zap_sos = 0;
        }
        else if ((x_trenutno - x_sosed) < 0) {
            zap_sos = 2;
        }
    }
    else {
        zap_sos = -1;
    }

    return zap_sos;
}

int main()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    //preberemo točke
    PreberiTocke("primer4mreza.txt");

    //ustvarimo matriko A in vektor b
    int size = X.size();
    int n_celice = celice.size();
    std::vector<std::vector<int> > A(size, std::vector<int>(size, 0));
    std::vector<double>  b(size, 0.0);

    
    //ustvraimo vektor rešitev z začetnimi vrednsotmi 100°C.
    std::vector<double>  T(size, 100.0);

    double dx = 0.25;
    double dy = 0.25;


    //določanje sosedov
    std::vector<std::vector<int>> sosedne_tocke;

    for (int vozlisce = 0; vozlisce < size; vozlisce++) {

        std::vector<int> sosedje = { -1, -1, -1, -1 };
        //std::cout << "vozlicse" << vozlisce << std::endl;

        for (int celica = 0; celica < celice.size(); celica++) {

            std::vector<double> trenutna_celica = celice[celica];

            if (vozlisce == trenutna_celica[0] || vozlisce == trenutna_celica[1] ||
                vozlisce == trenutna_celica[2] || vozlisce == trenutna_celica[3]) {

                for (int i = 0; i < 4; i++) {
                    int sos_vozlisce = trenutna_celica[i];
                    if (sos_vozlisce != vozlisce) {

                        double x_voz = X[vozlisce];
                        double y_voz = Y[vozlisce];
                        double x_sos = X[sos_vozlisce];
                        double y_sos = Y[sos_vozlisce];

                        int pozicija = preveri(x_voz, y_voz, x_sos, y_sos);

                        if (pozicija != -1) {

                            sosedje[pozicija] = sos_vozlisce;

                        }
                    }
                }
            }
        }

        sosedne_tocke.push_back(sosedje);

    }
    
    //Ustvarjanje matrike A in vektorja b

    for (int vozlisce = 0; vozlisce < size; vozlisce++) {


        std::vector<int> sosedi = sosedne_tocke[vozlisce];
        //definiramo sosede
        int levi_sosed = sosedi[0];
        int spodnji_sosed = sosedi[1];
        int desni_sosed = sosedi[2];
        int zgornji_sosed = sosedi[3];

        if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1) { //če je točka v sredini
            A[vozlisce][levi_sosed] = 1;
            A[vozlisce][spodnji_sosed] = 1;
            A[vozlisce][desni_sosed] = 1;
            A[vozlisce][zgornji_sosed] = 1;
            A[vozlisce][vozlisce] = -4;
        }
        else { // če je točka na robu

            int vrednost = 0;
            int vrsta = 0;


            for (int st_rp = 0; st_rp < vrsteRP.size(); st_rp++) {

                std::vector<int> tocke_v_RP = tockeRP[st_rp];

                for (int i = 0; i < tocke_v_RP.size(); i++) {
                    int tockaRP = tocke_v_RP[i];
                    if (vozlisce == tockaRP) {
                        vrednost = vrednostiRP[st_rp];
                        vrsta = vrsteRP[st_rp];
                        //std::cout << "tu sm 2" << vozlisce << std::endl;
                    }
                }
            }
            //std::cout << vrednost << " " << vrsta << "\n";

            if (vrsta == 1) {
                A[vozlisce][vozlisce] = 1;
                b[vozlisce] = vrednost;
            }
            else if (vrsta == 2) {
                
                int st_sos = 0;
                for (int i = 0; i < 4; i++) {
                    if (sosedi[i] != -1) {
                        st_sos += 1;
                    }
                }
                if (st_sos == 3) {

                    if (levi_sosed == -1) {
                        A[vozlisce][vozlisce] += -4;
                        A[vozlisce][desni_sosed] += 2;
                        A[vozlisce][spodnji_sosed] += 1;
                        A[vozlisce][zgornji_sosed] += 1;
                        b[vozlisce] = -2 * (vrednost * dx) / 4;
                    }

                    if (desni_sosed == -1) {
                        A[vozlisce][vozlisce] += -4;
                        A[vozlisce][levi_sosed] += 2;
                        A[vozlisce][spodnji_sosed] += 1;
                        A[vozlisce][zgornji_sosed] += 1;
                        b[vozlisce] = -2 * (vrednost * dx) / 4;
                    }

                    if (spodnji_sosed == -1) {
                        A[vozlisce][vozlisce] += -4;
                        A[vozlisce][zgornji_sosed] += 2;
                        A[vozlisce][levi_sosed] += 1;
                        A[vozlisce][desni_sosed] += 1;
                        b[vozlisce] = -2 * (vrednost * dy) / 4;
                    }
                    if (zgornji_sosed == -1) {
                        A[vozlisce][vozlisce] += -4;
                        A[vozlisce][spodnji_sosed] += 2;
                        A[vozlisce][levi_sosed] += 1;
                        A[vozlisce][desni_sosed] += 1;
                        b[vozlisce] = -2 * (vrednost * dy) / 4;
                    }
                }


            }
        }
    }

    //shranimo matriko in vektor v .txt
    std::ofstream matrikaID("matrika.txt");

    matrikaID << size << "\n";
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrikaID << A[i][j] << ",";
        }
        matrikaID << "\n";
    }

    matrikaID << "\n";
    matrikaID << size << "\n";
    for (int i = 0; i < size; i++)
    {
        matrikaID << b[i] << "\n";

    }

    //Iterativno reševanje sistema enačb
    int iter = 2000;

    //omp_set_num_threads(20);

    #pragma omp parallel for
    
    for (int iitt = 0; iitt < iter; iitt++)
    {
        //Gauss-Seidel metoda

        for (int jj = 0; jj < size; jj++)
        {
            double d = b[jj];

            for (int ii = 0; ii < size; ii++)
            {
                if (jj != ii) {

                    d -= (A[jj][ii] * T[ii]);

                }

                T[jj] = d / A[jj][jj];
            }

        }

    }


    std::ofstream fileID("rezultat_vtk.vtk");

    fileID << "# vtk DataFile Version 3.0\n";
    fileID << "Mesh_1\n";
    fileID << "ASCII\n";
    fileID << "DATASET UNSTRUCTURED_GRID\n";
    fileID << "POINTS " << size << " float\n";

    for (int koordinata_id = 0; koordinata_id < size; ++koordinata_id) {
        fileID << X[koordinata_id] << " " << Y[koordinata_id] << " 0\n";
    }

    fileID << "\n";
    fileID << "CELLS " << n_celice << " " << n_celice * 5 << "\n";

    for (int celica_id = 0; celica_id < n_celice; ++celica_id) {
        int vozl_id1 = celice[celica_id][0];
        int vozl_id2 = celice[celica_id][1];
        int vozl_id3 = celice[celica_id][2];
        int vozl_id4 = celice[celica_id][3];
        fileID << "4 " << vozl_id1 << " " << vozl_id2 << " " << vozl_id3 << " " << vozl_id4 << "\n";
    }

    fileID << "\n";
    fileID << "CELL_TYPES " << n_celice << "\n";

    for (int celica_id = 0; celica_id < n_celice; ++celica_id) {
        fileID << "9\n";
    }

    fileID << "\n";
    fileID << "POINT_DATA " << size << "\n";
    fileID << "SCALARS Temperature float 1\n";
    fileID << "LOOKUP_TABLE default\n";

    for (int koordinata_id = 0; koordinata_id < size; ++koordinata_id) {
        fileID << T[koordinata_id] << "\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_duration = end_time - start_time;
    std::cout << "Čas izvajanja programa " << time_duration.count() << " seconds" << std::endl;

}