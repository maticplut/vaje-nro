
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <regex>
#include <sstream>

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

        celice.push_back({ a, b ,c ,d});
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
    //preberemo točke
    PreberiTocke("primer4mreza.txt");

    //ustvarimo matriko A in vektor b
    int size = X.size();
    std::vector<std::vector<int> > A(size, std::vector<int>(size, 0));
    std::vector<double>  b(size, 0.0);

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

                for (int tockaRP = 0; tockaRP < tocke_v_RP.size(); tockaRP++) {
                    if (vozlisce == tockaRP) {
                        vrednost = vrednostiRP[st_rp];
                        vrsta = vrsteRP[st_rp];
                    }
                }
            }

            if (vrsta == 1) {
                A[vozlisce][vozlisce] = 1;
                b[vozlisce] = vrednost;
            }
            else if (vrsta == 2) {
                
                int st_sos = 0;
                for (int i = 1; i < 4; i++) {
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

}