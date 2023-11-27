#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <regex>

int N;

std::pair<std::vector<double>, std::vector<double>> PreberiTocke(std::string filename) {
    std::vector<double> X;
    std::vector<double> Y;


    std::ifstream indata;
    indata.open(filename);

    std::string N_tock_str;
    std::getline(indata, N_tock_str);
    int N_tock;
    N_tock = N = std::stoi(N_tock_str);
    std::cout << "Stevilo tock: " << N_tock << "\n";

    for (int i = 0; i < N_tock; i++) {
        std::vector<double> points = {};
        double x, y;

        std::string stevilka = "";

        std::string line;
        std::getline(indata, line);

        std::vector<int> point_IDs;
        std::string point_ID = "";

        line = std::regex_replace(line, std::regex(":"), " ");

        for (int ii = 0; ii < line.length(); ii++) {
            stevilka += line[ii];

            if (line[ii] == ' ') {
                double st;
                st = std::stod(stevilka);
                points.push_back(st);
                stevilka = "";
            }
            else if (ii == line.length() - 1) {
                double st;
                st = std::stod(stevilka);
                points.push_back(st);
                stevilka = "";
            }
        }

        X.push_back(points[1]);
        Y.push_back(points[2]);
       //std::cout << X[i] << ", " << Y[i] << std::endl;

    }
    return std::make_pair(X,Y);

}

int main()
{
    auto podatki = PreberiTocke("vhodna_datoteka.txt");

    std::vector<double> X = podatki.first;
    std::vector<double> Y = podatki.second;

    double h = X[50] - X[49];
    
    std::vector<double> dx = {};

    for (int i = 0; i < N; i++) {
        
        if (i == 0) {
            dx.push_back((-3 * Y[i] + 4 * Y[i + 1] - Y[i + 2]) / (2 * h));
        }
        else if (i == N - 1) {
            dx.push_back((3 * Y[i] - 4 * Y[i - 1] + Y[i - 2]) / (2 * h));
        }
        else {
            dx.push_back((Y[i + 1] - Y[i - 1]) / (2 * h));
        }
        std::cout << dx[i] << "\n";
    }
    
    const char* fileName = "diff_poly.txt";
    
    std::ofstream outputFile(fileName);

    outputFile << N << "\n";
    for (int i = 0; i < N; ++i) {
        outputFile << X[i] << " " << dx[i] << "\n";
    }
    outputFile.close();

    return 0;
}

