#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

extern "C" {
    void CPolynomial_GenLagrangePoly(int Nord, double* x_lgl, double* x, int Nx, double* l);
    void CPolynomial_GenDLagrangePoly_lglpt(int Nord, double* x_lgl, double* x, int Nx, double* dl);
    void CPolynomial_GenLegendrePoly(int Nord, double* x, int Nx, double* P);
}

std::vector<double> Polynomial_GenLagrangePoly(int Nord, std::vector<double> x_lgl, std::vector<double> x){
    int Nx = x.size();
    std::vector<double> l(Nx * (Nord+1));
    CPolynomial_GenLagrangePoly( Nord, x_lgl.data(), x.data(), Nx, l.data() );
    return l;
}

std::vector<double> Polynomial_GenDLagrangePoly_lglpt(int Nord, std::vector<double> x_lgl, std::vector<double> x){
    int Nx = x.size();
    std::vector<double> dl(Nx * (Nord+1));
    CPolynomial_GenDLagrangePoly_lglpt( Nord, x_lgl.data(), x.data(), Nx, dl.data() );
    return dl;
}

std::vector<double> Polynomial_GenLegendrePoly(int Nord, std::vector<double> x){
    int Nx = x.size();
    std::vector<double> P(Nx * (Nord+1));
    CPolynomial_GenLegendrePoly( Nord, x.data(), Nx, P.data() );
    return P;
}