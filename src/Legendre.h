#pragma once

#include <armadillo>

double Leg2D(double x, double y, int m, arma::vec &c);
arma::vec LegDerX(int m, arma::vec &c);
arma::vec LegDerY(int m, arma::vec &c);
