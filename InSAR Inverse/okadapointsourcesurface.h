#ifndef OKADAPOINTSOURCESURFACE_H_
#define OKADAPOINTSOURCESURFACE_H_

#include <cmath>

const double M_PI = 3.14159265358979323846;

const double lame_mu = 3.0e10; //Pa
const double lame_lambda = 3.0e10; //Pa
const double alpha = (lame_lambda + lame_mu)/(lame_lambda + 2.0*lame_mu); //unitless

double cos_ (double angle);
double sin_ (double angle);

double strike_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double strike_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double strike_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);

double dip_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double dip_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double dip_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);

double tensile_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double tensile_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double tensile_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);

double inflation_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double inflation_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double inflation_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);

#endif
