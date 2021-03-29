#ifndef OKADAPOINTSOURCESURFACE_H_
#define OKADAPOINTSOURCESURFACE_H_

double cos_ (double angle);
double sin_ (double angle);

double strike_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double strike_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double strike_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);
double strike_disp_z_nocall (double x, double y, double z, double a, double b, double c, double theta, double delt, double M);

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
