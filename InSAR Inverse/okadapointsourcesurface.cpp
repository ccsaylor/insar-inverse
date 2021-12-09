#include <cmath>
#include <iostream>
#include "okadapointsourcesurface.h"

double cos_ (double angle) {
	if (std::abs(cos(angle)) < 1e-8) {
		return 0;
	}
	else {
		return cos(angle);
	}
}

double sin_ (double angle) {
	if (std::abs(sin(angle)) < 1e-8) {
		return 0;
	}
	else {
		return sin(angle);
	}
}

double xx (double x, double y, double a, double b, double theta) {
	return (x - a)*cos_(theta) - (y - b)*sin_(theta);
}

double yy (double x, double y, double a, double b, double theta) {
		return (x - a)*sin_(theta) + (y - b)*cos_(theta);
}

double d (double z, double c) {
	return c - z;
}

double R (double x, double y, double z, double a, double b, double c, double theta) {
	return sqrt(pow(xx(x, y, a, b, theta), 2.0) + pow(yy(x, y, a, b, theta), 2.0) + pow(d(z, c), 2.0));
}

double p (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return yy(x, y, a, b, theta)*cos_(delt) + d(z, c)*sin_(delt);
}

double q (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return yy(x, y, a, b, theta)*sin_(delt) - d(z, c)*cos_(delt);
}

double s (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return p(x, y, z, a, b, c, theta, delt)*sin_(delt) + q(x, y, z, a, b, c, theta, delt)*cos_(delt);
}

double t (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return p(x, y, z, a, b, c, theta, delt)*cos_(delt) - q(x, y, z, a, b, c, theta, delt)*sin_(delt);
}

double A_3 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 3.0*pow(xx(x, y, a, b, theta), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double B_3 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 3.0*pow(yy(x, y, a, b, theta), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double C_3 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 3.0*pow(d(z, c), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double A_5 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 5.0*pow(xx(x, y, a, b, theta), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double B_5 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 5.0*pow(yy(x, y, a, b, theta), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double C_5 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 5.0*pow(d(z, c), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double A_7 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 7.0*pow(xx(x, y, a, b, theta), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double B_7 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 7.0*pow(yy(x, y, a, b, theta), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double C_7 (double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0 - 7.0*pow(d(z, c), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0);
}

double I_10(double x, double y, double z, double a, double b, double c, double theta) {
	return yy(x, y, a, b, theta)*(1.0/(R(x, y, z, a, b, c, theta)*pow(R(x, y, z, a, b, c, theta) + d(z, c), 2.0)) - pow(xx(x, y, a, b, theta), 2.0)*(3.0*R(x, y, z, a, b, c, theta) + d(z, c))/(pow(R(x, y, z, a, b, c, theta), 3.0)*pow(R(x, y, z, a, b, c, theta) + d(z, c), 3.0)));
}

double I_20(double x, double y, double z, double a, double b, double c, double theta) {
	return xx(x, y, a, b, theta)*(1.0/(R(x, y, z, a, b, c, theta)*pow(R(x, y, z, a, b, c, theta) + d(z, c), 2.0)) - pow(yy(x, y, a, b, theta), 2.0)*(3.0*R(x, y, z, a, b, c, theta) + d(z, c))/(pow(R(x, y, z, a, b, c, theta), 3.0)*pow(R(x, y, z, a, b, c, theta) + d(z, c), 3.0)));
}

double I_30(double x, double y, double z, double a, double b, double c, double theta) {
	return xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0) - I_20(x, y, z, a, b, c, theta);
}

double I_40(double x, double y, double z, double a, double b, double c, double theta) {
	return -xx(x, y, a, b, theta)*yy(x, y, a, b, theta)*(2.0*R(x, y, z, a, b, c, theta) + d(z, c))/(pow(R(x, y, z, a, b, c, theta), 3.0)*pow(R(x, y, z, a, b, c, theta) + d(z, c), 2.0));
}

double I_50(double x, double y, double z, double a, double b, double c, double theta) {
	return 1.0/(R(x, y, z, a, b, c, theta)*(R(x, y, z, a, b, c, theta) + d(z, c))) - pow(xx(x, y, a, b, theta), 2.0)*(2.0*R(x, y, z, a, b, c, theta) + d(z, c))/(pow(R(x, y, z, a, b, c, theta), 3.0)*pow(R(x, y, z, a, b, c, theta) + d(z, c), 2.0));
}

double u_ss_Ax0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/2.0*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 3.0) + alpha/2.0*3.0*pow(xx(x, y, a, b, theta), 2.0)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_ss_Ay0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/2.0*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0)*sin_(delt) + alpha/2.0*3.0*xx(x, y, a, b, theta)*yy(x, y, a, b, theta)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_ss_Az0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/2.0*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0)*cos_(delt) + alpha/2.0*3.0*xx(x, y, a, b, theta)*d(z, c)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_ss_Bx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -3.0*pow(xx(x, y, a, b, theta), 2.0)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) - (1.0 - alpha)/alpha*I_10(x, y, z, a, b, c, theta)*sin_(delt);
}

double u_ss_By0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -3.0*xx(x, y, a, b, theta)*yy(x, y, a, b, theta)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) - (1.0 - alpha)/alpha*I_20(x, y, z, a, b, c, theta)*sin_(delt);
}

double u_ss_Bz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -3.0*c*xx(x, y, a, b, theta)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) - (1.0 - alpha)/alpha*I_40(x, y, z, a, b, c, theta)*sin_(delt);
}

double u_ss_Cx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)*A_3(x, y, z, a, b, c, theta)/pow(R(x, y, z, a, b, c, theta), 3.0)*cos_(delt) + alpha*3.0*c*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0)*A_5(x, y, z, a, b, c, theta);
}

double u_ss_Cy0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)*3.0*xx(x, y, a, b, theta)*yy(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 5.0)*cos_(delt) + alpha*3.0*c*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 5.0)*(sin_(delt) - 5.0*yy(x, y, a, b, theta)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 2.0));
}

double u_ss_Cz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)*3.0*xx(x, y, a, b, theta)*yy(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 5.0)*sin_(delt) + alpha*3.0*c*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 5.0)*(cos_(delt) + 5.0*d(z, c)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 2.0));
}

double u_ds_Ax0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return alpha/2.0*3.0*xx(x, y, a, b, theta)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_ds_Ay0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/2.0*s(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 3.0) + alpha/2.0*3.0*yy(x, y, a, b, theta)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_ds_Az0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/2.0*t(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 3.0) + alpha/2.0*3.0*d(z, c)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_ds_Bx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -3.0*xx(x, y, a, b, theta)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) + (1.0 - alpha)/alpha*I_30(x, y, z, a, b, c, theta)*sin_(delt)*cos_(delt);
}

double u_ds_By0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -3.0*yy(x, y, a, b, theta)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) + (1.0 - alpha)/alpha*I_10(x, y, z, a, b, c, theta)*sin_(delt)*cos_(delt);
}

double u_ds_Bz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -3.0*c*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) + (1.0 - alpha)/alpha*I_50(x, y, z, a, b, c, theta)*sin_(delt)*cos_(delt);
}

double u_ds_Cx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)*3.0*xx(x, y, a, b, theta)*t(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) - alpha*15.0*c*xx(x, y, a, b, theta)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 7.0);
}

double u_ds_Cy0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/pow(R(x, y, z, a, b, c, theta), 3.0)*(cos_(2.0*delt) - 3.0*yy(x, y, a, b, theta)*t(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 2.0)) + alpha*3.0*c/pow(R(x, y, z, a, b, c, theta), 5.0)*(s(x, y, z, a, b, c, theta, delt) - 5.0*yy(x, y, a, b, theta)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 2.0));
}

double u_ds_Cz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)*A_3(x, y, z, a, b, c, theta)/pow(R(x, y, z, a, b, c, theta), 3.0)*sin_(delt)*cos_(delt) + alpha*3.0*c/pow(R(x, y, z, a, b, c, theta), 5.0)*(t(x, y, z, a, b, c, theta, delt) + 5.0*d(z, c)*p(x, y, z, a, b, c, theta, delt)*q(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 2.0));
}

double u_t_Ax0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/2.0*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0) - alpha/2.0*3.0*xx(x, y, a, b, theta)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_t_Ay0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/2.0*t(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 3.0) - alpha/2.0*3.0*yy(x, y, a, b, theta)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_t_Az0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/2.0*s(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 3.0) - alpha/2.0*3.0*d(z, c)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_t_Bx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return 3.0*xx(x, y, a, b, theta)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 5.0) - (1.0 - alpha)/alpha*I_30(x, y, z, a, b, c, theta)*pow(sin_(delt), 2.0);
}

double u_t_By0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return 3.0*yy(x, y, a, b, theta)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 5.0) - (1.0 - alpha)/alpha*I_10(x, y, z, a, b, c, theta)*pow(sin_(delt), 2.0);
}

double u_t_Bz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return 3.0*c*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 5.0) - (1.0 - alpha)/alpha*I_50(x, y, z, a, b, c, theta)*pow(sin_(delt), 2.0);
}

double u_t_Cx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)*3.0*xx(x, y, a, b, theta)*s(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 5.0) + alpha*15.0*c*xx(x, y, a, b, theta)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 7.0) - alpha*3.0*xx(x, y, a, b, theta)*z/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_t_Cy0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/pow(R(x, y, z, a, b, c, theta), 3.0)*(sin_(2.0*delt) - 3.0*yy(x, y, a, b, theta)*s(x, y, z, a, b, c, theta, delt)/pow(R(x, y, z, a, b, c, theta), 2.0)) + alpha*3.0*c/pow(R(x, y, z, a, b, c, theta), 5.0)*(t(x, y, z, a, b, c, theta, delt) - yy(x, y, a, b, theta) + 5.0*yy(x, y, a, b, theta)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0)) - alpha*3.0*yy(x, y, a, b, theta)*z/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_t_Cz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/pow(R(x, y, z, a, b, c, theta), 3.0)*(1.0 - A_3(x, y, z, a, b, c, theta)*pow(sin_(delt), 2.0)) - alpha*3.0*c/pow(R(x, y, z, a, b, c, theta), 5.0)*(s(x, y, z, a, b, c, theta, delt) - d(z, c) + 5.0*d(z, c)*pow(q(x, y, z, a, b, c, theta, delt), 2.0)/pow(R(x, y, z, a, b, c, theta), 2.0)) + alpha*3.0*d(z, c)*z/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_in_Ax0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/2.0*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0);
}

double u_in_Ay0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/2.0*yy(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0);
}

double u_in_Az0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return -(1.0 - alpha)/2.0*d(z, c)/pow(R(x, y, z, a, b, c, theta), 3.0);
}

double u_in_Bx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/alpha*xx(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0);
}

double u_in_By0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/alpha*yy(x, y, a, b, theta)/pow(R(x, y, z, a, b, c, theta), 3.0);
}

double u_in_Bz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)/alpha*d(z, c)/pow(R(x, y, z, a, b, c, theta), 3.0);
}

double u_in_Cx0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)*3.0*xx(x, y, a, b, theta)*d(z, c)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_in_Cy0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)*3.0*yy(x, y, a, b, theta)*d(z, c)/pow(R(x, y, z, a, b, c, theta), 5.0);
}

double u_in_Cz0 (double x, double y, double z, double a, double b, double c, double theta, double delt) {
	return (1.0 - alpha)*C_3(x, y, z, a, b, c, theta);
}

//Displacement equations in terms of seismic moment ###########################################################################

double strike_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_ss_Bx0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_ss_Ax0(x, y, z, a, b, c, theta, delt) - u_ss_Ax0(x, y, -z, a, b, c, theta, delt) + u_ss_Bx0(x, y, z, a, b, c, theta, delt) + z*u_ss_Cx0(x, y, z, a, b, c, theta, delt));
	}
}

double strike_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_ss_By0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_ss_Ay0(x, y, z, a, b, c, theta, delt) - u_ss_Ay0(x, y, -z, a, b, c, theta, delt) + u_ss_By0(x, y, z, a, b, c, theta, delt) + z*u_ss_Cy0(x, y, z, a, b, c, theta, delt));
	}
}

double strike_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_ss_Bz0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_ss_Az0(x, y, z, a, b, c, theta, delt) - u_ss_Az0(x, y, -z, a, b, c, theta, delt) + u_ss_Bz0(x, y, z, a, b, c, theta, delt) + z*u_ss_Cz0(x, y, z, a, b, c, theta, delt));
	}
}

double dip_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_ds_Bx0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_ds_Ax0(x, y, z, a, b, c, theta, delt) - u_ds_Ax0(x, y, -z, a, b, c, theta, delt) + u_ds_Bx0(x, y, z, a, b, c, theta, delt) + z*u_ds_Cx0(x, y, z, a, b, c, theta, delt));
	}
}

double dip_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_ds_By0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_ds_Ay0(x, y, z, a, b, c, theta, delt) - u_ds_Ay0(x, y, -z, a, b, c, theta, delt) + u_ds_By0(x, y, z, a, b, c, theta, delt) + z*u_ds_Cy0(x, y, z, a, b, c, theta, delt));
	}
}

double dip_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_ds_Bz0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_ds_Az0(x, y, z, a, b, c, theta, delt) - u_ds_Az0(x, y, -z, a, b, c, theta, delt) + u_ds_Bz0(x, y, z, a, b, c, theta, delt) + z*u_ds_Cz0(x, y, z, a, b, c, theta, delt));
	}
}

double tensile_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_t_Bx0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_t_Ax0(x, y, z, a, b, c, theta, delt) - u_t_Ax0(x, y, -z, a, b, c, theta, delt) + u_t_Bx0(x, y, z, a, b, c, theta, delt) + z*u_t_Cx0(x, y, z, a, b, c, theta, delt));
	}
}

double tensile_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_t_By0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_t_Ay0(x, y, z, a, b, c, theta, delt) - u_t_Ay0(x, y, -z, a, b, c, theta, delt) + u_t_By0(x, y, z, a, b, c, theta, delt) + z*u_t_Cy0(x, y, z, a, b, c, theta, delt));
	}
}

double tensile_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_t_Bz0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_t_Az0(x, y, z, a, b, c, theta, delt) - u_t_Az0(x, y, -z, a, b, c, theta, delt) + u_t_Bz0(x, y, z, a, b, c, theta, delt) + z*u_t_Cz0(x, y, z, a, b, c, theta, delt));
	}
}

double inflation_disp_x (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_in_Bx0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_in_Ax0(x, y, z, a, b, c, theta, delt) - u_in_Ax0(x, y, -z, a, b, c, theta, delt) + u_in_Bx0(x, y, z, a, b, c, theta, delt) + z*u_in_Cx0(x, y, z, a, b, c, theta, delt));
	}
}

double inflation_disp_y (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_in_By0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_in_Ay0(x, y, z, a, b, c, theta, delt) - u_in_Ay0(x, y, -z, a, b, c, theta, delt) + u_in_By0(x, y, z, a, b, c, theta, delt) + z*u_in_Cy0(x, y, z, a, b, c, theta, delt));
	}
}

double inflation_disp_z (double x, double y, double z, double a, double b, double c, double theta, double delt, double M) {
	if (z == 0) {
		return M/(2.0*M_PI*lame_mu)*u_in_Bz0(x, y, z, a, b, c, theta, delt);
	}
	else {
		return M/(2.0*M_PI*lame_mu)*(u_in_Az0(x, y, z, a, b, c, theta, delt) - u_in_Az0(x, y, -z, a, b, c, theta, delt) + u_in_Bz0(x, y, z, a, b, c, theta, delt) + z*u_in_Cz0(x, y, z, a, b, c, theta, delt));
	}
}
