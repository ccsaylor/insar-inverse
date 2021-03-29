#define _USE_MATH_DEFINES

#include <cmath>
//#include <array>
#include "okadapointsourcesurface.h"
#include "okadarectangularsourcesurface.h"

extern double lame_mu, lame_lambda;

double xs (double x, double y, double a, double b, double theta) {
	return (x - a)*cos_(theta) - (y - b)*sin_(theta);
}

double ys (double x, double y, double a, double b, double theta) {
	return (x - a)*sin_(theta) + (y - b)*cos_(theta);
}

double ps (double x, double y, double a, double b, double d, double theta, double delt) {
	return ys(x, y, a, b, theta)*cos_(delt) + d*sin_(delt);
}

double qs (double x, double y, double a, double b, double d, double theta, double delt) {
	return ys(x, y, a, b, theta)*sin_(delt) - d*cos_(delt);
}

double y_tilde (double x, double y, double a, double b, double d, double theta, double delt, double eta) {
	return eta*cos_(delt) + qs(x, y, a, b, d, theta, delt)*sin_(delt);
}

double d_tilde (double x, double y, double a, double b, double d, double theta, double delt, double eta) {
	return eta*sin_(delt) - qs(x, y, a, b, d, theta, delt)*cos_(delt);
}

double Rs (double x, double y, double a, double b, double d, double theta, double delt, double eta, double xi) {
	return sqrt(pow(xi, 2.0) + pow(eta, 2.0) + pow(qs(x, y, a, b, d, theta, delt), 2.0));
}

double I4 (double x, double y, double a, double b, double d, double theta, double delt, double eta, double xi) {
	if (cos_(delt) == 0) {
		return -lame_mu/(lame_lambda + lame_mu)*qs(x, y, a, b, d, theta, delt)/(Rs(x, y, a, b, d, theta, delt, eta, xi) + d_tilde(x, y, a, b, d, theta, delt, eta));
	}
	else {
		if (Rs(x, y, a, b, d, theta, delt, eta, xi) + eta == 0) {
			return lame_mu/(lame_lambda + lame_mu)/cos_(delt)*(log(Rs(x, y, a, b, d, theta, delt, eta, xi) + d_tilde(x, y, a, b, d, theta, delt, eta)) + sin_(delt)*log(Rs(x, y, a, b, d, theta, delt, eta, xi) - eta));
		}
		else {
			return lame_mu/(lame_lambda + lame_mu)/cos_(delt)*(log(Rs(x, y, a, b, d, theta, delt, eta, xi) + d_tilde(x, y, a, b, d, theta, delt, eta)) - sin_(delt)*log(Rs(x, y, a, b, d, theta, delt, eta, xi) + eta));
		}
	}
}

double u_z_ss (double x, double y, double a, double b, double d, double theta, double delt, double eta, double xi, double U1) {
	if (Rs(x, y, a, b, d, theta, delt, eta, xi) + eta == 0) {
		return I4(x, y, a, b, d, theta, delt, eta, xi)*sin_(delt);
	}
	else {
		return -U1/(2.0*M_PI)*(d_tilde(x, y, a, b, d, theta, delt, eta)*qs(x, y, a, b, d, theta, delt)/(Rs(x, y, a, b, d, theta, delt, eta, xi)*(Rs(x, y, a, b, d, theta, delt, eta, xi) + eta)) + qs(x, y, a, b, d, theta, delt)*sin_(delt)/(Rs(x, y, a, b, d, theta, delt, eta, xi) + eta) + I4(x, y, a, b, d, theta, delt, eta, xi)*sin_(delt));
	}
}

double disp_z_ss (double x, double y, double a, double b, double d, double theta, double delt, double U1, double L, double W) {
	return u_z_ss(x, y, a, b, d, theta, delt, ps(x, y, a, b, d, theta, delt), xs(x, y, a, b, theta), U1) - u_z_ss(x, y, a, b, d, theta, delt, ps(x, y, a, b, d, theta, delt) - W, xs(x, y, a, b, theta), U1) - u_z_ss(x, y, a, b, d, theta, delt, ps(x, y, a, b, d, theta, delt), xs(x, y, a, b, theta) - L, U1) + u_z_ss(x, y, a, b, d, theta, delt, ps(x, y, a, b, d, theta, delt) - W, xs(x, y, a, b, theta) - L, U1);
}
