import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
    
class PointOkada:

    def __init__(self, mu = 3.0e10, lam = 3.0e10):
        self.mu = mu
        self.lam = lam
        self.prefactor = self.mu/(self.mu + self.lam)
        self.coords = None
        self.dx = None
        self.dy = None
        self.dz = None
        
    def pfunc(self, y, d, delta):
        return y*np.cos(delta) + d*np.sin(delta)
        
    def qfunc(self, y, d, delta):
        return y*np.sin(delta) - d*np.cos(delta)
        
    def Rfunc(self, x, y, d):
        return np.sqrt(x**2 + y**2 + d**2)
    
    def str_ux(self, x, q, R, I, delta):
        return -1/(2*np.pi*self.mu)*(3*x**2*q/R**5 + I*np.sin(delta))
        
    def str_uy(self, x, y, q, R, I, delta):
        return -1/(2*np.pi*self.mu)*(3*x*y*q/R**5 + I*np.sin(delta))
    
    def str_uz(self, x, d, q, R, I, delta):
        return -1/(2*np.pi*self.mu)*(3*x*d*q/R**5 + I*np.sin(delta))
        
    def dip_ux(self, x, p, q, R, I, delta):
        return -1/(2*np.pi*self.mu)*(3*x*p*q/R**5 - I*np.sin(delta)*np.cos(delta))
        
    def dip_uy(self, y, p, q, R, I, delta):
        return -1/(2*np.pi*self.mu)*(3*y*p*q/R**5 - I*np.sin(delta)*np.cos(delta))
        
    def dip_uz(self, d, p, q, R, I, delta):
        return -1/(2*np.pi*self.mu)*(3*d*p*q/R**5 - I*np.sin(delta)*np.cos(delta))
        
    def I1(self, x, y, d, R):
        return self.prefactor*y*(1/(R*(R + d)**2) - x**2*(3*R + d)/(R**3*(R + d)**3))
        
    def I2(self, x, y, d, R):
        return self.prefactor*x*(1/(R*(R + d)**2) - y**2*(3*R + d)/(R**3*(R + d)**3))
        
    def I3(self, x, R, I):
        return self.prefactor*x/R**3 - I
        
    def I4(self, x, y, d, R):
        return self.prefactor*(-x*y*(2*R + d)/(R**3*(R + d)**2))
        
    def I5(self, x, d, R):
        return self.prefactor*(1/(R*(R + d)) - x**2*(2*R + d)/(R**3*(R + d)**2))
        
    def str_disp_x(self, x, y, d, delta):
        q = self.qfunc(y, d, delta)
        R = self.Rfunc(x, y, d)
        i1 = self.I1(x, y, d, R)
        return self.str_ux(x, q, R, i1, delta)
        
    def str_disp_y(self, x, y, d, delta):
        q = self.qfunc(y, d, delta)
        R = self.Rfunc(x, y, d)
        i2 = self.I2(x, y, d, R)
        return self.str_uy(x, y, q, R, i2, delta)
        
    def str_disp_z(self, x, y, d, delta):
        q = self.qfunc(y, d, delta)
        R = self.Rfunc(x, y, d)
        i4 = self.I4(x, y, d, R)
        return self.str_uz(x, d, q, R, i4, delta)
        
    def dip_disp_x(self, x, y, d, delta):
        p = self.pfunc(y, d, delta)
        q = self.qfunc(y, d, delta)
        R = self.Rfunc(x, y, d)
        i2 = self.I2(x, y, d, R)
        i3 = self.I3(x, R, i2)
        return self.dip_ux(x, p, q, R, i3, delta)
    
    def dip_disp_y(self, x, y, d, delta):
        p = self.pfunc(y, d, delta)
        q = self.qfunc(y, d, delta)
        R = self.Rfunc(x, y, d)
        i1 = self.I1(x, y, d, R)
        return self.dip_uy(y, p, q, R, i1, delta)
    
    def dip_disp_z(self, x, y, d, delta):
        p = self.pfunc(y, d, delta)
        q = self.qfunc(y, d, delta)
        R = self.Rfunc(x, y, d)
        i5 = self.I5(x, d, R)
        return self.dip_uz(d, p, q, R, i5, delta)
        
    def CalcDisplacement(self, dispfield, d, delta):
        
        for j in range(dispfield.ny):
            for i in range(dispfield.nx):
                x = dispfield.coords[j][i][0]
                y = dispfield.coords[j][i][1]
                
                dispfield.dx[j][i] = self.str_disp_x(x, y, d, delta)
                dispfield.dy[j][i] = self.str_disp_y(x, y, d, delta)
                dispfield.dz[j][i] = self.str_disp_z(x, y, d, delta)
                
        return dispfield
        
    def CalcDisplacementRotated(self, dispfield, d, delta, theta):
        
        for j in range(dispfield.ny):
            for i in range(dispfield.nx):
                x = dispfield.coords[j][i][0]
                y = dispfield.coords[j][i][1]
                
                newx = x*np.cos(theta) - y*np.sin(theta)
                newy = x*np.sin(theta) + y*np.cos(theta)
                
                dispfield.dx[j][i] = self.str_disp_x(newx, newy, d, delta)*np.cos(-theta) - self.str_disp_y(newx, newy, d, delta)*np.sin(-theta)
                dispfield.dy[j][i] = self.str_disp_x(newx, newy, d, delta)*np.sin(-theta) + self.str_disp_y(newx, newy, d, delta)*np.cos(-theta)
                dispfield.dz[j][i] = self.str_disp_z(newx, newy, d, delta)
                
        return dispfield
        
    def CalcDispRotTrans(self, dispfield, a, b, d, delta, theta):
        
        for j in range(dispfield.ny):
            for i in range(dispfield.nx):
                x = dispfield.coords[j][i][0]
                y = dispfield.coords[j][i][1]
                
                newx = (x - a)*np.cos(theta) - (y - b)*np.sin(theta)
                newy = (x - a)*np.sin(theta) + (y - b)*np.cos(theta)
                
                dispfield.dx[j][i] = self.str_disp_x(newx, newy, d, delta)*np.cos(-theta) - self.str_disp_y(newx, newy, d, delta)*np.sin(-theta)
                dispfield.dy[j][i] = self.str_disp_x(newx, newy, d, delta)*np.sin(-theta) + self.str_disp_y(newx, newy, d, delta)*np.cos(-theta)
                dispfield.dz[j][i] = self.str_disp_z(newx, newy, d, delta)
                
        return dispfield
        
    def CalcGreensRotTrans(self, x, y, a, b, d, strike, delta, rake):
        newx = (x - a)*np.cos(strike) - (y - b)*np.sin(strike)
        newy = (x - a)*np.sin(strike) + (y - b)*np.cos(strike)
        
        strdx = self.str_disp_x(newx, newy, d, delta)*np.cos(-strike) - self.str_disp_y(newx, newy, d, delta)*np.sin(-strike)
        strdy = self.str_disp_x(newx, newy, d, delta)*np.sin(-strike) + self.str_disp_y(newx, newy, d, delta)*np.cos(-strike)
        strdz = self.str_disp_z(newx, newy, d, delta)
        
        dipdx = self.dip_disp_x(newx, newy, d, delta)*np.cos(-strike) - self.dip_disp_y(newx, newy, d, delta)*np.sin(-strike)
        dipdy = self.dip_disp_x(newx, newy, d, delta)*np.sin(-strike) + self.dip_disp_y(newx, newy, d, delta)*np.cos(-strike)
        dipdz = self.dip_disp_z(newx, newy, d, delta)
        
        dx = strdx*np.cos(rake) + dipdx*np.sin(rake)
        dy = strdy*np.cos(rake) + dipdx*np.sin(rake)
        dz = strdz*np.cos(rake) + dipdz*np.sin(rake)
        
        return dx, dy, dz
        
    def CalcGreensX(self, x, y, a, b, d, strike, delta, rake):
        newx = (x - a)*np.cos(strike) - (y - b)*np.sin(strike)
        newy = (x - a)*np.sin(strike) + (y - b)*np.cos(strike)
        
        strdx = self.str_disp_x(newx, newy, d, delta)*np.cos(-strike) - self.str_disp_y(newx, newy, d, delta)*np.sin(-strike)
        dipdx = self.dip_disp_x(newx, newy, d, delta)*np.cos(-strike) - self.dip_disp_y(newx, newy, d, delta)*np.sin(-strike)
        
        dx = strdx*np.cos(rake) + dipdx*np.sin(rake)
        
        return dx

    def CalcGreensY(self, x, y, a, b, d, strike, delta, rake):
        newx = (x - a)*np.cos(strike) - (y - b)*np.sin(strike)
        newy = (x - a)*np.sin(strike) + (y - b)*np.cos(strike)
        
        strdy = self.str_disp_x(newx, newy, d, delta)*np.sin(-strike) + self.str_disp_y(newx, newy, d, delta)*np.cos(-strike)
        dipdy = self.dip_disp_x(newx, newy, d, delta)*np.sin(-strike) + self.dip_disp_y(newx, newy, d, delta)*np.cos(-strike)
        
        dy = strdy*np.cos(rake) + dipdy*np.sin(rake)
        
        return dy
        
    def CalcGreensZ(self, x, y, a, b, d, strike, delta, rake):
        newx = (x - a)*np.cos(strike) - (y - b)*np.sin(strike)
        newy = (x - a)*np.sin(strike) + (y - b)*np.cos(strike)
        
        strdz = self.str_disp_z(newx, newy, d, delta)
        dipdz = self.dip_disp_z(newx, newy, d, delta)
        
        dz = strdz*np.cos(rake) + dipdz*np.sin(rake)
        
        return dz
        
if __name__ == "__main__":

    x = 2
    y = 3
    d = 4
    delta = 70.0/180.0*np.pi
    
    test = PointOkada()
    strx = test.str_disp_x(x, y, d, delta)*test.mu
    stry = test.str_disp_y(x, y, d, delta)*test.mu
    strz = test.str_disp_z(x, y, d, delta)*test.mu
    
    print(strx, stry, strz)
    
    dipx = test.dip_disp_x(x, y, d, delta)*test.mu
    dipy = test.dip_disp_y(x, y, d, delta)*test.mu
    dipz = test.dip_disp_z(x, y, d, delta)*test.mu
    
    print(dipx, dipy, dipz)
