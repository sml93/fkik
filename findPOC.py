import math
import numpy as np
import casadi as cs

from casadi import *
from numpy import array
from scipy import linalg as lin
from scipy.integrate import odeint
from matplotlib import pyplot as plt


class findPOC():
  "Getting POC params"
  def __init__(self, ang, t, ht):
    self.vel_out = 146    # m/s
    self.m = 0.0005
    self.g = 9.81
    self.drag = 1.125
    self.D = 0.003
    self.C = 37
    self.c = self.C*(np.power(self.D,2))
    self.ang = ang
    self.t = t
    self.ht = ht

  def modelz(self, z, t):
    vz0 = z[0]
    dvdt = self.g - (self.c*vz0)/self.m
    dzdt = vz0
    return dvdt, dzdt

  def modelx(self, x, t):
    vx0 = x[0]
    dvdt = -(self.c*vx0)/self.m
    dxdt = vx0
    return dvdt, dxdt

  def solve(self, modelx, modelz, ang):
    "Solving ODE"
    vz0 = self.vel_out*np.cos(ang)
    vx0 = self.vel_out*np.sin(ang)
    x = [vx0, 0]
    z = [vz0, 0]

    solx = odeint(modelx, x, self.t)
    solz = odeint(modelz, z, self.t)

    "Getting all values for x"
    vx1 = solx[:, 0]
    x = solx[:, 1]

    "Getting all values for z and multiplying by -1"
    vz1 = solz[:, 0]*-1
    z = solz[:, 1]
    return x, z, vx1, vz1

  def getPOC(self, ang, t, ht):
    "Getting POC params"
    model = self.solve(self.modelx, self.modelz, ang)
    
    x = model[0]
    z = model[1]
    
    vx = model[2]
    vz = model[3]
    
    for i in range(len(z)):
      if ht >= z[i]:
        index = i
        # print(z[i-1], z[i], z[i+1])
        break

    time_taken = self.t[index]
    end_velx = vx[index]
    end_velz = vz[index]
    end_x = x[index]
    end_z = z[index]

    v_mag = np.sqrt(end_velx ** 2 + end_velz ** 2)

    beta_rad = math.acos(-end_velz / v_mag)
    beta_deg = beta_rad * (180 / np.pi)

    x_water = x[index] - x[0]
    z_water = z[0] - z[index]

    print('xwater: ', x_water)
    print('zwater: ', z_water)
    print("Time taken to hit roof:", time_taken, "s")
    print("x-distance: ", end_x, "metres")
    print("Beta: ", beta_deg, "degrees")
    print("z-distance: ", end_z, "metres")
    print("Vmag: ", v_mag)
    print("VelX: ", end_velx)

    # plt.plot(x[0:index+1], z[0:index+1])
    # plt.xlim(-0.5, 1)
    # plt.ylim(-2, 0.5)
    # plt.grid()
    # plt.show()
    return time_taken, x, z, x_water, z_water