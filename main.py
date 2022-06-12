import math
from re import T
import numpy as np
import casadi as cs

from casadi import *
from numpy import angle, array
from scipy import linalg as lin
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from hybrid_kinematics import kinematics


def close_event():
  plt.close()

def main(ang, t, ht_shel, deb_x, deb_y):
  run = kinematics(ang, t, ht_shel, deb_x, deb_y)
  link = run.fk(np.deg2rad(57.0), ang, 5.0)

  xcoords = [0, link[0], link[0] + link[1], link[0] + link[1] + link[2]]
  ycoords = [0, link[10], link[10], link[10]]
  zcoords = [0, link[3], link[3] - link[4], link[3] - link[4] - link[5]]

  xcoords_water = [link[0] + link[1] + link[2], link[0] + link[1] + link[2] + link[6]]
  ycoords_water = [link[10],link[10]]
  zcoords_water = [link[3] - link[4] - link[5], link[3] - link[4] - link[5] - link[7]]

  x_val = [0, link[0], link[0] + link[1], link[0] + link[1] + link[2], link[0] + link[1] + link[2] + link[6]]
  y_val = [0, link[10], link[10], link[10], link[10]]
  z_val = [0, link[3], link[3] - link[4], link[3] - link[4] - link[5], link[3] - link[4] - link[5] - link[7]]

  fig = plt.figure()
  timer = fig.canvas.new_timer(interval=1000)
  timer.add_callback(close_event)
  ax = fig.gca(projection='3d')
  ax.scatter3D(0, 0, 0, cmap='Greens', c='r')
  ax.scatter3D(x_val, y_val, z_val, cmap='Greens', c='r')
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')

  ax.plot3D(x_val, y_val, z_val)
  ax.plot3D(xcoords, ycoords, zcoords, '-', c='gray')
  ax.plot3D(xcoords_water, ycoords_water, zcoords_water, '-', c='blue')
  timer.start()
  plt.show()

if __name__ == "__main__":
  t = np.linspace(0,5,1000)
  ang = np.deg2rad(5)
  ht_shel = 2.75
  # ht_shelt = 5
  # ht_shel = np.linspace(2.75,5,10)
  deb_x = np.linspace(3.0, 20.0, 10)
  deb_y = np.linspace(1.1, 18.1, 10)
  # deb_x = 3.0
  # deb_y = 1.1
  ht = 5.0
  # print(range(len(ht_shel)-1))
  for i in range(len(deb_x)):
    main(ang, t, ht_shel, deb_x[i], deb_y[i])