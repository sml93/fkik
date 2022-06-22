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

from hybrid_kinematics_force import kinematics_force


def close_event():
  plt.close()

def main(ang, t, ht_shel, deb_x, deb_y):
  run = kinematics_force(ang, t, ht_shel, deb_x, deb_y)
  link = run.fk(np.deg2rad(45.0), ang, 5.5)
  print('FK_xpoc: ', link[8])
  print('FK_ypoc: ', link[10])
  print('FK_zpoc: ', link[9])

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
  timer = fig.canvas.new_timer(interval=int)
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
  ## timer.start()
  ## plt.show()


  '''Initializing the IK method and parsing the solution into the FK to check'''
  inverse_sol = run.ik(ht_limit, deb_x, deb_y, ht_shel, vel_out)
  link = run.fkik(inverse_sol[0][0][0], inverse_sol[1][0][0], np.linalg.norm(inverse_sol[2][0][0]), inverse_sol[3][0][0])
  print('IK_xpoc: ', link[8])
  print('IK_ypoc: ', link[10])
  print('IK_zpoc: ', link[9])
  print()
  print()

  xcoords = [0, link[0], link[0] + link[1], link[0] + link[1] + link[2]]
  ycoords = [0, link[10], link[10], link[10]]
  zcoords = [0, link[3], link[3] - link[4], link[3] - link[4] - link[5]]

  xcoords_water = [link[0] + link[1] + link[2], link[0] + link[1] + link[2] + link[6]]
  ycoords_water = [link[10],link[10]]
  zcoords_water = [link[3] - link[4] - link[5], link[3] - link[4] - link[5] - link[7]]

  xIK_val = [0, link[0], link[0] + link[1], link[0] + link[1] + link[2], link[0] + link[1] + link[2] + link[6]]
  yIK_val = [0, link[10], link[10], link[10], link[10]]
  zIK_val = [0, link[3], link[3] - link[4], link[3] - link[4] - link[5], link[3] - link[4] - link[5] - link[7]]

  fig = plt.figure()
  timer = fig.canvas.new_timer(interval=int)
  timer.add_callback(close_event)
  ax = fig.gca(projection='3d')
  ax.scatter3D(0, 0, 0, cmap='Greens', c='k')
  ax.scatter3D(xIK_val, yIK_val, zIK_val, cmap='Greens', c='k')
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')

  ax.plot3D(xIK_val, yIK_val, zIK_val)
  ax.plot3D(xcoords, ycoords, zcoords, '-', c='gray')
  ax.plot3D(xcoords_water, ycoords_water, zcoords_water, '-', c='blue')
  plt.ylim([-deb_y-0.2, deb_y+0.2])
  timer.start()
  # plt.show()


if __name__ == "__main__":
  int = 8000
  t = np.linspace(0,5,1000)

  # """ Once """
  # ang = np.deg2rad(20)
  # ht_shel = 2.75            # m
  # vel_out = 152             # m/s

  # ht_limit = 6.0            # m
  # deb_x = 3.5               # m
  # deb_y = 0.5               # m

  # main(ang, t, ht_shel, deb_x, deb_y)


  """ Multiple """
  ang = np.deg2rad(20)     # rad
  ht_shel = 2.75            # m
  vel_out = 152             # m/s

  ht_limit = 6.0            # m
  deb_x = np.linspace(3.5, 20.0, 10)
  deb_y = np.linspace(0.5, 18.0, 10)

  for i in range(len(deb_x)):
    main(ang, t, ht_shel, deb_x[i], deb_y[i])