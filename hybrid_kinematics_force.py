import math
from re import T
import numpy as np
import casadi as cs

from casadi import *
from numpy import angle, array
from scipy import linalg as lin
from scipy.integrate import odeint
from matplotlib import pyplot as plt

from findPOC import findPOC
from Aroi import getAroi


class kinematics_force():
  def setUp(self):
    super(kinematics_force, self).setUp()
    self.init_hyperparameters()

  def tearDown(self):
    super(kinematics_force, self).tearDown()

  def __init__(self, angle, t, ht_shel, deb_x, deb_y):
    self.ybs = np.deg2rad(73.333)
    self.lbs = 0.052276
    self.lsn = 0.12
    self.density = 995.65
    self.Aroi = float(getAroi()[1.5]['0'])

    # self.vel_out = 146    # m/s
    self.g = 9.81
    # self.m = 0.0005
    self.deb_x = deb_x
    self.deb_y = deb_y
    self.ht_shel = ht_shel
    self.angle = angle
    self.t = t

    # self.drag = 1.125    
    # self.D = 0.003
    # self.C = 37
    # self.c = self.C*(np.power(self.D,2))

    self.psi_ob = np.arctan(self.deb_y/self.deb_x)
    self.theta_b = np.deg2rad(20)

  def transform_ob(self, gamma_ob, psi_ob, lob):
    hob = [[np.cos(psi_ob), np.sin(psi_ob), 0, lob * np.cos(gamma_ob) * np.cos(psi_ob)],
            [np.sin(psi_ob), np.cos(psi_ob), 0, lob * np.cos(gamma_ob) * np.sin(psi_ob)],
            [0, 0, 1, lob * np.sin(gamma_ob) * np.cos(psi_ob)],
            [0, 0, 0, 1]]

    p_b = [[0], [0], [0], [1]]

    p_o = np.dot(hob, p_b)
    # print(p_o)
    return p_o[0], p_o[1], p_o[2]

  def transform_bs(self, theta, gamma_ob, psi_ob, lob):
      hbs = [[np.cos(theta), 0, np.sin(theta), self.lbs*np.cos(self.ybs)*np.cos(theta)],
              [0, 1, 0, lob*np.cos(gamma_ob)*np.sin(psi_ob)],
              [-np.sin(theta), 0, np.cos(theta), self.lbs*np.sin(self.ybs)*np.cos(theta)],
              [0, 0, 0, 1]]

      p_s = [[0], [0], [0], [1]]
      p_b = np.dot(hbs, p_s)

      return p_b[0], p_b[1], p_b[2]

  def transform_sn(self, theta, alpha_p, gamma_ob, psi_ob, lob):
      if alpha_p == 0:
          hsn = [[np.cos(theta), 0, np.sin(theta), self.lsn*np.cos((alpha_p-theta)-np.pi/2)*np.cos(theta)],
                  [0, 1, 0, lob * np.cos(gamma_ob) * np.sin(psi_ob)],
                  [-np.sin(theta), 0, np.cos(theta), self.lsn*np.cos(theta)],
                  [0, 0, 0, 1]]
      else:
          hsn = [[np.cos(theta), 0, np.sin(theta), self.lsn * np.cos((alpha_p-theta) - np.pi / 2) * np.cos(theta)],
                  [0, 1, 0, lob * np.cos(gamma_ob) * np.sin(psi_ob)],
                  [-np.sin(theta), 0, np.cos(theta), self.lsn * np.sin(np.pi / 2 - (alpha_p-theta)) * np.cos(theta)],
                  [0, 0, 0, 1]]

      p_n = [[0], [0], [0], [1]]
      p_s = np.dot(hsn, p_n)
      return p_s[0], p_s[1], p_s[2]


  def fk(self, y_ob, alpha_p, l_ob):
    "Getting the FK of the system, getting xob/yob/zob values"
    pos_ob = self.transform_ob(y_ob, self.psi_ob, l_ob)
    pos_bs = self.transform_bs(self.theta_b, y_ob, self.psi_ob, l_ob)
    pos_sn = self.transform_sn(self.theta_b, alpha_p, y_ob, self.psi_ob, l_ob)

    "X params"
    x_ob = pos_ob[0][0]   # 7
    x_bs = pos_bs[0][0]   # 8
    # x_sn0 = self.lsn*np.cos(alpha_p-np.pi/2)
    x_sn = pos_sn[0][0]

    "Y params"
    y_ob = pos_ob[1][0]
    y_bs = pos_bs[1][0]
    y_sn = pos_sn[1][0]
    y_poc = pos_ob[1][0]

    "Z params"
    z_ob = pos_ob[2][0]   # 2
    z_bs = pos_bs[2][0]   # 3
    z_sn = pos_sn[2][0]
    if alpha_p == 0:
      z_sn1 = self.lsn    # 4
    else:
      z_sn1 = self.lsn*np.sin(np.pi/2-alpha_p)    # 4
    
    height = z_ob - (z_bs + z_sn + self.ht_shel)
    poc = findPOC(alpha_p-self.theta_b, self.t, height)
    poc = poc.getPOC(alpha_p-self.theta_b, self.t, height)
    x_water = poc[3]      # 6
    z_water = poc[4]      # 7

    x_poc = x_ob + x_bs + x_sn + x_water    # 8
    z_poc = z_ob - (z_bs + z_sn + z_water)  # 9

    return x_ob, x_bs, x_sn, z_ob, z_bs, z_sn, x_water, z_water, x_poc, z_poc, y_ob, y_bs, y_sn, y_poc


  def fkik(self, y_ob, alpha_p, l_ob, psi_ob):
    "Getting FK of system"
    "Getting xob/yob/zob values"
    pos_body = self.transform_ob(y_ob, psi_ob, l_ob)

    "X params"
    x_ob = pos_body[0][0]
    # x_ob = l_ob*np.cos(y_ob)                    # 7
    print("x_ob: ", x_ob)
    x_bs = self.lbs*np.cos(self.ybs)*np.cos(self.theta_b)            # 8
    print("x_bs: ", x_bs)
    x_sn = self.lsn*np.cos(alpha_p-np.pi/2)*np.cos(self.theta_b)     # 9
    print("x_sn: ", x_sn)

    "Y params"
    y_ob = pos_body[1][0]
    print('y_ob: ', y_ob)
    y_bs = pos_body[1][0]
    y_sn = pos_body[1][0]
    y_poc = pos_body[1][0]

    "Z params"
    z_ob = pos_body[2][0]
    # z_ob = l_ob*np.sin(y_ob)                    # 2
    print('z_ob: ', z_ob)
    z_bs = self.lbs*np.sin(self.ybs)*np.cos(self.theta_b)            # 3
    print('z_bs: ', z_bs)
    if alpha_p == 0:
        z_sn = self.lsn*np.cos(self.theta_b)                         # 4
    else:
        # print('alp_p:', alpha_p)
        z_sn = self.lsn*np.sin(np.pi/2-alpha_p)*np.cos(self.theta_b) # 4
    print('z_sn: ', z_sn)

    "Calling Function for calculating water params from dynamics"
    height = z_ob - (z_bs + z_sn + self.ht_shel)
    # print('ht: ', height)
    poc = findPOC(alpha_p, self.t, height)
    poc = poc.getPOC(alpha_p, self.t, height)
    # poc = self.findPOC(alpha_p, self.t, height)
    x_water = poc[3]                            # 10
    z_water = poc[4]                            # 5

    x_poc = x_ob + x_bs + x_sn + x_water        # 11
    z_poc = z_ob - (z_bs + z_sn + z_water)      # 6

    return x_ob, x_bs, x_sn, z_ob, z_bs, z_sn, x_water, z_water, x_poc, z_poc, y_ob, y_bs, y_sn, y_poc
            #0     #1    #2    #3    #4    #5      #6       #7      #8     #9    #10   #11   #12   #13


  def ik(self, ht, deb_x, deb_y, deb_z, vel_out):
    "Getting IK of system"
    self.vel_out = vel_out           # setting vel_out
    x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z'); p = SX.sym('p'); # to include dso or force as a limiting constraint

    "Defining objective function variables"
    zpoc = deb_z
    zob = z*sin(x)*cos(p)
    zbs = self.lbs*sin(self.ybs)*cos(self.theta_b)
    zsn = self.lsn*sin(np.pi/2-y)*cos(self.theta_b)
    zwater = zob - (zbs + zsn + zpoc)

    ypoc = deb_y
    yob = z*cos(x)*sin(p)
    ybs = z*cos(x)*sin(p)
    ysn = z*cos(x)*sin(p)
    ywater = z*cos(x)*sin(p)

    xpoc = deb_x
    xob = z*cos(x)*cos(p)
    xbs = self.lbs*cos(self.ybs)*cos(self.theta_b)
    xsn = self.lsn*cos(y)*cos(self.theta_b)
    xwater = xpoc - (xob + xbs + xsn)

    # standoff distance that takes into account of the 3d vector of the traj of the water
    d = np.sqrt(np.power(xwater,2) + np.power(ywater,2) + np.power(zwater,2))

    fpoc = (self.density*d*self.Aroi)/cos(y)

    "Defining weights"
    q = 1.0         #0.6   #0.2     #increasing q, decreases xob, increases alpha_prime      #0.9
    r = 1.0         #1.2   #0.6     #increasing r, decrease xob, decreases alpha_prime       #0.162
    s = 1.0         #0.8   #0.11    #decreasing s, increases xpoc                            #0.132
    o = 1.0

    "x = y_ob, y = alpha_p, z = l_ob, p = psi_ob"
    # obj = q*zwater + r*ywater + s*xwater + o*z
    obj = y + z -fpoc

    ineq = cs.vertcat(
                      z*sin(x),         #constrains the UAV altitude
                      xob,              #constrains the UAV xpoc
                      x,                #constrains the UAV y_ob
                      y,                #constrains nozzle angle, alpha_prime
                      yob,              #constrains the UAV ypoc
                     )

    nlp = {'x': vertcat(x, y, z, p), 'f': obj, 'g': ineq}
    opts = {'ipopt.max_iter': 2000, 'ipopt.acceptable_tol': 1e-20}
    S = nlpsol('S', 'ipopt', nlp, opts)

    initial_guess = [1, 0.0174532, 5, 0.174532]

    "x = y_ob, y = alpha_p, z = l_ob, p = psi_ob"
    # lower_bound = [0.0, deb_x/4, 0.01, ht, deb_x-0.1, deb_z, 0.0, 0.0]
    # upper_bound = [0.0, deb_x/2, 45*3.1415926/180, 10.0, deb_x+0.1, deb_z, 0.1, 0.1]
    # lower_bound = [1.4, 0, 1*3.1415926/180, ht-0.1, deb_x-0.1, deb_z-0.1, 0.0, 0.0]
    # upper_bound = [1.5, deb_z, 45*3.1415926/180, 1.5*ht, deb_x, deb_z, 0.1, 0.1]

    lower_bound = [ht, deb_x-0.1, np.deg2rad(10), np.deg2rad(0), deb_y-0.1]
    upper_bound = [ht, deb_x-0.1, np.deg2rad(90), np.deg2rad(20), deb_y]

    # lower_bound = [ht, np.deg2rad(10), np.deg2rad(0), deb_y]
    # upper_bound = [ht, np.deg2rad(90), np.deg2rad(70), deb_y]
    sol = S(x0=initial_guess, lbg=lower_bound, ubg=upper_bound)
    sol_opt = sol['x']

    x = array(sol_opt[0])
    y = array(sol_opt[1])
    z = array(sol_opt[2])
    p = array(sol_opt[3])

    print('sol_opt:', sol_opt)
    print('yob: %f' % np.rad2deg(x), 'alp_prime: %f' % np.rad2deg(y), 'lob: %f' % z,
          'psi: %f' % np.rad2deg(p))
    # print('yob: %f' % (x*180/3.1415926), 'alp_prime: %f' % (y*180/3.1415926), 'lob: %f' % z,
    #       'psi: %f' % (p*180/3.1415926))
    print('ik error func: ', sol['f'])

    return x, y, z, p