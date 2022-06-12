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


class kinematics():
  # def setUp(self):
  #   super(kinematics, self).setUp()
  #   self.init_hyperparameters()

  # def tearDown(self):
  #   super(kinematics, self).tearDown()

  def __init__(self, angle, t, ht_shel, deb_x, deb_y):
    self.ybs = np.deg2rad(73.333)
    self.lbs = 0.052276
    self.lsn = 0.12

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
    self.theta_b = np.deg2rad(0)

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
    poc = findPOC(alpha_p-self.theta_b, self.t, self.ht_shel)
    poc = poc.getPOC(alpha_p-self.theta_b, self.t, self.ht_shel)
    x_water = poc[3]      # 6
    z_water = poc[4]      # 7

    x_poc = x_ob + x_bs + x_sn + x_water    # 8
    z_poc = z_ob - (z_bs + z_sn + z_water)  # 9

    return x_ob, x_bs, x_sn, z_ob, z_bs, z_sn, x_water, z_water, x_poc, z_poc, y_ob, y_bs, y_sn, y_poc
