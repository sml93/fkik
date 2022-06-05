import math
import numpy as np
import casadi as cs

from casadi import *
from solve import solve
from numpy import array
from scipy import linalg as lin
from scipy.integrate import odeint
from matplotlib import pyplot as plt


class kinematics():
  def setUp(self):
    super(kinematics, self).setUp()
    self.init_hyperparameters()

  def tearDown(self):
    super(kinematics, self).tearDown()

  def init_hyperparameters(self):
    self.ybs = np.deg2rad(73.333)
    self.lbs = 0.052276
    self.lns = 0.12

    # self.vel_out = 146    # m/s
    self.g = 9.81
    # self.m = 0.0005
    self.deb_x = 0
    self.deb_y = 0

    # self.drag = 1.125    
    # self.D = 0.003
    # self.C = 37
    # self.c = self.C*(np.power(self.D,2))

    self.psi_ob = np.arctan(self.deb_y/self.deb_x)
    self.theta_b = np.deg2rad(0)
