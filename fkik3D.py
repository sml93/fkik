import math
import numpy as np
import sympy as sp
import casadi as cs
from casadi import *
from numpy import array
from sympy import Symbol
from scipy import linalg as lin
from matplotlib import animation
from mpl_toolkits import mplot3d
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from sympy.printing.latex import LatexPrinter, print_latex
from sympy.core.function import UndefinedFunction, Function

t = np.linspace(0,2,1000)


class kinematics():
    def __init__(self, angle, t, ht_shel):
        self.ybs = np.deg2rad(73.333)
        self.lbs = 0.052276
        self.lsn = 0.12
        self.vel_out = 146
        self.g = 9.81
        self.drag = 1.125
        self.z_0 = 0
        self.m = 0.0005
        self.D = 0.003
        self.C = 37
        self.c = self.C*(np.power(self.D,2))
        self.psi_ob = np.arctan(deb_y/deb_x)
        self.theta_b = np.deg2rad(0)
        # self.psi_ob = np.deg2rad(20)
        print(self.psi_ob*180/np.pi)

    def modelz(self, z, t):
        vz0 = z[0]
        dvdt = self.g - self.c*vz0/self.m
        dzdt = vz0
        return dvdt, dzdt

    def modelx(self, x, t):
        vx0 = x[0]
        # print(vx0)
        dvdt = -self.c*vx0/self.m
        dxdt = vx0
        return dvdt, dxdt

    def homotrans_ob(self, gamma_ob, psi_ob, lob):
        hob = [[np.cos(psi_ob), np.sin(psi_ob), 0, lob * np.cos(gamma_ob) * np.cos(psi_ob)],
               [np.sin(psi_ob), np.cos(psi_ob), 0, lob * np.cos(gamma_ob) * np.sin(psi_ob)],
               [0, 0, 1, lob * np.sin(gamma_ob) * np.cos(psi_ob)],
               [0, 0, 0, 1]]

        p_b = [[0], [0], [0], [1]]

        p_o = np.dot(hob, p_b)
        # print(p_o)
        return p_o[0], p_o[1], p_o[2]

    def homotrans_bs(self, theta, gamma_ob, psi_ob, lob):
        hbs = [[np.cos(theta), 0, np.sin(theta), self.lbs*np.cos(self.ybs)*np.cos(theta)],
               [0, 1, 0, lob*np.cos(gamma_ob)*np.sin(psi_ob)],
               [-np.sin(theta), 0, np.cos(theta), self.lbs*np.sin(self.ybs)*np.cos(theta)],
               [0, 0, 0, 1]]

        p_s = [[0], [0], [0], [1]]
        p_b = np.dot(hbs, p_s)

        return p_b[0], p_b[1], p_b[2]

    def homotrans_sn(self, theta, alpha_p, gamma_ob, psi_ob, lob):
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

    def solve(self, modelz, modelx, ang, t):
        "Solving ODE"
        # ang = np.deg2rad(ang)
        print('ang(rad):', ang)
        vz0 = self.vel_out*np.cos(ang)
        vx0 = self.vel_out*np.sin(ang)
        z = [vz0, 0]
        x = [vx0, 0]

        sol = odeint(modelz, z, t)
        sol1 = odeint(modelx, x, t)

        "Getting all values for x"
        vx1 = sol1[:, 0]
        x = sol1[:, 1]

        "Getting all values for z and multiplying by -1"
        vz1 = sol[:, 0]*-1
        z = sol[:, 1]*-1
        return x,z,vx1,vz1

    def findPOC(self, ang, t, ht):
        "Getting POC params"
        ht = ht*-1
        model = self.solve(self.modelz, self.modelx, ang, t)

        x = model[0]
        z = model[1]

        vx = model[2]
        vz = model[3]

        for i in range(len(z)):
            if ht >= z[i]:
                index = i
                # print(z[i-1], z[i], z[i+1])
                break

        time_taken = t[index]
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

    def fk(self, y_ob, alpha_p, l_ob):
        "Getting FK of system"
        "Getting xob/yob/zob values"
        pos_ob = self.homotrans_ob(y_ob, self.psi_ob, l_ob)
        pos_bs = self.homotrans_bs(self.theta_b, y_ob, self.psi_ob, l_ob)
        pos_sn = self.homotrans_sn(self.theta_b, alpha_p, y_ob, self.psi_ob, l_ob)

        "X params"
        x_ob = pos_ob[0][0]                           # 7
        print("x_ob: ", x_ob)
        x_bs = pos_bs[0][0]                           # 8
        print("x_bs: ", x_bs)
        x_sn0 = self.lsn*np.cos(alpha_p-np.pi/2)      # 9
        print("x_sn0: ", x_sn0)
        x_sn = pos_sn[0][0]
        print("x_sn: ", x_sn)

        "Y params"
        y_ob = pos_ob[1][0]
        print('y_ob: ', y_ob)
        y_bs = pos_bs[1][0]
        print('y_bs: ', y_bs)
        y_sn = pos_sn[1][0]
        y_poc = pos_ob[1][0]

        "Z params"
        z_ob = pos_ob[2][0]                           # 2
        print('z_ob: ', z_ob)
        z_bs = pos_bs[2][0]                           # 3
        print('z_bs: ', z_bs)
        z_sn = pos_sn[2][0]
        if alpha_p == 0:
            z_sn1 = self.lsn                          # 4
        else:
            z_sn1 = self.lsn*np.sin(np.pi/2-alpha_p)  # 4
        print('z_sn: ', z_sn)
        print('z_sn1: ', z_sn1)

        "Calling Function for calculating water params from dynamics"
        height = z_ob - (z_bs + z_sn + ht_shel)
        # print('ht: ', height)
        poc = self.findPOC(alpha_p-self.theta_b, t, height)
        x_water = poc[3]                            # 6
        z_water = poc[4]                            # 7

        x_poc = x_ob + x_bs + x_sn + x_water        # 8
        z_poc = z_ob - (z_bs + z_sn + z_water)      # 9

        return x_ob, x_bs, x_sn, z_ob, z_bs, z_sn, x_water, z_water, x_poc, z_poc, y_ob, y_bs, y_sn, y_poc
               #0     #1    #2    #3    #4    #5      #6       #7      #8     #9    #10   #11   #12   #13

    def fkik(self, y_ob, alpha_p, l_ob, psi_ob):
        "Getting FK of system"
        "Getting xob/yob/zob values"
        pos_body = self.homotrans_ob(y_ob, psi_ob, l_ob)

        "X params"
        x_ob = pos_body[0][0]
        # x_ob = l_ob*np.cos(y_ob)                    # 7
        print("x_ob: ", x_ob)
        x_bs = self.lbs*np.cos(self.ybs)*np.cos(self.theta_b)            # 8
        print("x_bs: ", x_bs)
        x_sn = self.lsn*np.cos(alpha_p-np.pi/2)*np.cos(self.theta_b)     # 9
        print(x_sn)

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
        height = z_ob - (z_bs + z_sn + ht_shel)
        # print('ht: ', height)
        poc = self.findPOC(alpha_p, t, height)
        x_water = poc[3]                            # 10
        z_water = poc[4]                            # 5

        x_poc = x_ob + x_bs + x_sn + x_water        # 11
        z_poc = z_ob - (z_bs + z_sn + z_water)      # 6

        return x_ob, x_bs, x_sn, z_ob, z_bs, z_sn, x_water, z_water, x_poc, z_poc, y_ob, y_bs, y_sn, y_poc
               #0     #1    #2    #3    #4    #5      #6       #7      #8     #9    #10   #11   #12   #13

    def ik(self, ht, deb_x, deb_y, deb_z, vel_out):
        "Getting IK of system"
        self.vel_out = vel_out           # setting vel_out
        x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z'); p = SX.sym('p');

        "Defining objective function variables"
        zpoc = deb_z
        zob = z*sin(x)*cos(p)
        zbs = self.lbs*sin(self.ybs)*cos(self.theta_b)
        zsn = self.lsn*sin(pi/2-y)*cos(self.theta_b)
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

        "Defining weights"
        q = 1.0         #0.6   #0.2     #increasing q, decreases xob, increases alpha_prime      #0.9
        r = 1.0         #1.2   #0.6     #increasing r, decrease xob, decreases alpha_prime       #0.162
        s = 1.0         #0.8   #0.11    #decreasing s, increases xpoc                            #0.132
        o = 1.0

        "x = y_ob, y = alpha_p, z = l_ob, p = psi_ob"
        obj = q*zwater + o*ywater + r*xwater + s*z

        ineq = cs.vertcat(
                          (z*sin(x)),       #constrains the UAV altitude
                          xob,              #constrains the UAV xposition
                          x,                #constrains the UAV y_ob
                          y,                #constrains nozzle angle
                          yob,              #constrains the UAV yposition
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
        lower_bound = [ht, deb_x-0.1, 50*3.1415926/180, -20*3.1415926/180, deb_y]
        upper_bound = [ht, deb_x-0.1, 90*3.1415926/180, 20*3.1415926/180, deb_y]
        # print(sin(20))
        sol = S(x0=initial_guess, lbg=lower_bound, ubg=upper_bound)
        sol_opt = sol['x']

        x = array(sol_opt[0])
        y = array(sol_opt[1])
        z = array(sol_opt[2])
        p = array(sol_opt[3])

        print('sol_opt:', sol_opt)
        print('yob: %f' % (x*180/3.1415926), 'alp_prime: %f' % (y*180/3.1415926), 'lob: %f' % z,
              'psi: %f' % (p*180/3.1415926))
        print('ik error func: ', sol['f'])

        return x, y, z, p


def main(ang, t, ht_shel):
    run = kinematics(ang, t, ht_shel)
    # link = run.fk(np.deg2rad(60), ang, ht/np.sin(np.deg2rad(60)))  #guessing the joint angles: for gamma_ob, alpha_p, psi_ob, l_ob
    link = run.fk(np.deg2rad(57.0), ang, 5.0)
    # link = run.fk(0.964468, 0.110656, 5)
    print('xpoc: ', link[8])
    print('ypoc: ', link[13])
    print('zpoc: ', link[9])

    # ycoords = [0, link[10], link[10] + link[11], link[10] + link[11] + link[12],
    #            link[10] + link[11] + link[12] + link[13]]
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
    ax = plt.axes(projection='3d')
    ax.scatter3D(0, 0, 0, cmap='Greens', c='r')
    ax.scatter3D(x_val, y_val, z_val, cmap='Greens', c='r')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # ax.set_xlim3d(0.0, 4.0)
    # ax.set_ylim3d(0.0, 3.0)
    # ax.set_zlim3d(0.0, 6.0)
    ax.plot3D(x_val, y_val, z_val)

    # x = [0, xob[0]]
    # y = [0, yob[0]]
    # z = [0, zob[0]]
    ax.plot3D(xcoords, ycoords, zcoords, '-', c='gray')
    ax.plot3D(xcoords_water, ycoords_water, zcoords_water, '-', c='blue')

    # xcoords = [0, link[0], link[0] + link[1], link[0] + link[1] + link[2], link[0] + link[1] + link[2] + link[6]]
    # zcoords = [0, link[3], link[3] - link[4], link[3] - link[4] - link[5], link[3] - link[4] - link[5] - link[7]]
    #
    # plt.scatter(0, 0, label='origin')
    # plt.scatter(link[0], link[3], label='body', color='r')
    # plt.scatter(link[0] + link[1], link[3] - link[4], label='servo', color='g')
    # plt.scatter(link[0] + link[1] + link[2], link[3] - link[4] - link[5], label='nozzle', color='orange')
    # plt.scatter(link[0] + link[1] + link[2] + link[6], link[3] - link[4] - link[5] - link[7], label='poc', color='black')
    # plt.plot(xcoords, zcoords, '-')
    # plt.legend()
    # plt.grid()
    # plt.show()

    '''Initializing the IK method and parsing the solution into the FK to check'''
    inverse_sol = run.ik(ht, deb_x, deb_y, 2.75, 150)
    # print('norm:', (inverse_sol[1][0][0]))
    link = run.fkik(inverse_sol[0][0][0], inverse_sol[1][0][0], np.linalg.norm(inverse_sol[2][0][0]), inverse_sol[3][0][0])
    print('IK_Xpoc: ', link[8])
    print('IK_Ypoc: ', link[10])
    print('IK_Zpoc: ', link[9])
    #
    # xcoords = [0, link[0], link[0] + link[1], link[0] + link[1] + link[2], link[0] + link[1] + link[2] + link[6]]
    # zcoords = [0, link[3], link[3] - link[4], link[3] - link[4] - link[5], link[3] - link[4] - link[5] - link[7]]
    #
    # plt.scatter(0, 0, label='origin_ikfk')
    # plt.scatter(link[0], link[3], label='body_ikfk', color='r')
    # plt.scatter(link[0] + link[1], link[3] - link[4], label='servo_ikfk', color='g')
    # plt.scatter(link[0] + link[1] + link[2], link[3] - link[4] - link[5], label='nozzle_ikfk', color='orange')
    # plt.scatter(link[0] + link[1] + link[2] + link[6], link[3] - link[4] - link[5] - link[7], label='poc_ikfk', color='black')
    # plt.plot(xcoords, zcoords, '-')
    # plt.legend()
    # plt.grid()
    #
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
    ax = plt.axes(projection='3d')
    ax.scatter3D(0, 0, 0, cmap='Greens', c='k')
    ax.scatter3D(xIK_val, yIK_val, zIK_val, cmap='Greens', c='k')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # ax.set_xlim3d(0.0, 8.0)
    # ax.set_ylim3d(-4.0, 4.0)
    # ax.set_zlim3d(0.0, 8.0)
    ## ax.set_xlim3d(0.0, 4.0)
    ## ax.set_ylim3d(0.0, 3.0)
    ## ax.set_zlim3d(0.0, 6.0)
    ax.plot3D(xIK_val, yIK_val, zIK_val)

    # x = [0, xob[0]]
    # y = [0, yob[0]]
    # z = [0, zob[0]]
    ax.plot3D(xcoords, ycoords, zcoords, '-', c='gray')
    # ax.set_xlim3d(0.0, 8.0)
    # ax.set_ylim3d(-4.0, 4.0)
    # ax.set_zlim3d(0.0, 8.0)
    ax.plot3D(xcoords_water, ycoords_water, zcoords_water, '-', c='blue')
    plt.show()


if __name__ == "__main__":
    ang = np.deg2rad(5)
    ht_shel = 2.75
    deb_x = 3.0
    deb_y = 1.1
    ht = 5.0
    # deb_x = float(raw_input("deb_x: "))
    # deb_y = float(raw_input("deb_y: "))
    # ht = float(raw_input("height: "))
    main(ang, t, ht_shel)

