import numpy as np

def get_GLP(NumGLP): # 获取高斯积分点和权重
    if(NumGLP == 5):
        point = np.array([-0.9061798459386640, -0.5384693101056831 ,0, 0.5384693101056831 , 0.9061798459386640])
        weight = np.array([0.2369268850561891, 0.4786286704993665, 0.5688888888888889,  0.4786286704993665,  0.2369268850561891])
        return point, weight
    
def rho0(x):
    if x < 0.5:
        return 1
    else:
        return 0.1

def u0(x):
    if x < 0.5:
        return 0
    else:
        return 0

def p0(x):
    if x < 0.5:
        return 1
    else:
        return 0.125
    
def f1(u):
    return u[1]

def f2(u):
    gammar = 1.4
    rho = u[0]
    u = u[1] / u[0]
    p = (gammar - 1) * (u[2] - 0.5 * u[1] * u[1] / u[0])
    return rho * u * u + p

def f3(u):
    gammar = 1.4
    u = u[1] / u[0]
    E = u[2]
    p = (gammar - 1) * (u[2] - 0.5 * u[1] * u[1] / u[0])
    return u * (E + p)
    
def wavespeed(u):
    gammar = 1.4
    p =  (gammar - 1) * (u[2] - 0.5 * u[1] * u[1] / u[0])
    c = np.sqrt(gammar * p / u[0])
    SR = u[1] / u[0] + c
    SL = u[1] / u[0] - c
    return SL, SR

def get_u(uh,phi_Gauss):
    u_num = np.zeros((NumEq, Nx, NumGLP))
    NumEq = uh.shape[0]
    Nx = uh.shape[1]
    dimPk = uh.shape[2]
    NumGLP = 2 * dimPk - 1

    # for i in range(NumEq):
    #     for j in range(Nx):
    #         for k in range(NumGLP):
    #             for l in range(dimPk):
    #                 u_num[i,j,k] += uh[i,j,l] * phi_Gauss[k,l]
    for i in range(NumEq):
        u_num[i] = np.dot(uh[i], phi_Gauss.T)
    
    return u_num


def Lh(uh,bcL,bcR,phi_gauss,phi_grad_gauss,weight,point,phi_l,phi_r):
    NumEq = uh.shape[0]
    Nx = uh.shape[1]
    dimPk = uh.shape[2]
    NumGLP = 2 * dimPk - 1

    uhb = np.zeros((NumEq, Nx + 2, dimPk))
    uhb[:,1:-1,:] = uh
    if bcL == 1: # periodic boundary condition
        uhb[:,0,:] = uh[:,-1,:]
    elif bcL == 2: 
        uhb[:,0,:] = uh[:,0,:]
    
    if bcR == 1: # periodic boundary condition
        uhb[:,-1,:] = uh[:,0,:]
    elif bcR == 2:
        uhb[:,-1,:] = uh[:,-1,:]

    # 计算内部积分项
    uh_gauss = get_u(uh) # 获取高斯点的数值
    
    f1_gauss = np.zeros((NumEq, Nx, NumGLP))
    f2_gauss = np.zeros((NumEq, Nx, NumGLP))
    f3_gauss = np.zeros((NumEq, Nx, NumGLP))
    for i in range(Nx):
        for k in range(NumGLP):
            u = uh_gauss[:,i,k]
            f1_gauss[:,i,k] = f1(u)
            f2_gauss[:,i,k] = f2(u)
            f3_gauss[:,i,k] = f3(u)
    
    weight_matrix = np.diag(weight)
    du = np.zeros((NumEq, Nx, dimPk))
    for i in range(NumEq):
        du[0] = 0.5 * f1_gauss @ weight_matrix @ phi_grad_gauss
        du[1] = 0.5 * f2_gauss @ weight_matrix @ phi_grad_gauss
        du[2] = 0.5 * f3_gauss @ weight_matrix @ phi_grad_gauss

    # 数值通量
    uhR = np.zeros((NumEq, Nx + 1))
    uhL = np.zeros((NumEq, Nx + 1))
    for i in range(NumEq):
        for j in range(Nx + 1):
            for k in range(dimPk):
                uhR[i,j] += uhb[i,j,k] * phi_r[k]
                uhL[i,j] += uhb[i,j+1,k] * phi_l[k]

    for i in range(Nx+1):
        uR = uhL[:,i]
        uL = uhR[:,i]
        