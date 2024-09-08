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
        return 0.125

def u0(x):
    if x < 0.5:
        return 0
    else:
        return 0

def p0(x):
    if x < 0.5:
        return 1
    else:
        return 0.1
    
def f1(U):
    return U[1]

def f2(U):
    gammar = 1.4
    rho = U[0]
    u = U[1] / U[0]
    p = (gammar - 1) * (U[2] - 0.5 * U[1] * U[1] / U[0])
    return rho * u * u + p

def f3(U):
    gammar = 1.4
    u = U[1] / U[0]
    E = U[2]
    p = (gammar - 1) * (U[2] - 0.5 * U[1] * U[1] / U[0])
    return u * (E + p)
    
def wavespeed(u):
    gammar = 1.4
    p =  (gammar - 1) * (u[2] - 0.5 * u[1] * u[1] / u[0])
    c = np.sqrt(gammar * p / u[0])
    SR = u[1] / u[0] + c
    SL = u[1] / u[0] - c
    return SR, SL

def get_u(uh,phi_Gauss):
    NumEq = uh.shape[0]
    Nx = uh.shape[1]
    dimPk = uh.shape[2]
    NumGLP = 2 * dimPk - 1
    u_num = np.zeros((NumEq, Nx, NumGLP))

    # for i in range(NumEq):
    #     for j in range(Nx):
    #         for k in range(NumGLP):
    #             for l in range(dimPk):
    #                 u_num[i,j,k] += uh[i,j,l] * phi_Gauss[k,l]
    for i in range(NumEq):
        u_num[i] = np.dot(uh[i], phi_Gauss.T)
    
    return u_num

def LF_Flux(uR,uL,fR,fL,SR,SL):
    alpha = max(abs(SR),abs(SL))
    fhat = 0.5*(fR + fL) - 0.5*alpha*(uR - uL)
    return fhat



def Lh(uh,bcL,bcR,phi_Gauss,phi_grad_gauss,weight,phi_l,phi_r,flux_type,hx,mass):
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
    uh_gauss = get_u(uh,phi_Gauss) # 获取高斯点的数值
    
    f1_gauss = np.zeros((Nx, NumGLP))
    f2_gauss = np.zeros((Nx, NumGLP))
    f3_gauss = np.zeros((Nx, NumGLP))
    for i in range(Nx):
        for k in range(NumGLP):
            u = uh_gauss[:,i,k]
            f1_gauss[i,k] = f1(u)
            f2_gauss[i,k] = f2(u)
            f3_gauss[i,k] = f3(u)
    
    weight_matrix = np.diag(weight)
    du = np.zeros((NumEq, Nx, dimPk))
    # for i in range(NumEq):
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

    fR = np.zeros(NumEq)
    fL = np.zeros(NumEq)
    fhat = np.zeros((NumEq, Nx + 1))
    for i in range(Nx+1):
        uR = uhL[:,i]
        uL = uhR[:,i]
        fR[0] = f1(uR)
        fR[1] = f2(uR)
        fR[2] = f3(uR)
        fL[0] = f1(uL)
        fL[1] = f2(uL)
        fL[2] = f3(uL)
        SLmax,SLmin = wavespeed(uL)
        SRmax,SRmin = wavespeed(uR)
        SR = max(SLmax,SRmax)
        SL = min(SLmin,SRmin)

        if flux_type == 1:
            fhat[:,i] = LF_Flux(uR,uL,fR,fL,SR,SL)


    for i in range(NumEq):
        for j in range(Nx):
            for k in range(dimPk):
                du[i,j,k] -= (1/hx) * (phi_r[k]*fhat[i,j+1] - phi_l[k]*fhat[i,j])
    
    for d in range(dimPk):
        du[:,:,d] /= mass[d]
    
    return du


def TVD_Limiter_P1(uh,bcL,bcR,hx):
    NumEq = uh.shape[0]
    Nx = uh.shape[1]
    dimPk = uh.shape[2]

    uhmod = np.zeros_like(uh)
    uhmod[:,:,0] = uh[:,:,0]

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

    deltaU = np.zeros(NumEq)
    deltaURM = np.zeros(NumEq)
    deltaULM = np.zeros(NumEq)
    
    for i in range(Nx):
        for n in range(NumEq):
            deltaU[n] = uh[n,i,1]
            deltaURM[n] = uhb[n,i+2,0] - uhb[n,i+1,1]
            deltaULM[n] = uhb[n,i+1,0] - uhb[n,i,1]

        v = uh[1,i,0] / uh[0,i,0]
        p = get_pressure(uh[0,i,0],uh[1,i,0],uh[2,i,0])
        c = np.sqrt(1.4 * p / uh[0,i,0])
        H = (uh[2,i,0] + p) / uh[0,i,0] # H = (E + p) / rho

        R = np.array([[1,1,1],[v-c,v,v+c],[H-c*v,0,H+c*v]])
        L = np.linalg.inv(R)

        deltaU1 = R @ minmod(L @ deltaU, L @ deltaURM, L @ deltaULM, hx)

        if np.linalg.norm(deltaU1 - deltaU) > 1e-6:
            for n in range(NumEq):
                uhmod[n,i,1] = deltaU1[n]
                uhmod[n,i,2] = 0

    uh = uhmod
    return uh

def get_pressure(rho,rhou,E):
    gammar = 1.4
    p = (gammar - 1) * (E - 0.5 * rhou * rhou / rho)
    return p

def minmod(a,b,c,hx):
    NumEq = a.shape[0]
    a1 = np.zeros(NumEq)
    M = 1

    for i in range(NumEq):
        if np.abs(a[i]) < M * hx * hx:
            a1[i] = a[i]
        else:
            if np.sign(a[i]) == np.sign(b[i]) and np.sign(b[i]) == np.sign(c[i]):
                a1[i] = np.sign(a[i]) * min(np.abs(a[i]),np.abs(b[i]),np.abs(c[i]))
            else:
                a1[i] = 0

    return a1