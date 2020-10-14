# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 21:10:34 2020

@author: clase
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sympy as sym


def sol_onda(desp_inicial, vel_inicial, carga, npuntos=100,
             niteraciones=200, longitud=1.0, area = 5, tension = 5, neumann=True):
     

    
    x = np.linspace(0, longitud, npuntos)   #Discretización de los nodos de la cuerda
    dx = x[1] - x[0]        #Tamañp del diferencial entre los nodos

    miu = np.linspace(7000,7200,npuntos)    #Vector de densidad variable
    v = np.sqrt(tension/(miu * area))       #Vector de velocidad dependiente de la densidad variable 
    dt = 0.0000001                          # Tamaño del diferencial de tiempo sin parámetro de seguridad de estabilidad
    t = np.linspace(0, dt*niteraciones, niteraciones)
    alpha = (dt/dx)
    
    # SOLUCIÓN
    u = np.zeros((niteraciones, npuntos))  #Creación de la matriz de solución

    # CONDICIONES INICIALES
    u[0, :] = desp_inicial
            
    u[1, 1:-1] = 0.5*alpha**2 * ((v[1: -1] + v[: -2])*(u[0, 2:]-u[0, 1:-1])     #Esquema para el nodo fantasma para determinar el primer valor depués de la condición inicial
                      - (v[1: -1] + v[: -2])*(u[0, 1:-1])-u[0, 0:-2]) +  \
                      - dt*vel_inicial(x)[:1:-1] + dt**2 * carga(x, 0)[1:-1]
        
                      
    #ESQUEMA EXPLICÍTO POR DIFERENCIAS FINITAS CENTRADAS
    
    for cont_t in range(2, niteraciones):
        q = carga(x, cont_t * dt)
        u[cont_t, 1:-1] = alpha**2 * 0.5*((v[1: -1] + v[2:])*(u[cont_t - 1, 2:] - u[cont_t - 1, 1:-1])  #
                - (v[1: -1] + v[: -2])*(u[cont_t-1 , 1:-1]-u[cont_t - 1, :-2]))\
                  + 2*u[cont_t - 1, 1:-1] - u[cont_t - 2, 1:-1]\
                  + dt**2 * q[1:-1]
                  

    # CONDICIÓN DE NEUMANN EN EL LADO DERECHO
        if neumann:
            u[cont_t, -1] = 2*alpha**2 * (u[cont_t - 1, -2] 
                                  - u[cont_t - 1, -1]) + 2*u[cont_t - 1, -1]\
                                  - u[cont_t - 2, -1] + dt**2 * q[-1]
                                                        

    return x, t, u, v
    

# Condiciones iniciales
def desp_inicial(x):                #Definición del desplazamiento inicial de la cuerda
    """Desplazamiento inicial de la cuerda"""
    return np.exp(-1000*(x - longitud/2)**2)


def vel_inicial(x):                 #Velocidad inicial como un vector de ceros
    """Velocidad inicial de la cuerda"""
    return np.zeros_like(x)

# Carga
def carga(x, t):                    #Función fuente o de carga 
    return 3.0 * np.ones_like(x)*t


    
#PARAMETROS DE ENTRADA

npuntos = 100
niteraciones = 200
diametro = 1.66e-3  # m
area = (np.pi/4) * diametro**2 # m^2
tension = 100 # N
longitud = 1

x, t, u, v = sol_onda(desp_inicial, vel_inicial, carga,
                            longitud=longitud,
                          cfl=1.0, neumann=True)

#ANIMACIÓN DE LA SOLUCIÓN 

def update(data, line):         # Actualizar animacion
    line.set_ydata(data)
    return line,

max_val = max(np.max(u), -np.min(u))
fig0, ax0 = plt.subplots(figsize=(6, 3))
line0, = ax0.plot(x, u[0, :])
ani = animation.FuncAnimation(fig0, update, u,interval=niteraciones,    # Animacion
                              repeat=True, fargs=(line0,))
plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Cuerda vibrando')
ax0.set_ylim(-1.2*max_val, 1.2*max_val)
plt.show()


#SOLUCIONES MANUFACTURADAS####### DEBIDO A LA CONDICIÓN DE VELOCIDAD VARIABLE LA FUNCIÓN DE CARGA 
#MANUFACTURADA RESULTA COMO UN VECTOR VARIABLE, POR LO CUAL NO FUIMOS CAPACES DE RESOLVER ESTE MÉTODO

# x = np.linspace(0, longitud, npuntos)
# T= 3
# time = np.linspace(0, niteraciones, T) 


# xi, tau = sym.symbols("xi tau")
# L, v = sym.symbols("L v")

# miu = np.linspace(7000,7200,npuntos)
# v1 = np.sqrt(tension/(miu * area))
 
# u_manu = xi*sym.exp(-3*xi)*(xi - L)*sym.sin(3*xi)**2*sym.cos(sym.pi*v/L*tau)
# carga_manu = sym.simplify(sym.diff(u_manu, xi, 2) - miu*sym.diff(u_manu, tau, 2)/v**2) #Fuente en términos de la soluición manufacturada, con la densidad variable (miu)


# u_manu_fun = sym.lambdify((xi, tau), u_manu.subs({L: 1.0, v: 1.0}), "numpy")
# carga = sym.lambdify((xi, tau), carga_manu.subs({v: 1.0, L: 1.0}), "numpy")
# u_inicial = sym.lambdify((xi), u_manu.subs({v: 1.0, L: 1.0, tau: 0}), "numpy")
# v_inicial = sym.lambdify((xi), 0*xi + u_manu.diff(tau).subs({v: 1.0, L: 1.0, tau:0}), "numpy")

# sym.plot(u_manu.subs({L:1, tau: 0}), (xi, 0, 1));



# x, t, u = sol_onda(u_inicial, v_inicial, carga, vel_onda=1.0)




#ERROR

# log_error = []
# for cont in [100, 500, 1000, 5000, 10000]:
#     x, t, u = sol_onda(u_inicial, v_inicial, carga, npuntos=cont)
#     X, T = np.meshgrid(x, t)
#     solucion_manu = u_manu_fun(X, T)
#     error_fun = solucion - solucion_manu
#     error_abs = np.linalg.norm(error_fun)
#     log_error.append(np.log10(error_abs))
    
# plt.figure()
# log_dx = np.log10(1/np.array([100, 500, 1000, 5000, 10000]))
# plt.plot(log_dx, log_error, marker="o")
# plt.xlabel(r"$\log\Delta x$")
# plt.ylabel(r"$\log\Vert u_\mathrm{manu} - u\Vert$");