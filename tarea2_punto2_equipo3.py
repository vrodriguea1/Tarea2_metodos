import meshio 
import numpy as np
import sympy as sy
from sympy import Matrix

## En esta sección se abre y se lee el documento 
malla= meshio.read('concentric_cylinders.vtk')
pts=malla.points
celdita=malla.cells 
celdita2=celdita[0]
elem=celdita2[1]
Potencial=malla.point_data
Potencial=Potencial['electrostatic_potential']
 
## Ctes 
epsilon=8.8541878128E-12
delthaV=1
Utot=0  ## Energía 

## Funciones base  
## se definen las variables simbolicas 
s = sy.Symbol('s') 
r = sy.Symbol('r')
N0 = (1/4)*(1-r)*(1-s)
N1 = (1/4)*(1+r)*(1-s)
N2 = (1/4)*(1+r)*(1+s)
N3 = (1/4)*(1-r)*(1+s)
Nt = np.array([N0, N1, N2, N3]) 
dNR=sy.diff(Nt,r); 
dNS=sy.diff(Nt,s);

## Definición del mallado (respecto a puntos y nodos)

def mapeo(l):
    nodes=elem[l]
    pts0=pts[nodes[0]]
    pts1=pts[nodes[1]]
    pts2=pts[nodes[2]]
    pts3=pts[nodes[3]]
    ptsx=np.array([pts0[0], pts1[0], pts2[0], pts3[0]])
    ptsy=np.array([pts0[1], pts1[1], pts2[1], pts3[1]])
    X=np.dot(Nt,ptsx)
    Y=np.dot(Nt,ptsy)
    return ptsx, ptsy, X, Y, nodes

## Calculo del Jacobiano 
    
def J(ptsx,ptsy):
    J0=np.dot(dNR,ptsx)
    J1=np.dot(dNR,ptsy)
    J2=np.dot(dNS,ptsx)
    J3=np.dot(dNS,ptsy)
    Jp=J=np.array([J0, J1, J2, J3])
    J=np.array([[J0, J1],[J2, J3]])
    detJ=J0*J3-J1*J2 ## Determinante del jacobiano 
    return Jp, J, detJ

## Potencial eléctrico

def potencial(nodes):
    P0=Potencial[nodes[0]]
    P1=Potencial[nodes[1]]
    P2=Potencial[nodes[2]]
    P3=Potencial[nodes[3]]
    P=np.array([P0, P1, P2, P3])
    return P

## Calculo del gradiente 

def GRA (Jp, dNR, dNS, pot, detJ):
    Grax=pot[0]*(dNR[0]*Jp[3] - dNS[0]*Jp[1]) + pot[1]*(dNR[1]*Jp[3] - dNS[1]*Jp[1]) + pot[2]*(dNR[2]*Jp[3] - dNS[2]*Jp[1]) + pot[3]*(dNR[3]*Jp[3] - dNS[3]*Jp[1])
    Gray= - pot[0]*(dNR[0]*Jp[2] - dNS[0]*Jp[0]) -pot[1]*(dNR[1]*Jp[2] - dNS[1]*Jp[0]) - pot[2]*(dNR[2]*Jp[2] - dNS[2]*Jp[0]) - pot[3]*(dNR[3]*Jp[2] - dNS[3]*Jp[0])
    Grax=(1/detJ)*Grax
    Gray=(1/detJ)*Gray
    Gradiente=np.array([Grax,Gray])
    return Gradiente

## Campo eléctrico, usando el gradiente  
    
def C_E (Gradiente):
    E=-Matrix(Gradiente)
    magnitudE = (E[0]**2 + E[1]**2)**(0.5)
    magnitudE2=magnitudE**2
    return magnitudE2

## Sumatoria de la cuadratura de Gauss, para hallar la energía potencial eléctrica

def U_E(magnitudE2, detJ):
    R = 0.5773502691896258
    S=R
    gs=magnitudE2 * detJ
    G0=gs.subs([(r,-R),(s,-S)])
    G1=gs.subs([(r,R),(s,-S)])
    G2=gs.subs([(r,R),(s,S)])
    G3=gs.subs([(r,-R),(s,S)])
    G=G0+G1+G2+G3
    U=0.5*epsilon*G
    return U

## Calculo energía total 
    
for elementos in range(0, len(elem)):
    mapp= mapeo(elementos)
    Jac= J(mapp[0], mapp[1])
    POT = potencial(mapp[4])
    Grad = GRA(Jac[0], dNR, dNS, POT, Jac[2])
    magnitudE= C_E(Grad)
    U= U_E( magnitudE, Jac[2])
    Utot= U +Utot
    
    
capacitancia=2*Utot/delthaV   ## Calculamos la capacitancia 
CapTeo = 2.266180070913597 * epsilon   # Capacitancia teórica 
error = abs((CapTeo - capacitancia) / CapTeo )* 100   ## Calculo del error respecto a un valoe teórico 
## Visulización de resultados  
print('La capacitancia obtenida es: ', capacitancia,'La capacitancia teórica es: ', CapTeo, ' Porcentaje de error: ', error )


