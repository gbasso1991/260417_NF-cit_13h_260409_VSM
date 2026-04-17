#%% VSM NE@citrato 260203 & NF@citrato 260203 Febrero 2026
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import os
from sklearn.metrics import r2_score 
from mlognormfit import fit3
from mvshtools import mvshtools as mt
import re
from uncertainties import ufloat
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
#%% Funciones
def lineal(x,m,n):
    return m*x+n

def coercive_field(H, M):
    """
    Devuelve los valores de campo coercitivo (Hc) donde la magnetización M cruza por cero.
    
    Parámetros:
    - H: np.array, campo magnético (en A/m o kA/m)
    - M: np.array, magnetización (en emu/g)
    
    Retorna:
    - hc_values: list de valores Hc (puede haber más de uno si hay múltiples cruces por cero)
    """
    H = np.asarray(H)
    M = np.asarray(M)
    hc_values = []

    for i in range(len(M)-1):
        if M[i]*M[i+1] < 0:  # Cambio de signo indica cruce por cero
            # Interpolación lineal entre (H[i], M[i]) y (H[i+1], M[i+1])
            h1, h2 = H[i], H[i+1]
            m1, m2 = M[i], M[i+1]
            hc = h1 - m1 * (h2 - h1) / (m2 - m1)
            hc_values.append(hc)

    return hc_values
#%% NF@cit-260409
data_NF = np.loadtxt('NF@cit_13h_260409.txt', skiprows=12)
H_NF = data_NF[:, 0]  # Gauss
m_NF = data_NF[:, 1]*1.00435414  #con correccion de Flavio

conc_NF = 32.8 #mg/mL Magnetita

conc_NF_mm = conc_NF/1000 #mg np/mg de solvente
masa_NF = 1#(0.1172-0.0666)*conc_NF_mm
m_NF_norm = m_NF/masa_NF #emu/g

fig, a = plt.subplots( figsize=(8, 6), constrained_layout=True)
a.plot(H_NF, m_NF, '.-', label='NF@cit_13h')
a.set_ylabel('m (emu)')
a.legend()
a.grid()
a.set_title('NF@cit_13h')
a.set_xlabel('H (G)')
a.set_ylabel('m (emu)')
plt.show()
#%% NF Normalizada por masa
fig2, b = plt.subplots( figsize=(8, 6), constrained_layout=True)

b.plot(H_NF, m_NF_norm,'.-', label=f'NE core (norm con m = {masa_NF:.1e} g)')

b.set_ylabel('m (emu/g)')
b.legend()
b.grid()
b.set_title('NF@cit_13h - Normalizada por masa')
b.set_xlabel('H (G)')
b.set_ylabel('m (emu/g)')

axins = inset_axes(b,
    width="40%",
    height="40%",
    loc='lower right',
    bbox_to_anchor=(-0.01, 0.08, 0.98, 1), 
    bbox_transform=b.transAxes,
    borderpad=0) 


# Volver a graficar las curvas en el inset
axins.plot(H_NF, m_NF_norm,'.-')

# Definir región de zoom (AJUSTAR ESTOS VALORES)
axins.set_xlim(-20,20)   # rango eje X
axins.set_ylim(-4, 4)   # rango eje Y

axins.grid()

plt.savefig('NF@cit_13h.png',dpi=300)
plt.show()
