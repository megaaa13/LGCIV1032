"""
rupture.py

Based on the lua code in the course "LGCIV1032 - Structures en béton armé"
Author : Jean-François Cap

Translation of the .lua code in python by Martin G.
"""

from math import ceil, pi, sqrt, pow, cos

import matplotlib.pyplot as plt
import numpy as np

# Classes for data
class Reinforcement:
    def __init__(self, y, A):
        self.y = y
        self.A = A

class Section:
    def __init__(self, b, fck, h, reinforcements, dx=None,
                 fyk=500, alpha_cc=0.85, 
                 gamma_c=1.5, gamma_s=1.15, 
                 euk=5/100, k=1.08,
                 Es=200e3):
        self.b = b
        self.h = h
        self.dx = dx
        self.fck = fck
        self.fyk = fyk
        self.alpha_cc = alpha_cc
        self.gamma_c = gamma_c
        self.gamma_s = gamma_s
        self.euk = euk
        self.k = k
        self.Es = Es
        self.reinforcements = reinforcements

    def dxOrH(self):
        return self.dx if self.dx is not None else self.h / 1000
    
class Data:
    def __init__(self, section: Section, N=0, NEd=None, y0=None):
        self.section = section
        self.N = N
        self.h = self.section.h
        self.y0 = y0
        self.NEd = NEd
    
    def y0orH(self):
        return self.y0 if self.y0 is not None else self.h / 2

# Functions
def interpolate(x1, y1, x2, y2):
    def function(x):
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    return function

def primitive(f, fromm, dx):
    def function(x):
        integral = 0
        steps = ceil((x - fromm) / dx)
        if steps == 0:
            return 0
        dxx = (x - fromm) / steps
        for i in range(1, steps + 1):
            integral += f(fromm + (i-0.5) * dxx) * dxx
        return integral
    return function

def map(table, f, init):
    res = init
    for i, v in enumerate(table):
        res = f(i, v, res)
    return res

def searchValueInRows(rows, key, val):
    for i in range(1, len(rows)):
        v1 = rows[i][key]
        v2 = rows[i+1][key]
        if v1 <= val and val <= v2:
            res = {}
            alpha = 1 - (val - v1) / (v2 - v1)
            for k, va in rows[i].items():
                vb = rows[i+1][k]
                if isinstance(va, float) and isinstance(vb, float):
                    res[k] = alpha*va + (1-alpha)*vb
            return res

def StrainStressConcrete(fck, alpha_cc, gamma_c):
    ec2 = 2/1000 if fck <= 50 else ((2 + 0.085 * pow(fck, 0.53)) / 1000)
    ecu2 = 3.5/1000 if fck <= 50 else ((2.6 + 35 * pow((90-fck)/100, 4)) / 1000)
    n = 2 if fck <= 50 else (1.4 + 23.4 * pow((90-fck)/100, 4))
    fcd = alpha_cc * fck / gamma_c
    def function(epsilon):
        if epsilon <= 0:
            return 0
        elif epsilon <= ec2:
            return (1-pow(1-epsilon/ec2, n))*fcd
        elif epsilon <= ecu2:
            return fcd
        else:
            return fcd
        
    return function, ec2, ecu2

def StrainStressSteel(fyk, gamma_s, euk, k, E):
    fyd, ftd, eyd, eud = fyk / gamma_s, fyk * k / gamma_s, fyk / (gamma_s * E), euk * 0.8
    dfd = (ftd - fyd) / (euk - eyd)
    sud = fyd + dfd * (eud - eyd)
    def function(epsilon):
        if epsilon < - eud:
            return - sud
        elif epsilon < - eyd:
            return - fyd + (epsilon + eyd) * dfd
        elif epsilon < eyd:
            return epsilon * E
        elif epsilon <= eud:
            return fyd + (epsilon - eyd) * dfd
        else:
            return sud
    
    return function, eyd, eud

def UltimeStrains(d, h, ec2, ecu2, eyd, eud):
    return [
        lambda y: -eud, # 1
        interpolate(0, 0, d, -eud), # 2
        interpolate(0, ec2, d, -eud), # 2'
        interpolate(0, ecu2, d, -eud), # 3
        interpolate(0, ecu2, 0.15*d, 0), # x/d = 0.15
        interpolate(0, ecu2, 0.25*d, 0), # x/d = 0.25
        interpolate(0, ecu2, 0.35*d, 0), # x/d = 0.35
        interpolate(0, ecu2, 0.45*d, 0), # x/d = 0.45
        interpolate(0, ecu2, d, -eyd), # 4
        interpolate(0, ecu2, d, -eyd/2), # 4/5
        interpolate(0, ecu2, d, 0), # 5
        interpolate(0, ecu2, h, 0), # 6
        lambda y: ec2 # 7
    ]
    
def Interaction(data: Data):
    section = data.section
    fck, alpha_cc, gamma_c = section.fck, section.alpha_cc, section.gamma_c
    fyk, gamma_s, euk = section.fyk, section.gamma_s, section.euk
    k, Es, reinf = section.k, section.Es, section.reinforcements
    b, y0, dx = section.b, data.y0orH(), section.dxOrH()
    h, d = section.h, map(reinf, lambda i, r, d: max(r.y, d), 0)
    #? strain stress relations
    sigma_c, ec2, ecu2 = StrainStressConcrete(fck, alpha_cc, gamma_c)
    sigma_s, eyd, eud = StrainStressSteel(fyk, gamma_s, euk, k, Es)
    #? build strain(y) functions
    strains = UltimeStrains(d, h, ec2, ecu2, eyd, eud)
    #? Compute forces by stresses integration
    interaction = []
    for i in range(len(strains)):
        strain = strains[i]
        ec, es = strain(0), strain(d)
        if abs(ec-es) < 1e-6:
            value = 1e10
        else:
            value = ec/(ec-es)*d
        x = min(max(0, value), h)
        Fc = primitive(lambda y: sigma_c(strain(y))*b(y)*1000, 0, dx)(x)
        Mc = primitive(lambda y: sigma_c(strain(y))*b(y)*(y0-y)*1000, 0, dx)(x)
        Fs = map(reinf, lambda i, r, F: F+r.A*sigma_s(strain(r.y))*1000, 0)
        Ms = map(reinf, lambda i, r, M: M+r.A*sigma_s(strain(r.y))*1000*(y0-r.y), 0)
        NRd, MRd = Fc + Fs, Mc + Ms
        interaction.append({
            'ec': ec, 'es': es, 'x': x,
            'Fc': Fc, 'Mc': Mc, 'Fs': Fs, 'Ms': Ms,
            'NRd': NRd, 'MRd': MRd,
            'h': h, 'd': d
        })
    return interaction

#? ULS Resistance check

def Resistance(data: Data):
    interaction = Interaction(data)
    print("ULS Interaction")
    abs = np.zeros(len(interaction)*2) # Pour plot le diagramme d'interaction
    abs1 = np.zeros(len(interaction)*2) # Pour plot le diagramme d'interaction
    ord = np.zeros(len(interaction)*2) # Pour plot le diagramme d'interaction
    abs2 = np.zeros(len(interaction)*2) # Pour plot le diagramme d'interaction
    for i, s in enumerate(interaction):
        abs[i], abs[len(interaction) + i] = s['NRd'], s['NRd'] # Pour plot le diagramme d'interaction
        abs
        ord[i], ord[len(interaction) + i] = s['MRd'], -s['MRd'] # Pour plot le diagramme d'interaction
        print(f"ec={s['ec']*100:6.3f}% es={s['es']*100:6.3f}% x/d={s['x']/s['d']:6.3f} Fc={s['Fc']:8.1f}kN Fs={s['Fs']:8.1f}kN NRd={s['NRd']:8.1f}kN MRd={s['MRd']:8.1f}kNm")
    #Plot du diagramme d'interaction
    plt.plot(abs[:len(interaction)], ord[:len(interaction)], 'blue', linestyle='dashed')
    plt.plot(abs[len(interaction):], ord[len(interaction):], 'blue', linestyle='dashed')
    plt.scatter(abs, ord, color='salmon')
    if data.NEd is not None:
        print("ULS Check")
        NEd, NRdmin, NRdmax = data.NEd, interaction[0]['NRd'], interaction[-1]['NRd']
        if (NEd < NRdmin):
            print(f"NEd={NEd:8.1f}kN < NRdmin={NRdmin:8.1f}kN")
        elif (NEd > NRdmax):
            print(f"NEd={NEd:8.1f}kN > NRdmax={NRdmax:8.1f}kN")
        else:
            s = searchValueInRows(interaction, 'NRd', NEd)
            print(f"MRd(NEd={NEd:8.1f}) = {s['MRd']:8.1f}kNm")
        print(f"ec={s['ec']*100:6.3f}% es={s['es']*100:6.3f}%, x/d={s['x']/s['d']:6.3f}, Fc={s['Fc']:8.1f}kN, Fs={s['Fs']:8.1f}kN, NRd={s['NRd']:8.1f}kN, MRd={s['MRd']:8.1f}kNm")
        
#? Calcul section rectangulaire
data = Data(
    NEd=400, # Résistance en flexion composée NEd = 400kN
    section=Section(
        h=0.5, # Hauteur 50cm
        b=lambda y: 0.3, # Largeur constante 30cm
        fck=30,
        fyk=500,
        reinforcements=[
            Reinforcement(0.45, 5 * 16**2 * pi / 4 * 1e-6) # 5 barres de 16 en y = 45cm
        ]
    )
)
# Resistance(data)

#? Calcul section circulaire
r = 0.25
rs, n, phi = 0.2, 10, 16 # 10 barres de 16mm distribuée sur 0.2m de rayon

data = Data(
    NEd=300,  
    section=Section(
        fck=30,
        fyk=500, 
        b=lambda y: 2*sqrt(r**2-(r-y)**2), 
        h=2*r, 
        reinforcements=[
            Reinforcement(r- rs * cos(i*2*pi/n), phi**2 * pi / 4 * 1e-6) for i in range(n)
        ]
    )
)

# Resistance(data)

#? Section en T
data = Data(
    NEd=1200,
    section=Section(
        fck=30,
        fyk=500,
        b=lambda y: 1.2 if y <= 0.12 else 0.2 if y <= 0.6 else 0.4,
        h=0.75,
        reinforcements=[
            Reinforcement(0.04, 1664e-6),
            Reinforcement(0.68, 5485e-6)
        ]
    )
)

# Resistance(data)

#? Béton sup, section en T
data = Data(
    NEd=1200,
    section=Section(
        fck=35,
        fyk=500,
        b=lambda y: 1.2 if y <= 0.12 else 0.2 if y <= 0.6 else 0.4,
        h=0.75,
        reinforcements=[
            Reinforcement(0.04, 8 * 12**2 * pi / 4 * 1e-6),
            Reinforcement(0.68, 5356e-6)
        ]
    )
)

# Resistance(data)

#? Test exo 4.1
data = Data(
    section=Section(
        fck=40,
        fyk=500,
        b=lambda y: 0.8,
        h=0.88,
        reinforcements=[
            Reinforcement(0.08, 12 * 40**2 * pi / 4 * 1e-6),
            Reinforcement(0.88-0.08, 12 * 40**2 * pi / 4 * 1e-6)
        ]
    )
)
Resistance(data)


#Valeurs exo 8.3, flambement
N = [6375, 5175, 4575]
M = [2497, 9227, 8993]
plt.scatter(N, M)
plt.grid()
plt.show()

#? Nouvelle section (8.3.4)
data = Data(
    section=Section(
        fck=40,
        fyk=500,
        b=lambda y: 0.8,
        h=1,
        reinforcements=[
            Reinforcement(0.08, 16 * 40**2 * pi / 4 * 1e-6),
            Reinforcement(1-0.08, 16 * 40**2 * pi / 4 * 1e-6)
        ]
    )
)
Resistance(data)
#Valeurs exo 8.3, flambement
N = [6375, 5175, 4575]
M = [2497, 9227, 8993]
plt.scatter(N, M)
plt.grid()
plt.show()