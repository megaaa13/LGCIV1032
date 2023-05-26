"""
flexion_composee.py

Based on the lua code in the course "LGCIV1032 - Structures en béton armé"
Author : Jean-François Cap

Translation of the .lua code in python by Martin G.
"""

from math import ceil, pi, sqrt, cos

# Classes for data
class Section:
    def __init__(self, b, m, h, reinforcements, dx=None):
        self.b = b
        self.m = m
        self.dx = dx
        self.h = h
        self.reinforcements = reinforcements

    def dxOrH(self):
        return self.dx if self.dx is not None else self.h / 1000
    
class Data:
    def __init__(self, section: Section, N=0, M=0, y0=None, h=None, itermax=100, precision=1e-4):
        self.sections = section
        self.N = N
        self.M = M
        self.y0 = y0
        self.h = section.h
        self.itermax = itermax
        self.precision = precision
    
    def y0orH(self):
        return self.y0 if self.y0 is not None else self.h / 2

# Functions
def map(table, f, init):
    res = init
    for i, v in enumerate(table):
        res = f(i, v, res)
    return res

def primitive(f, fromm, dx):
    def function(x):
        integral = 0
        steps = ceil((x - fromm) / dx)
        dxx = (x - fromm) / steps
        for i in range(1, steps + 1):
            integral += f(fromm + (i-0.5) * dxx) * dxx
        return integral
    return function

def interpolate(x1, y1, x2, y2):
    def function(x):
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    return function

def elastic_prop(section: Section):
    b = section.b
    m = section.m
    dx = section.dxOrH()
    reinf = section.reinforcements

    Ac = primitive(lambda y: b(y), 0, dx)
    Acy = primitive(lambda y: b(y) * y, 0, dx)
    Acyy = primitive(lambda y: b(y) * y * y, 0, dx)

    As = map(reinf, lambda i, r, A: A+r['A'], 0)
    Asy = map(reinf, lambda i, r, Ay: Ay+r['A']*r['y'], 0)
    Asyy = map(reinf, lambda i, r, Ayy: Ayy+r['A']*r['y']*r['y'], 0)

    def func(y):
        At = m*As+Ac(y)
        Aty = m*Asy+Acy(y)
        Atyy = m*Asyy+Acyy(y)
        return At, Aty/At, Atyy-Aty**2/At
    return func

def elastic(data: Data):
    section = data.sections
    N = data.N
    M = data.M
    y0 = data.y0orH()
    itermax = data.itermax
    precision = data.precision
    print(f'N = {N:5g} M = {M:5g} y0 = {y0:5g}')
    h = section.h
    m = section.m
    x = h
    prop = elastic_prop(section)
    x0, inter, MGt, At, yGt, It = 0, 0, 0, 0, 0, 0
    while (abs(x - x0)/h > precision or inter > itermax):
        inter += 1
        x0 = x
        At, yGt, It = prop(x0)
        MGt = N * (yGt-y0) + M
        if MGt >= 0:
            x = max(min(yGt+(N*It)/(MGt*At), h), 0)
        else:
            raise Exception('MGt < 0')
        print(f'[{inter}] x = {x:6.4f} yGt = {yGt:6.4f} At = {At:8.5g} It = {It:8.5e} MGt = {MGt:8.5g}')
    if (inter < itermax):
        print(f"Concrete stress = {N/At+MGt/It*yGt:.1f}")
        for i, r in enumerate(section.reinforcements):
            print(f"Reinforcement [{i+1}] stress (y={r['y']:6.3f}) = {m*(N/At+MGt/It*(yGt-r['y'])):6g}")
    else:
        print("No convergence")

# Calcul sur section rectangulaire

data = Data(
    N=90, 
    M=130, 
    section=Section(
        m=6.06, 
        b=lambda y: 0.3, 
        h=0.5, 
        reinforcements=[
            { 'y': 0.45, 'A': 5*16**2 * pi / 4 * 1e-6 },
            { 'y': 0.05, 'A': 4*12**2 * pi / 4 * 1e-6 }
        ]
    )
)

# elastic(data)

# Calcul sur section circulaire
r = 0.25
rs, n, phi = 0.2, 10, 16 # 10 barres de 16mm distribuée sur 0.2m de rayon

reinforcements = []
for i in range(n):
    reinforcements.append({
        'y': r- rs * cos(i*2*pi/n),
        'A': phi**2 * pi / 4 * 1e-6
    })

data = Data(
    N=200, 
    M=90, 
    section=Section(
        m=6.06, 
        b=lambda y: 2*sqrt(r**2-(r-y)**2), 
        h=2*r, 
        reinforcements=reinforcements
    )
)

elastic(data)

