import math
from Bio.SeqUtils import MeltingTemp as mt

# Molar gas constant (cal/C.mol)
R = 1.987

# Temperature (C)
T0 = -273.15

# Temperature (C)
t = -21.6

# strand concentration
y = 0.00000005

#
R_lny_4 = R * math.log(y/4)

#
T0_t = T0 + t

''' Nucleotides '''
bases = {'A':0, 'C':1, 'G':2, 'T':3}

"""
Enthalpy. 

Breslauer et al. 1986

"""
BRESLAUER_ENTHALPY_NN = [[9.1, 6.5, 7.8, 8.6],
                         [5.8, 11, 11.9, 7.8],
                         [5.6, 11.1, 11, 6.5],
                         [6, 5.6, 5.8, 9.1]]

"""
Entropy. 

Breslauer et al. 1986

"""
BRESLAUER_ENTROPY_NN = [[24, 17.3, 20.8, 23.9],
                        [12.9, 26.6, 27.8, 20.8],
                        [13.5, 26.7, 26.6, 17.3],
                        [16.9, 13.5, 12.9, 24]]


"""
Enthalpy. 

Sugimoto et al. 1996

"""
SUGIMOTO_ENTHALPY_NN = [[-8,   -9.4,  -6.6,  -5.6],
                        [-8.2, -10.9, -11.8, -6.6],
                        [-8.8, -10.3, -10.9, -9.4],
                        [-6.6, -8.8,  -8.2,  -8]]

"""
Entropy. 

Sugimoto et al. 1996

"""
SUGIMOTO_ENTROY_NN = [[-21.9, -25.5, -16.4, -15.2],
                      [-21,   -28.4, -29,   -16.4],
                      [-23.5, -26.4, -28.4, -25.5],
                      [-18.4, -23.5, -21,   -21.9]]


"""
Enthalpy. 

Santa Lucia. 1998

"""
SANTA_LUCIA_ENTHALPY_NN = [[-7.9, -8.4, -7.8,  -7.2],
                           [-8.5, -8,   -10.6, -7.8],
                           [-8.2, -9.8, -8,    -8.4],
                           [-7.2, -8.2, -8.5,  -7.9]]

"""
Entropy. 

Santa Lucia. 1998

"""
SANTA_LUCIA_ENTROPY_NN = [[-22.2, -22.4, -21,   -20.4],
                          [-22.7, -19.9, -27.2, -21],
                          [-22.2, -24.4, -19.9, -22.4],
                          [-21.3, -22.2, -22.7, -22.2]]


'''
 GC content
'''
def GC(s):  
    length = 0
    content = 0
    for i in s:
        length += 1
        if i in {'G','C'}:
            content += 1
    return (float(content)/length) * 100

'''
Melting temperature
'''
def Tm1(p):
    return mt.Tm_NN(p, dnac1=100, dnac2=100, Mg=1.5, dNTPs=0.2, saltcorr=6)


#    return (dH(p)/(dS(p) + R_lny_4)) + T0_t

'''
Melting temperature
'''
def Tm2(s):  
    at = 0
    gc = 0
    for i in s:
        if i in {'G','C'}:
            gc += 1
        if i in {'A','T'}:
            at += 1
    return (4 * gc + 2 * at)

'''
Enthalpy
'''
def dH(s, matrix=BRESLAUER_ENTHALPY_NN):
    n = len(s)
    result = 0
    for i in range(n - 1):
        result += matrix[bases[s[i]]][bases[s[i + 1]]]
    return result

'''
Entropy
'''
def dS(s, matrix=BRESLAUER_ENTROPY_NN):
    n = len(s)
    result = 0
    for i in range(n - 1):
        result += matrix[bases[s[i]]][bases[s[i + 1]]]
    return result
