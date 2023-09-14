amideangle = 120
amidebond  = 1.334
nhbond     = 1.099
pobond     = 1.56
pnbond     = 1.56

ang2bohr = 1.8897259886
eh2kjmol = 2625.5002
nm2bohr  = ang2bohr*10.0

covrad = {
   'H' : 0.32,                                                                         'He': 0.37,
   'Li': 1.30, 'Be': 0.94, 'B' : 0.77, 'C' : 0.75, 'N' : 0.71, 'O' : 0.64, 'F' : 0.64, 'Ne': 0.62,
   'Na': 1.60, 'Mg': 1.40, 'Al': 1.24, 'Si': 1.14, 'P' : 1.09, 'S' : 1.04, 'Cl': 1.00, 'Ar': 1.01,
   'K' : 2.00, 'Ca': 1.74, 
   'Sc': 1.70, 'Ti': 1.60, 'V' : 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
                           'Ga': 1.23, 'Ge': 1.20, 'As': 1.20, 'Se': 1.18, 'Br': 1.17, 'Kr': 1.16,
   'Rb': 2.15, 'Sr': 1.90, 
   'Y' : 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
                           'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.37, 'I' : 1.36, 'Xe': 1.36,
   'Cs': 2.44, 'Ba': 2.15,
   'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
   'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W' : 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
   '                        Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50
}

mass = {
        "H": 1.00797,
        "C": 12.011,
        "N": 14.0067,
        "O": 15.9994,
        "S": 32.06
        }

ztosym = {
    1 : 'H' ,                                                           2 : 'He', 
    3 : 'Li', 4 : 'Be', 5 : 'B' , 6 : 'C' , 7 : 'N', 8 : 'O', 9 : 'F' , 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 
    19: 'K' , 20: 'Ca',                                       35: 'Br', 
                                                              53: 'I' }

symtoz = { 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O':8, 'f':9, 'Ne':10,
          'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18,
          'K': 19, 'Ca':20, 'Br':35, 'I':53 }

nwchem_input = """
echo
title {}

start
memory {} mb

geometry noautoz noautosym
 load "{}"
end


xtb
 acc 1.0
end

driver
 clear
 maxiter 100
end

task xtb optimize ignore

basis "ao basis" spherical
 * library {}
end

basis "cd basis" spherical
 * library {}
end

dft
 adft
 mult {}
 xc {}
 movecs input atomic output "dftmos.movecs"
end

driver
  maxiter 20
end


task dft optimize ignore

esp
end

task esp
"""
