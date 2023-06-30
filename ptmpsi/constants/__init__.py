amideangle = 120
amidebond  = 1.334
nhbond     = 1.099
pobond     = 1.56
pnbond     = 1.56

ang2bohr = 1.8897259886
eh2kjmol = 2625.5002
nm2bohr  = ang2bohr*10.0

covrad = {
        "H": 0.31,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "S": 1.05,
        }

mass = {
        "H": 1.00797,
        "C": 12.011,
        "N": 14.0067,
        "O": 15.9994,
        "S": 32.06
        }

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
