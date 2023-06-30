import numpy as np
from ptmpsi.residues.template import Template
from ptmpsi.constants import amidebond
from ptmpsi.math import nerf

ABA = Template()
ABA.name = "ABA"
ABA.fullname = "Gamma-aminobutyric acid"
ABA.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CA" ],
["H", "HA2"],
["H", "HA3"],
["C", "C"  ],
["O", "O"  ]
])
ABA.backbone = [ABA.find("N"),ABA.find("CA"),ABA.find("C")]
ABA.coordinates = np.array([
[-2.016,   3.670,   8.556],
[-3.018,   3.793,   8.545],
[-1.221,   3.575,   9.792],
[-0.704,   2.615,   9.797],
[-0.479,   4.374,   9.795],
[-2.047,   3.677,  11.040],
[-2.562,   4.638,  11.040],
[-2.788,   2.878,  11.040],
[-1.221,   3.575,  12.289],
[-0.488,   4.381,  12.288],
[-0.696,   2.620,  12.288],
[-2.047,   3.677,  13.550],
[-3.228,   3.822,  13.550]
    ])
ABA.cattach = nerf(ABA.find_coord("O"), ABA.find_coord("CA"),
        ABA.find_coord("C"), amidebond, 120, 180)
ABA.nattach = nerf(ABA.find_coord("H"),ABA.find_coord("CG"),
        ABA.find_coord("N"), amidebond, 120, 180)


EAC = Template()
EAC.name = "Epsilon-aminocaproic acid"
EAC.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CE" ],
["H", "HE2"],
["H", "HE3"],
["C", "CD" ],
["H", "HD2"],
["H", "HD3"],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CA" ],
["H", "HA2"],
["H", "HA3"],
["C", "C"  ],
["O", "O"  ]
])
EAC.coordinates = np.array([
[6.500,   07.56,   04.85],
[05.50,   07.51,   04.86],
[07.23,   07.60,   06.11],
[07.85,   06.82,   06.13],
[07.77,   08.45,   06.13],
[06.28,   07.55,   07.31],
[05.67,   08.34,   07.28],
[05.75,   06.71,   07.28],
[07.08,   07.59,   08.62],
[07.70,   06.81,   08.65],
[07.61,   08.44,   08.65],
[06.13,   07.54,   09.82],
[05.51,   08.33,   09.79],
[05.60,   06.70,   09.79],
[06.92,   07.58,   11.13],
[07.54,   06.80,   11.17],
[07.45,   08.43,   11.17],
[06.01,   07.54,   12.32],
[04.78,   07.47,   12.22]
    ])
EAC.cattach = nerf(EAC.find_coord("O"),EAC.find_coord("CA"),EAC.find_coord("C"),
        amidebond, 120, 180)
EAC.nattach = nerf(EAC.find_coord("H"),EAC.find_coord("CE"),EAC.find_coord("N"),
        amidebond, 120, 180)



QCS = Template()
QCS.name = "QCS"
QCS.fullname = ""
QCS.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["C", "CD" ],
["O", "OE1"],
["N", "NE2"],
["H", "H21"],
["H", "H22"],
["C", "C"  ],
["O", "O"  ],
])
QCS.backbone = [QCS.find("N"), QCS.find("CA"), QCS.find("C")]
QCS.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.960286, 6.239388, 2.628674],
[3.189373, 5.830009, 3.509413],
[4.474326, 7.396902, 2.722621],
[4.250869, 7.995623, 3.530396],
[5.114681, 7.736951, 1.991040],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06]
    ])
QCS.cattach = nerf(QCS.find_coord("O"), QCS.find_coord("CA"),
        QCS.find_coord("C"), amidebond, 120, 180)
QCS.nattach = nerf(QCS.find_coord("H"), QCS.find_coord("CA"),
        QCS.find_coord("N"),amidebond,120,180)
QCS.chi1 = np.array([
    QCS.find("N"), QCS.find("CA"), QCS.find ("CB"), QCS.find("SG") ])
QCS.chi2 = np.array([
    QCS.find("CA"), QCS.find("CB"), QCS.find ("SG"), QCS.find("CD") ])



SMC = Template()
SMC.name = "SMC"
SMC.fullname = "S-methyl cysteine"
SMC.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["C", "CS" ],
["H", "HS1"],
["H", "HS2"],
["H", "HS3"],
["C", "C"  ],
["O", "O"  ],
])
SMC.backbone = [SMC.find("N"), SMC.find("CA"), SMC.find("C")]
SMC.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.289524, 6.708167, 1.878455],
[3.718115, 7.636795, 1.510345],
[2.289867, 6.583612, 1.470036],
[3.236208, 6.738847, 2.963701],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
SMC.cattach = nerf(SMC.find_coord("O"), SMC.find_coord("CA"),
        SMC.find_coord("C"), amidebond, 120, 180)
SMC.nattach = nerf(SMC.find_coord("H"), SMC.find_coord("CA"),
        SMC.find_coord("N"),amidebond,120,180)
SMC.chi1 = np.array([
    SMC.find("N"), SMC.find("CA"), SMC.find ("CB"), SMC.find("SG") ])
SMC.chi2 = np.array([
    SMC.find("CA"), SMC.find("CB"), SMC.find ("SG"), SMC.find("CS") ])



XCN = Template()
XCN.name = "XCN"
XCN.fullname = "S-cyano cysteine"
XCN.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["C", "CS" ],
["N", "NC" ],
["C", "C"  ],
["O", "O"  ],
])
XCN.backbone = [XCN.find("N"), XCN.find("CA"), XCN.find("C")]
XCN.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.960286, 6.239388, 2.628674],
[3.724793, 6.916212, 3.515271],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
XCN.cattach = nerf(XCN.find_coord("O"), XCN.find_coord("CA"),
        XCN.find_coord("C"), amidebond, 120, 180)
XCN.nattach = nerf(XCN.find_coord("H"), XCN.find_coord("CA"),
        XCN.find_coord("N"), amidebond,120,180)
XCN.chi1 = np.array([
    XCN.find("N"), XCN.find("CA"), XCN.find ("CB"), XCN.find("SG") ])
XCN.chi2 = np.array([
    XCN.find("CA"), XCN.find("CB"), XCN.find ("SG"), XCN.find("CS") ])



SNC = Template()
SNC.name = "SNC"
SNC.fullname = "S-nitroso cysteine"
SNC.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["N", "ND" ],
["O", "OE" ],
["C", "C"  ],
["O", "O"  ],
])
SNC.backbone = [SNC.find("N"), SNC.find("CA"), SNC.find("C")]
SNC.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[5.584222, 5.598917, 2.513037],
[6.023476, 4.605091, 3.297235],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
SNC.cattach = nerf(SNC.find_coord("O"), SNC.find_coord("CA"),
        SNC.find_coord("C"), amidebond, 120, 180)
SNC.nattach = nerf(SNC.find_coord("H"), SNC.find_coord("CA"),
        SNC.find_coord("N"), amidebond,120,180)
SNC.chi1 = np.array([
    SNC.find("N"), SNC.find("CA"), SNC.find ("CB"), SNC.find("SG") ])
SNC.chi2 = np.array([
    SNC.find("CA"), SNC.find("CB"), SNC.find ("SG"), SNC.find("ND") ])


CSD = Template()
CSD.name = "CSD"
CSD.fullname = "S-cysteine sulfinate"
CSD.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["O", "OD1"],
["O", "OD2"],
["C", "C"  ],
["O", "O"  ],
])
CSD.backbone = [CSD.find("N"), CSD.find("CA"), CSD.find("C")]
CSD.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.981981, 6.181260, 2.550249],
[5.248347, 5.802042, 0.293524],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
CSD.cattach = nerf(CSD.find_coord("O"), CSD.find_coord("CA"),
        CSD.find_coord("C"), amidebond, 120, 180)
CSD.nattach = nerf(CSD.find_coord("H"), CSD.find_coord("CA"),
        CSD.find_coord("N"), amidebond,120,180)
CSD.chi1 = np.array([
    CSD.find("N"), CSD.find("CA"), CSD.find ("CB"), CSD.find("SG") ])
CSD.chi2 = np.array([
    CSD.find("CA"), CSD.find("CB"), CSD.find ("SG"), CSD.find("OD1") ])



OCS = Template()
OCS.name = "OCS"
OCS.fullname = "S-cysteine sulfonate (cysteic acid)"
OCS.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["O", "OD1"],
["O", "OD2"],
["O", "OD3"],
["C", "C"  ],
["O", "O"  ],
])
OCS.backbone = [OCS.find("N"), OCS.find("CA"), OCS.find("C")]
OCS.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.823978, 5.968695, 2.631676],
[3.901581, 6.133817, 0.172615],
[5.785907, 5.065102, 1.157005],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
OCS.cattach = nerf(OCS.find_coord("O"), OCS.find_coord("CA"),
        OCS.find_coord("C"), amidebond, 120, 180)
OCS.nattach = nerf(OCS.find_coord("H"), OCS.find_coord("CA"),
        OCS.find_coord("N"), amidebond,120,180)
OCS.chi1 = np.array([
    OCS.find("N"), OCS.find("CA"), OCS.find ("CB"), OCS.find("SG") ])
OCS.chi2 = np.array([
    OCS.find("CA"), OCS.find("CB"), OCS.find ("SG"), OCS.find("OD1") ])




CSO = Template()
CSO.name = "CSO"
CSO.fullname = "S-hydroxy cysteine (Cysteine sulfenic acid)"
CSO.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["O", "OD" ],
["H", "HD" ],
["C", "C"  ],
["O", "O"  ],
])
CSO.backbone = [CSO.find("N"), CSO.find("CA"), CSO.find("C")]
CSO.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.775740, 6.034772, 2.757401],
[4.160258, 6.900628, 2.827677],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
CSO.cattach = nerf(CSO.find_coord("O"), CSO.find_coord("CA"),
        CSO.find_coord("C"), amidebond, 120, 180)
CSO.nattach = nerf(CSO.find_coord("H"), CSO.find_coord("CA"),
        CSO.find_coord("N"), amidebond,120,180)
CSO.chi1 = np.array([
    CSO.find("N"), CSO.find("CA"), CSO.find ("CB"), CSO.find("SG") ])
CSO.chi2 = np.array([
    CSO.find("CA"), CSO.find("CB"), CSO.find ("SG"), CSO.find("OD") ])




CSS = Template()
CSS.name = "CSS"
CSS.fullname = "S-mercaptocysteine"
CSS.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["S", "SG" ],
["S", "SD" ],
["H", "HD" ],
["C", "C"  ],
["O", "O"  ],
])
CSS.backbone = [CSS.find("N"), CSS.find("CA"), CSS.find("C")]
CSS.coordinates = np.array([
[3.325770, 1.547909, -1.607257E-06],
[3.909407, 0.723611, -2.739904E-06],
[3.970048, 2.845795, -1.312144E-07],
[3.671663, 3.400129, -0.889820],
[3.576965, 3.653838, 1.232143],
[2.496995, 3.801075, 1.241379],
[3.877484, 3.115795, 2.131197],
[4.309573, 5.303523, 1.366036],
[3.647106, 6.210976, 3.092670],
[4.193526, 7.441404, 3.192535],
[5.485541, 2.705207, -4.398851E-06],
[6.008824, 1.593175, -8.449829E-06],
        ])
CSS.cattach = nerf(CSS.find_coord("O"), CSS.find_coord("CA"),
        CSS.find_coord("C"), amidebond, 120, 180)
CSS.nattach = nerf(CSS.find_coord("H"), CSS.find_coord("CA"),
        CSS.find_coord("N"), amidebond,120,180)
CSS.pka = 2.04
CSS.chi1 = np.array([
    CSS.find("N"), CSS.find("CA"), CSS.find ("CB"), CSS.find("SG") ])
CSS.chi2 = np.array([
    CSS.find("CA"), CSS.find("CB"), CSS.find ("SG"), CSS.find("SD") ])



CGL = Template()
CGL.name = "CGL"
CGL.fullname = "S-glutathionyl cysteine"
CGL.elements = np.array([
["N", "N"    ],
["H", "H"    ],
["C", "CA"   ],
["H", "HA"   ],
["C", "CB"   ],
["H", "HB2"  ],
["H", "HB3"  ],
["S", "SG"   ],
["N", "N2"   ],
["H", "H21"  ],
["H", "H22"  ],
["H", "H23"  ],
["C", "CA2"  ],
["H", "HA2"  ],
["C", "CB2"  ],
["H", "HB21" ],
["H", "HB22" ],
["C", "CG2"  ],
["H", "HG21" ],
["H", "HG22" ],
["C", "C2"   ],
["O", "OC21" ],
["O", "OC22" ],
["C", "CD2"  ],
["O", "OE2"  ],
["N", "N3"   ],
["H", "H3"   ],
["C", "CA3"  ],
["H", "HA3"  ],
["C", "CB3"  ],
["H", "HB31" ],
["H", "HB32" ],
["S", "SG3"  ],
["C", "C3"   ],
["O", "O3"   ],
["N", "N4"   ],
["H", "H4"   ],
["C", "CA4"  ],
["H", "HA41" ],
["H", "HA42" ],
["C", "C4"   ],
["O", "OC41" ],
["O", "OC42" ],
["C", "C"    ],
["O", "O"    ],
])
CGL.backbone = [CGL.find("N"), CGL.find("CA"), CGL.find("C")]
CGL.coordinates = np.array([
[3.419,   1.565,  -0.055], 
[4.016,   0.752,  -0.065], 
[4.040,   2.874,  -0.044], 
[3.730,   3.432,  -0.928], 
[3.637,   3.663,   1.197], 
[2.555,   3.791,   1.210], 
[3.949,   3.122,   2.089], 
[4.342,   5.324,   1.345], 
[0.252,   4.644,   8.310], 
[0.041,   5.512,   7.823], 
[0.161,   4.594,   9.244], 
[0.112,   3.797,   7.817], 
[1.724,   4.294,   8.419], 
[1.962,   4.626,   9.411], 
[2.505,   5.022,   7.332], 
[3.317,   4.383,   6.967], 
[1.855,   5.185,   6.461], 
[3.125,   6.366,   7.734], 
[2.406,   7.169,   7.699], 
[3.437,   6.309,   8.786], 
[1.896,   2.744,   8.384], 
[2.958,   2.250,   8.819], 
[0.806,   2.214,   7.985], 
[4.297,   6.827,   6.924], 
[5.444,   6.498,   7.236], 
[4.126,   7.504,   5.760], 
[3.134,   7.651,   5.609], 
[4.938,   7.651,   4.508], 
[6.030,   7.483,   4.453], 
[4.194,   7.441,   3.193], 
[3.361,   8.153,   3.126], 
[4.875,   7.653,   2.361], 
[3.668,   6.203,   3.082], 
[4.978,   9.215,   4.533], 
[3.980,   9.814,   4.106], 
[5.934,   9.938,   5.183], 
[6.622,   9.334,   5.645], 
[5.842,  11.260,   5.839], 
[5.574,  11.062,   6.881], 
[5.024,  11.817,   5.414], 
[7.083,  12.207,   5.805], 
[7.760,  12.238,   6.860], 
[7.198,  12.905,   4.747], 
[5.559,   2.759,  -0.050], 
[6.101,   1.657,  -0.062], 
])
CGL.cattach = nerf(CGL.find_coord("O"), CGL.find_coord("CA"),
        CGL.find_coord("C"), amidebond, 120, 180)
CGL.nattach = nerf(CGL.find_coord("H"), CGL.find_coord("CA"),
        CGL.find_coord("N"), amidebond,120,180)
CGL.chi1 = np.array([
    CGL.find("N"), CGL.find("CA"), CGL.find ("CB"), CGL.find("SG") ])
CGL.chi2 = np.array([
    CGL.find("CA"), CGL.find("CB"), CGL.find ("SG"), CGL.find("SG3") ])



IYY = Template()
IYY.name = "IYY"
IYY.fullname = "S-cysteinyl cysteine"
IYY.elements = np.array([
["N", "N"    ],
["H", "H"    ],
["C", "CA"   ],
["H", "HA"   ],
["C", "CB"   ],
["H", "HB2"  ],
["H", "HB3"  ],
["S", "SG"   ],
["S", "SG2"  ],
["C", "CB2"  ],
["H", "HB21" ],
["H", "HB22" ],
["C", "CA2"  ],
["H", "HA2"  ],
["C", "C2"   ],
["O", "OC21" ],
["O", "OC22" ],
["N", "N2"   ],
["H", "H21"  ],
["H", "H22"  ],
["H", "H23"  ],
["C", "C"    ],
["O", "O"    ],
])
IYY.backbone = [IYY.find("N"), IYY.find("CA"), IYY.find("C")]
IYY.coordinates = np.array([
   [ 3.325770 ,      1.547909 ,     -0.000002],
   [ 3.909407 ,      0.723611 ,     -0.000003],
   [ 3.970048 ,      2.845795 ,     -0.000000],
   [ 3.671663 ,      3.400129 ,     -0.889820],
   [ 3.576965 ,      3.653838 ,      1.232143],
   [ 2.496995 ,      3.801075 ,      1.241379],
   [ 3.877484 ,      3.115795 ,      2.131197],
   [ 4.309573 ,      5.303523 ,      1.366036],
   [ 3.647106 ,      6.210976 ,      3.092670],
   [ 4.193526 ,      7.441404 ,      3.192535],
   [ 3.374135 ,      8.167196 ,      3.117361],
   [ 4.881162 ,      7.633216 ,      2.361455],
   [ 4.937839 ,      7.650927 ,      4.507672],
   [ 4.965952 ,      8.757173 ,      4.526293],
   [ 6.482664 ,      7.413117 ,      4.429727],
   [ 7.214125 ,      8.412117 ,      4.493729],
   [ 7.011519 ,      6.295140 ,      4.504299],
   [ 4.126262 ,      7.504348 ,      5.759700],
   [ 3.330130 ,      8.130878 ,      5.791605],
   [ 3.775557 ,      6.559438 ,      5.824886],
   [ 4.706110 ,      7.704350 ,      6.562118],
   [ 5.485541 ,      2.705207 ,     -0.000004],
   [ 6.008824 ,      1.593175 ,     -0.000008],
])
IYY.cattach = nerf(IYY.find_coord("O"), IYY.find_coord("CA"),
        IYY.find_coord("C"), amidebond, 120, 180)
IYY.nattach = nerf(IYY.find_coord("H"), IYY.find_coord("CA"),
        IYY.find_coord("N"), amidebond,120,180)
IYY.chi1 = np.array([
    IYY.find("N"), IYY.find("CA"), IYY.find ("CB"), IYY.find("SG") ])
IYY.chi2 = np.array([
    IYY.find("CA"), IYY.find("CB"), IYY.find ("SG"), IYY.find("SG2") ])

