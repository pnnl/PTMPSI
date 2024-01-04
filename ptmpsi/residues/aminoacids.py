import numpy as np
from .template import Template
from ptmpsi.constants import amidebond

aminolist = [
            "ACE",
            "ALA",
            "ARG",
            "ASH",
            "ASN",
            "ASP",
            "CYM",
            "CYS",
            "CYX",
            "GLH",
            "GLN",
            "GLU",
            "GLY",
            "HID",
            "HIE",
            "HIP",
            "HIS",
            "ILE",
            "LEU",
            "LYN",
            "LYS",
            "MET",
            "NHE",
            "NME",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL"
        ]

ACE = Template()
ACE.name     = 'ACE'
ACE.elements = np.array([
#Elements  AMBER  GROMACS   NWCHEM
["H",      "HH31"],
["C",      "CH3" ],
["H",      "HH32"],
["H",      "HH33"],
["C",      "C"   ],
["O",      "O"   ]
])
ACE.coordinates = np.array([
 [2.000001, 1.000000, -1.346410E-06],
 [2.000001, 2.090000,  1.211769E-07],
 [1.486264, 2.453849,  0.889824    ],
 [1.486259, 2.453852, -0.889820    ],
 [3.427420, 2.640795, -2.981008E-06],
 [4.390580, 1.877406, -6.602501E-06]
])
ACE.charges = np.array([
 0.1123, -0.3622, 0.1123, 0.1123, 0.5972, -0.5679   
])
ACE.backbone = np.array([0,1,4])
co = ACE.coordinates[5] - ACE.coordinates[4]
cc = ACE.coordinates[1] - ACE.coordinates[4]
bisector = np.linalg.norm(co)*cc + np.linalg.norm(cc)*co
ACE.cattach = -amidebond*bisector/np.linalg.norm(bisector) + ACE.coordinates[4]


ALA = Template()
ALA.name = 'ALA'
ALA.elements = np.array([
["N", "N"  ], 
["H", "H"  ], 
["C", "CA" ], 
["H", "HA" ], 
["C", "CB" ], 
["H", "HB1"], 
["H", "HB2"], 
["H", "HB3"], 
["C", "C"  ], 
["O", "O"  ]
])
ALA.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.075059, 4.623017,  1.205786    ],
 [2.496995, 3.801075,  1.241379    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
ALA.charges = np.array([
-0.415700, 0.271900, 0.03370, 0.08230, -0.18250, 0.0603, 0.0603, 0.0603, 0.5973, -0.5679
    ])
ALA.backbone = np.array([0,2,8])
ALA.cattach = np.array([6.204455,3.702003,-2.420137E-06])
nh = ALA.coordinates[1] - ALA.coordinates[0]
can = ALA.coordinates[2] - ALA.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
ALA.nattach = -amidebond*bisector/np.linalg.norm(bisector) + ALA.coordinates[0]



ARG = Template()
ARG.name = 'ARG'
ARG.elements = np.array([
[ "N",   "N"   ],
[ "H",   "H"   ],
[ "C",   "CA"  ],
[ "H",   "HA"  ],
[ "C",   "CB"  ],
[ "H",   "HB2" ],
[ "H",   "HB3" ],
[ "C",   "CG"  ],
[ "H",   "HG2" ],
[ "H",   "HG3" ],
[ "C",   "CD"  ],
[ "H",   "HD2" ],
[ "H",   "HD3" ],
[ "N",   "NE"  ],
[ "H",   "HE"  ],
[ "C",   "CZ"  ],
[ "N",   "NH1" ],
[ "H",   "HH11"],
[ "H",   "HH12"],
[ "N",   "NH2" ],
[ "H",   "HH21"],
[ "H",   "HH22"],
[ "C",   "C"   ],
[ "O",   "O"   ]
])
ARG.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.274186, 5.009602,  1.194577    ],
 [5.354271, 4.863178,  1.185788    ],
 [3.973781, 5.548460,  0.295972    ],
 [3.881105, 5.817645,  2.426721    ],
 [2.801135, 5.964881,  2.435959    ],
 [4.181626, 5.279602,  3.325774    ],
 [4.540320, 7.142723,  2.424483    ],
 [5.151805, 7.375492,  1.655065    ],
 [4.364284, 8.040989,  3.389382    ],
 [3.575026, 7.807606,  4.434133    ],
 [3.088949, 6.925423,  4.508848    ],
 [3.465367, 8.513631,  5.147998    ],
 [5.006254, 9.201287,  3.286991    ],
 [5.604855, 9.375325,  2.492329    ],
 [4.892216, 9.903045,  4.004368    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
ARG.backbone = np.array([0,2,-2])
ARG.charges = np.array([
-0.347900,
 0.274700,
-0.263700,
 0.156000,
-0.000700,
 0.032700,
 0.032700,
 0.039000,
 0.028500,
 0.028500,
 0.048600,
 0.068700,
 0.068700,
-0.529500,
 0.345600,
 0.807600,
-0.862700,
 0.447800,
 0.447800,
-0.862700,
 0.447800,
 0.447800,
 0.734100,
-0.589400
])
nh = ARG.coordinates[1] - ARG.coordinates[0]
can = ARG.coordinates[2] - ARG.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
ARG.nattach = -amidebond*bisector/np.linalg.norm(bisector) + ARG.coordinates[0]
ARG.cattach = np.array([6.204455,3.702003,-2.420137E-06])
ARG.chi1 = np.array([
    ARG.find("N"), ARG.find("CA"), ARG.find ("CB"), ARG.find("CG") ])
ARG.chi2 = np.array([
    ARG.find("CA"), ARG.find("CB"), ARG.find ("CG"), ARG.find("CD") ])




ASH = Template()
ASH.name = "ASH"
ASH.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["O", "OD1"],
["O", "OD2"],
["H", "HD2"],
["C", "C"  ],
["O", "O"  ]
])
ASH.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.275101, 5.011380,  1.194527    ],
 [3.669108, 5.954940,  0.620011    ],
 [5.407731, 5.091879,  1.740667    ],
 [5.742902, 5.987179,  1.652920    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
ASH.backbone = np.array([0,2,11])
nh = ASH.coordinates[1] - ASH.coordinates[0]
can = ASH.coordinates[2] - ASH.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
ASH.nattach = -amidebond*bisector/np.linalg.norm(bisector) + ASH.coordinates[0]
ASH.cattach = np.array([6.204455,3.702003,-2.420137E-06])
ASH.chi1 = np.array([
    ASH.find("N"), ASH.find("CA"), ASH.find ("CB"), ASH.find("CG") ])
ASH.chi2 = np.array([
    ASH.find("CA"), ASH.find("CB"), ASH.find ("CG"), ASH.find("OD1") ])
ASH.charges = np.array([
-0.415700,
 0.271900,
 0.034100,
 0.086400,
-0.031600,
 0.048800,
 0.048800,
 0.646200,
-0.555400,
-0.637600,
 0.474700,
 0.597300,
-0.567900
])

ASN = Template()
ASN.name = "ASN"
ASN.elements = np.array([
    ["N", "N"   ],
    ["H", "H"   ],
    ["C", "CA"  ],
    ["H", "HA"  ],
    ["C", "CB"  ],
    ["H", "HB2" ],
    ["H", "HB3" ],
    ["C", "CG"  ],
    ["O", "OD1" ],
    ["N", "ND2" ],
    ["H", "HD21"],
    ["H", "HD22"],
    ["C", "C"   ],
    ["O", "O"   ]    
    ])
ASN.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.253700, 5.017112,  1.232144    ],
 [5.005299, 5.340406,  0.315072    ],
 [3.984885, 5.817909,  2.265917    ],
 [4.408015, 6.733702,  2.314743    ],
 [3.359611, 5.504297,  2.994464    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
ASN.backbone = np.array([0,2,12])
nh = ASN.coordinates[1] - ASN.coordinates[0]
can = ASN.coordinates[2] - ASN.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
ASN.nattach = -amidebond*bisector/np.linalg.norm(bisector) + ASN.coordinates[0]
ASN.cattach = np.array([6.204455,3.702003,-2.420137E-06])
ASN.chi1 = np.array([
    ASN.find("N"), ASN.find("CA"), ASN.find ("CB"), ASN.find("CG") ])
ASN.chi2 = np.array([
    ASN.find("CA"), ASN.find("CB"), ASN.find ("CG"), ASN.find("OD1") ])
ASN.charges = np.array([
-0.415700,
 0.271900,
 0.014300,
 0.104800,
-0.204100,
 0.079700,
 0.079700,
 0.713000,
-0.593100,
-0.919100,
 0.419600,
 0.419600,
 0.597300,
-0.567900
])


ASP = Template()
ASP.name = "ASP"
ASP.elements = np.array([
    ["N","N"  ],
    ["H","H"  ],
    ["C","CA" ],
    ["H","HA" ],
    ["C","CB" ],
    ["H","HB2"],
    ["H","HB3"],
    ["C","CG" ],
    ["O","OD1"],
    ["O","OD2"],
    ["C","C"  ],
    ["O","O"  ]
    ])
ASP.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.275101, 5.011380,  1.194527    ],
 [3.669108, 5.954940,  0.620011    ],
 [5.407731, 5.091879,  1.740667    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
ASP.backbone = np.array([0,2,-2])
nh = ASP.coordinates[1] - ASP.coordinates[0]
can = ASP.coordinates[2] - ASP.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
ASP.nattach = -amidebond*bisector/np.linalg.norm(bisector) + ASP.coordinates[0]
ASP.cattach = np.array([6.204455,3.702003,-2.420137E-06])
ASP.chi1 = np.array([
    ASP.find("N"), ASP.find("CA"), ASP.find ("CB"), ASP.find("CG") ])
ASP.chi2 = np.array([
    ASP.find("CA"), ASP.find("CB"), ASP.find ("CG"), ASP.find("OD1") ])
ASP.charges = np.array([
-0.516300,
 0.293600,
 0.038100,
 0.088000,
-0.030300,
-0.012200,
-0.012200,
 0.799400,
-0.801400,
-0.801400,
 0.536600,
-0.581900
])

CYM = Template()
CYM.name = "CYM"
CYM.elements = np.array([
    ["N", "N"  ],
    ["H", "H"  ],
    ["C", "CA" ],
    ["H", "HA" ],
    ["C", "CB" ],
    ["H", "HB3"],
    ["H", "HB2"],
    ["S", "SG" ],
    ["C", "C"  ],
    ["O", "O"  ]  
    ])
CYM.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [3.877484, 3.115795,  2.131197    ],
 [2.496995, 3.801075,  1.241379    ],
 [4.309573, 5.303523,  1.366036    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
CYM.backbone = np.array([0,2,-2])
nh = CYM.coordinates[1] - CYM.coordinates[0]
can = CYM.coordinates[2] - CYM.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
CYM.nattach = -amidebond*bisector/np.linalg.norm(bisector) + CYM.coordinates[0]
CYM.cattach = np.array([6.204455,3.702003,-2.420137E-06])
CYM.chi1 = np.array([
    CYM.find("N"), CYM.find("CA"), CYM.find ("CB"), CYM.find("SG") ])
CYM.chi2 = None
CYM.charges = np.array([
-0.415700,
 0.271900,
-0.035100,
 0.050800,
-0.241300,
 0.112200,
 0.112200,
-0.884400,
 0.597300,
-0.567900
])

CYS = Template()
CYS.name = "CYS"
CYS.elements = np.array([
    ["N", "N"  ],
    ["H", "H"  ],
    ["C", "CA" ],
    ["H", "HA" ],
    ["C", "CB" ],
    ["H", "HB2"],
    ["H", "HB3"],
    ["S", "SG" ],
    ["H", "HG" ],
    ["C", "C"  ],
    ["O", "O"  ]   
    ])
CYS.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.309573, 5.303523,  1.366036    ],
 [3.725392, 5.622018,  2.517640    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
CYS.backbone = np.array([0,2,-2])
nh = CYS.coordinates[1] - CYS.coordinates[0]
can = CYS.coordinates[2] - CYS.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
CYS.nattach = -amidebond*bisector/np.linalg.norm(bisector) + CYS.coordinates[0]
CYS.cattach = np.array([6.204455,3.702003,-2.420137E-06])
CYS.chi1 = np.array([
    CYS.find("N"), CYS.find("CA"), CYS.find ("CB"), CYS.find("SG") ])
CYS.chi2 = np.array([
    CYS.find("CA"), CYS.find("CB"), CYS.find ("SG"), CYS.find("HG") ])
CYS.charges = np.array([
-0.415700,
 0.271900,
 0.021300,
 0.112400,
-0.123100,
 0.111200,
 0.111200,
-0.311900,
 0.193300,
 0.597300,
-0.567900
])

CYX = Template()
CYX.name = "CYX"
CYX.elements = np.array([
    ["N", "N"  ],
    ["H", "H"  ],
    ["C", "CA" ],
    ["H", "HA" ],
    ["C", "CB" ],
    ["H", "HB2"],
    ["H", "HB3"],
    ["S", "SG" ],
    ["C", "C"  ],
    ["O", "O"  ]   
    ])
CYX.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.309573, 5.303523,  1.366036    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
CYX.backbone = np.array([0,2,-2])
nh = CYX.coordinates[1] - CYX.coordinates[0]
can = CYX.coordinates[2] - CYX.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
CYX.nattach = -amidebond*bisector/np.linalg.norm(bisector) + CYX.coordinates[0]
CYX.cattach = np.array([6.204455,3.702003,-2.420137E-06])
CYX.chi1 = np.array([
    CYX.find("N"), CYX.find("CA"), CYX.find ("CB"), CYX.find("SG") ])
CYX.chi2 = None
CYX.charges = np.array([
-0.415700,
 0.271900,
 0.042900,
 0.076600,
-0.079000,
 0.091000,
 0.091000,
-0.108100,
 0.597300,
-0.567900])

GLH = Template()
GLH.name = "GLH"
GLH.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["C", "CD" ],
["O", "OE1"],
["O", "OE2"],
["H", "HE2"],
["C", "C"  ],
["O", "O"  ]
])
GLH.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.267328, 4.996267,  1.194946    ],
 [5.347413, 4.849843,  1.186158    ],
 [3.966923, 5.535124,  0.296342    ],
 [3.873732, 5.805369,  2.428706    ],
 [4.594590, 5.679012,  3.454376    ],
 [2.855965, 6.542070,  2.333721    ],
 [2.710526, 6.996624,  3.166684    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
GLH.backbone = np.array([0,2,-2])
nh = GLH.coordinates[1] - GLH.coordinates[0]
can = GLH.coordinates[2] - GLH.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
GLH.nattach = -amidebond*bisector/np.linalg.norm(bisector) + GLH.coordinates[0]
GLH.cattach = np.array([6.204455,3.702003,-2.420137E-06])
GLH.chi1 = np.array([
    GLH.find("N"), GLH.find("CA"), GLH.find ("CB"), GLH.find("CG") ])
GLH.chi2 = np.array([
    GLH.find("CA"), GLH.find("CB"), GLH.find ("CG"), GLH.find("CD") ])
GLH.charges = np.array([
-0.415700,
 0.271900,
 0.014500,
 0.077900,
-0.007100,
 0.025600,
 0.025600,
-0.017400,
 0.043000,
 0.043000,
 0.680100,
-0.583800,
-0.651100,
 0.464100,
 0.597300,
-0.567900])

GLN = Template()
GLN.name = "GLN"
GLN.elements = np.array([
["N", "N"   ],
["H", "H"   ],
["C", "CA"  ],
["H", "HA"  ],
["C", "CB"  ],
["H", "HB2" ],
["H", "HB3" ],
["C", "CG"  ],
["H", "HG2" ],
["H", "HG3" ],
["C", "CD"  ],
["O", "OE1" ],
["N", "NE"  ],
["H", "HE21"],
["H", "HE22"],
["C", "C"   ],
["O", "O"   ],
])
GLN.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.274186, 5.009602,  1.194577    ],
 [5.354271, 4.863178,  1.185788    ],
 [3.973781, 5.548460,  0.295972    ],
 [3.906976, 5.848443,  2.410302    ],
 [3.138962, 5.408349,  3.262893    ],
 [4.458856, 7.061523,  2.488333    ],
 [4.248434, 7.659045,  3.274966    ],
 [5.084281, 7.376210,  1.760379    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
GLN.backbone = np.array([0,2,-2])
nh = GLN.coordinates[1] - GLN.coordinates[0]
can = GLN.coordinates[2] - GLN.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
GLN.nattach = -amidebond*bisector/np.linalg.norm(bisector) + GLN.coordinates[0]
GLN.cattach = np.array([6.204455,3.702003,-2.420137E-06])
GLN.chi1 = np.array([
    GLN.find("N"), GLN.find("CA"), GLN.find ("CB"), GLN.find("CG") ])
GLN.chi2 = np.array([
    GLN.find("CA"), GLN.find("CB"), GLN.find ("CG"), GLN.find("CD") ])
GLN.charges = np.array([
-0.415700,
 0.271900,
-0.003100,
 0.085000,
-0.003600,
 0.017100,
 0.017100,
-0.064500,
 0.035200,
 0.035200,
 0.695100,
-0.608600,
-0.940700,
 0.425100,
 0.425100,
 0.597300,
-0.567900])

GLU = Template()
GLU.name = "GLU"
GLU.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["C", "CD" ],
["O", "OE1"],
["O", "OE2"],
["C", "C"  ],
["O", "O"  ]
])
GLU.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.267328, 4.996267,  1.194946    ],
 [5.347413, 4.849843,  1.186158    ],
 [3.966923, 5.535124,  0.296342    ],
 [3.873732, 5.805369,  2.428706    ],
 [4.594590, 5.679012,  3.454376    ],
 [2.855965, 6.542070,  2.333721    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
GLU.backbone = np.array([0,2,-2])
nh = GLU.coordinates[1] - GLU.coordinates[0]
can = GLU.coordinates[2] - GLU.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
GLU.nattach = -amidebond*bisector/np.linalg.norm(bisector) + GLU.coordinates[0]
GLU.cattach = np.array([6.204455,3.702003,-2.420137E-06])
GLU.chi1 = np.array([
    GLU.find("N"), GLU.find("CA"), GLU.find ("CB"), GLU.find("CG") ])
GLU.chi2 = np.array([
    GLU.find("CA"), GLU.find("CB"), GLU.find ("CG"), GLU.find("CD") ])
GLU.charges = np.array([
-0.516300,
 0.293600,
 0.039700,
 0.110500,
 0.056000,
-0.017300,
-0.017300,
 0.013600,
-0.042500,
-0.042500,
 0.805400,
-0.818800,
-0.818800,
 0.536600,
-0.581900])

GLY = Template()
GLY.name = "GLY"
GLY.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA2"],
["H", "HA3"],
["C", "C"  ],
["O", "O"  ]
])
GLY.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.671668, 3.400125,  0.889824    ],
 [5.483710, 2.686702, -4.438857E-06],
 [5.993369, 1.568360, -8.469843E-06]
])
GLY.backbone = np.array([0,2,-2])
nh = GLY.coordinates[1] - GLY.coordinates[0]
can = GLY.coordinates[2] - GLY.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
GLY.nattach = -amidebond*bisector/np.linalg.norm(bisector) + GLY.coordinates[0]
GLY.cattach = np.array([6.204455,3.702003,-2.420137E-06])
GLY.charges = np.array([
-0.415700,
 0.271900,
-0.025200,
 0.069800,
 0.069800,
 0.597300,
-0.567900])

HID = Template()
HID.name = "HID"
HID.elements = np.array([
["N","N"  ],
["H","H"  ],
["C","CA" ],
["H","HA" ],
["C","CB" ],
["H","HB2"],
["H","HB3"],
["C","CG" ],
["N","ND1"],
["H","HD1"],
["C","CE1"],
["H","HE1"],
["N","NE2"],
["C","CD2"],
["H","HD2"],
["C","C"  ],
["O","O"  ]
])
HID.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.200813, 5.026064,  1.321087    ],
 [3.942782, 5.885086,  2.382972    ],
 [3.339725, 5.691913,  3.169805    ],
 [4.624274, 6.997642,  2.182500    ],
 [4.563048, 7.811875,  2.904563    ],
 [5.294011, 6.891451,  1.061663    ],
 [5.058974, 5.678868,  0.492453    ],
 [5.537741, 5.417846, -0.451343    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
HID.backbone = np.array([0,2,-2])
nh = HID.coordinates[1] - HID.coordinates[0]
can = HID.coordinates[2] - HID.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
HID.nattach = -amidebond*bisector/np.linalg.norm(bisector) + HID.coordinates[0]
HID.cattach = np.array([6.204455,3.702003,-2.420137E-06])
HID.chi1 = np.array([
    HID.find("N"), HID.find("CA"), HID.find ("CB"), HID.find("CG") ])
HID.chi2 = np.array([
    HID.find("CA"), HID.find("CB"), HID.find ("CG"), HID.find("ND1") ])
HID.charges = np.array([
-0.415700,
 0.271900,
 0.018800,
 0.088100,
-0.046200,
 0.040200,
 0.040200,
-0.026600,
-0.381100,
 0.364900,
 0.205700,
 0.139200,
-0.572700,
 0.129200,
 0.114700,
 0.597300,
-0.567900])


HIE = Template()
HIE.name = "HIE"
HIE.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["N", "ND1"],
["C", "CE1"],
["H", "HE1"],
["N", "NE2"],
["H", "HE2"],
["C", "CD2"],
["H", "HD2"],
["C", "C"  ],
["O", "O"  ] 
])
HIE.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.200813, 5.026064,  1.321087    ],
 [3.942782, 5.885086,  2.382972    ],
 [4.624274, 6.997642,  2.182500    ],
 [4.563048, 7.811875,  2.904563    ],
 [5.294011, 6.891451,  1.061663    ],
 [5.896297, 7.605085,  0.676854    ],
 [5.058974, 5.678868,  0.492453    ],
 [5.537741, 5.417846, -0.451343    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
HIE.backbone = np.array([0,2,-2])
nh = HIE.coordinates[1] - HIE.coordinates[0]
can = HIE.coordinates[2] - HIE.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
HIE.nattach = -amidebond*bisector/np.linalg.norm(bisector) + HIE.coordinates[0]
HIE.cattach = np.array([6.204455,3.702003,-2.420137E-06])
HIE.chi1 = np.array([
    HIE.find("N"), HIE.find("CA"), HIE.find ("CB"), HIE.find("CG") ])
HIE.chi2 = np.array([
    HIE.find("CA"), HIE.find("CB"), HIE.find ("CG"), HIE.find("ND1") ])
HIE.charges = np.array([
-0.415700,
 0.271900,
-0.058100,
 0.136000,
-0.007400,
 0.036700,
 0.036700,
 0.186800,
-0.543200,
 0.163500,
 0.143500,
-0.279500,
 0.333900,
-0.220700,
 0.186200,
 0.597300,
-0.567900])


HIP = Template()
HIP.name = "HIP"
HIP.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["N", "ND1"],
["H", "HD1"],
["C", "CE1"],
["H", "HE1"],
["N", "NE2"],
["H", "HE2"],
["C", "CD2"],
["H", "HD2"],
["C", "C"  ],
["O", "O"  ]
])
HIP.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.200813, 5.026064,  1.321087    ],
 [3.942782, 5.885086,  2.382972    ],
 [3.339725, 5.691913,  3.169805    ],
 [4.624274, 6.997642,  2.182500    ],
 [4.563048, 7.811875,  2.904563    ],
 [5.294011, 6.891451,  1.061663    ],
 [5.896297, 7.605085,  0.676854    ],
 [5.058974, 5.678868,  0.492453    ],
 [5.537741, 5.417846, -0.451343    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
HIP.backbone = np.array([0,2,-2])
nh = HIP.coordinates[1] - HIP.coordinates[0]
can = HIP.coordinates[2] - HIP.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
HIP.nattach = -amidebond*bisector/np.linalg.norm(bisector) + HIP.coordinates[0]
HIP.cattach = np.array([6.204455,3.702003,-2.420137E-06])
HIP.chi1 = np.array([
    HIP.find("N"), HIP.find("CA"), HIP.find ("CB"), HIP.find("CG") ])
HIP.chi2 = np.array([
    HIP.find("CA"), HIP.find("CB"), HIP.find ("CG"), HIP.find("ND1") ])
HIP.charges = np.array([
-0.347900,
 0.274700,
-0.135400,
 0.121200,
-0.041400,
 0.081000,
 0.081000,
-0.001200,
-0.151300,
 0.386600,
-0.017000,
 0.268100,
-0.171800,
 0.391100,
-0.114100,
 0.231700,
 0.734100,
-0.589400])

ILE = Template()
ILE.name = "ILE"
ILE.elements = np.array([
["N", "N"   ],
["H", "H"   ],
["C", "CA"  ],
["H", "HA"  ],
["C", "CB"  ],
["H", "HB"  ],
["C", "CG2" ],
["H", "HG21"],
["H", "HG22"],
["H", "HG23"],
["C", "CG1" ],
["H", "HG12"],
["H", "HG13"],
["C", "CD1" ],
["H", "HD11"],
["H", "HD12"],
["H", "HD13"],
["C", "C"   ],
["O", "O"   ]
])
ILE.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.552136, 3.620733,  1.245168    ],
 [2.470128, 3.752486,  1.245640    ],
 [3.970045, 2.845728,  2.490296    ],
 [5.052053, 2.713974,  2.490763    ],
 [3.671561, 3.399208,  3.380615    ],
 [3.485650, 1.869275,  2.490737    ],
 [4.230204, 4.986694,  1.245169    ],
 [5.312310, 4.855746,  1.245164    ],
 [3.931820, 5.541027,  0.355348    ],
 [3.812294, 5.761632,  2.490339    ],
 [4.110777, 5.208104,  3.380628    ],
 [4.296689, 6.738085,  2.490833    ],
 [2.730286, 5.893383,  2.490813    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
ILE.backbone = np.array([0,2,-2])
nh = ILE.coordinates[1] - ILE.coordinates[0]
can = ILE.coordinates[2] - ILE.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
ILE.nattach = -amidebond*bisector/np.linalg.norm(bisector) + ILE.coordinates[0]
ILE.cattach = np.array([6.204455,3.702003,-2.420137E-06])
ILE.chi1 = np.array([
    ILE.find("N"), ILE.find("CA"), ILE.find ("CB"), ILE.find("CG1") ])
ILE.chi2 = np.array([
    ILE.find("CA"), ILE.find("CB"), ILE.find ("CG1"), ILE.find("CD1") ])
ILE.charges = np.array([
-0.415700,
 0.271900,
-0.059700,
 0.086900,
 0.130300,
 0.018700,
-0.320400,
 0.088200,
 0.088200,
 0.088200,
-0.043000,
 0.023600,
 0.023600,
-0.066000,
 0.018600,
 0.018600,
 0.018600,
 0.597300,
-0.567900
])

LEU = Template()
LEU.name = "LEU"
LEU.elements = np.array([
["N", "N"   ],
["H", "H"   ],
["C", "CA"  ],
["H", "HA"  ],
["C", "CB"  ],
["H", "HB2" ],
["H", "HB3" ],
["C", "CG"  ],
["H", "HG"  ],
["C", "CD1" ],
["H", "HD11"],
["H", "HD12"],
["H", "HD13"],
["C", "CD2" ],
["H", "HD21"],
["H", "HD22"],
["H", "HD23"],
["C", "C"   ],
["O", "O"   ]
])
LEU.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.274186, 5.009602,  1.194577    ],
 [5.354271, 4.863178,  1.185788    ],
 [3.853429, 5.762895, -0.062857    ],
 [2.773449, 5.910113, -0.054557    ],
 [4.351513, 6.732052, -0.090203    ],
 [4.134159, 5.185704, -0.943846    ],
 [3.881105, 5.817645,  2.426721    ],
 [4.181626, 5.279602,  3.325774    ],
 [4.379198, 6.786825,  2.400363    ],
 [2.801135, 5.964881,  2.435959    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
LEU.backbone = np.array([0,2,-2])
nh = LEU.coordinates[1] - LEU.coordinates[0]
can = LEU.coordinates[2] - LEU.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
LEU.nattach = -amidebond*bisector/np.linalg.norm(bisector) + LEU.coordinates[0]
LEU.cattach = np.array([6.204455,3.702003,-2.420137E-06])
LEU.chi1 = np.array([
    LEU.find("N"), LEU.find("CA"), LEU.find ("CB"), LEU.find("CG") ])
LEU.chi2 = np.array([
    LEU.find("CA"), LEU.find("CB"), LEU.find ("CG"), LEU.find("CD1") ])
LEU.charges = np.array([
-0.415700,
 0.271900,
-0.051800,
 0.092200,
-0.110200,
 0.045700,
 0.045700,
 0.353100,
-0.036100,
-0.412100,
 0.100000,
 0.100000,
 0.100000,
-0.412100,
 0.100000,
 0.100000,
 0.100000,
 0.597300,
-0.567900])

LYN = Template()
LYN.name = "LYN"
LYN.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["C", "CD" ],
["H", "HD2"],
["H", "HD3"],
["C", "CE" ],
["H", "HE2"],
["H", "HE3"],
["N", "NZ" ],
["H", "HZ2"],
["H", "HZ3"],
["C", "C"  ],
["O", "O"  ]
])
LYN.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.274186, 5.009602,  1.194577    ],
 [5.354271, 4.863178,  1.185788    ],
 [3.973781, 5.548460,  0.295972    ],
 [3.881105, 5.817645,  2.426721    ],
 [2.801135, 5.964881,  2.435959    ],
 [4.181626, 5.279602,  3.325774    ],
 [4.578325, 7.173410,  2.389153    ],
 [5.658410, 7.026987,  2.380363    ],
 [4.277917, 7.712267,  1.490550    ],
 [4.199422, 7.952309,  3.576860    ],
 [4.661186, 8.850226,  3.551979    ],
 [3.198675, 8.088466,  3.584971    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
LYN.backbone = np.array([0,2,-2])
nh = LYN.coordinates[1] - LYN.coordinates[0]
can = LYN.coordinates[2] - LYN.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
LYN.nattach = -amidebond*bisector/np.linalg.norm(bisector) + LYN.coordinates[0]
LYN.cattach = np.array([6.204455,3.702003,-2.420137E-06])
LYN.chi1 = np.array([
    LYN.find("N"), LYN.find("CA"), LYN.find ("CB"), LYN.find("CG") ])
LYN.chi2 = np.array([
    GLN.find("CA"), GLN.find("CB"), GLN.find ("CG"), GLN.find("CD") ])
LYN.charges = np.array([
-0.415700,
 0.271900,
-0.072060,
 0.099400,
-0.048450,
 0.034000,
 0.034000,
 0.066120,
 0.010410,
 0.010410,
-0.037680,
 0.011550,
 0.011550,
 0.326040,
-0.033580,
-0.033580,
-1.035810,
 0.386040,
 0.386040,
 0.597300,
-0.567900])


LYS = Template()
LYS.name = "LYS"
LYS.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["C", "CD" ],
["H", "HD2"],
["H", "HD3"],
["C", "CE" ],
["H", "HE2"],
["H", "HE3"],
["N", "NZ" ],
["H", "HZ1"],
["H", "HZ2"],
["H", "HZ3"],
["C", "C"  ],
["O", "O"  ]
])
LYS.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.274186, 5.009602,  1.194577    ],
 [5.354271, 4.863178,  1.185788    ],
 [3.973781, 5.548460,  0.295972    ],
 [3.881105, 5.817645,  2.426721    ],
 [2.801135, 5.964881,  2.435959    ],
 [4.181626, 5.279602,  3.325774    ],
 [4.578325, 7.173410,  2.389153    ],
 [5.658410, 7.026987,  2.380363    ],
 [4.277917, 7.712267,  1.490550    ],
 [4.199422, 7.952309,  3.576860    ],
 [4.478085, 7.453366,  4.409628    ],
 [4.661186, 8.850226,  3.551979    ],
 [3.198675, 8.088466,  3.584971    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
LYS.backbone = np.array([0,2,-2])
nh = LYS.coordinates[1] - LYS.coordinates[0]
can = LYS.coordinates[2] - LYS.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
LYS.nattach = -amidebond*bisector/np.linalg.norm(bisector) + LYS.coordinates[0]
LYS.cattach = np.array([6.204455,3.702003,-2.420137E-06])
LYS.chi1 = np.array([
    LYS.find("N"), LYS.find("CA"), LYS.find ("CB"), LYS.find("CG") ])
LYS.chi2 = np.array([
    LYS.find("CA"), LYS.find("CB"), LYS.find ("CG"), LYS.find("CD") ])
LYS.charges = np.array([
-0.347900,
 0.274700,
-0.240000,
 0.142600,
-0.009400,
 0.036200,
 0.036200,
 0.018700,
 0.010300,
 0.010300,
-0.047900,
 0.062100,
 0.062100,
-0.014300,
 0.113500,
 0.113500,
-0.385400,
 0.340000,
 0.340000,
 0.340000,
 0.734100,
-0.589400])

MET = Template()
MET.name = "MET"
MET.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["H", "HG2"],
["H", "HG3"],
["S", "SD" ],
["C", "CE" ],
["H", "HE1"],
["H", "HE2"],
["H", "HE3"],
["C", "C"  ],
["O", "O"  ]
])
MET.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.274186, 5.009602,  1.194577    ],
 [5.354271, 4.863178,  1.185788    ],
 [3.973781, 5.548460,  0.295972    ],
 [3.817309, 5.981266,  2.651708    ],
 [4.753212, 7.463128,  2.340949    ],
 [4.433582, 7.904044,  1.396741    ],
 [4.585907, 8.175299,  3.148985    ],
 [5.814074, 7.218763,  2.286554    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
MET.backbone = np.array([0,2,-2])
nh = MET.coordinates[1] - MET.coordinates[0]
can = MET.coordinates[2] - MET.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
MET.nattach = -amidebond*bisector/np.linalg.norm(bisector) + MET.coordinates[0]
MET.cattach = np.array([6.204455,3.702003,-2.420137E-06])
MET.chi1 = np.array([
    MET.find("N"), MET.find("CA"), MET.find ("CB"), MET.find("CG") ])
MET.chi2 = np.array([
    MET.find("CA"), MET.find("CB"), MET.find ("CG"), MET.find("SD") ])
MET.charges = np.array([
-0.415700,
 0.271900,
-0.023700,
 0.088000,
 0.034200,
 0.024100,
 0.024100,
 0.001800,
 0.044000,
 0.044000,
-0.273700,
-0.053600,
 0.068400,
 0.068400,
 0.068400,
 0.597300,
-0.567900])


NHE = Template()
NHE.name = "NHE"
NHE.elements = np.array([
["N", "N"  ],
["H", "HN1"],
["H", "HN2"]
])
NHE.coordinates = np.array([
 [2.193696, 1.597759, -1.607204E-06],
 [3.034948, 1.038834, -2.739928E-06],
 [2.250077, 2.606184, -5.030662E-07]
])
nh1 = NHE.coordinates[1] - NHE.coordinates[0]
nh2 = NHE.coordinates[2] - NHE.coordinates[0]
bisector = np.linalg.norm(nh1)*nh2 + np.linalg.norm(nh2)*nh1
NHE.nattach = -amidebond*bisector/np.linalg.norm(bisector) + NHE.coordinates[0]
NHE.charges = np.array([
-0.463000,
 0.231500,
 0.231500])

NME = Template()
NME.name = "NME"
NME.elements = np.array([
["N", "N"   ],
["H", "H"   ],
["C", "CH3" ],
["H", "HH31"],
["H", "HH32"],
["H", "HH33"]
])
NME.coordinates = np.array([
 [3.325770, 1.547909, -1.607257E-06],
 [3.909407, 0.723611, -2.739904E-06],
 [3.970048, 2.845795, -1.312144E-07],
 [3.211504, 3.628554,  2.347971E-06],
 [4.591993, 2.943271,  0.889822    ],
 [4.591988, 2.943275, -0.889825    ]
])
NME.backbone = np.array([0,2,3])
nh = NME.coordinates[1] - NME.coordinates[0]
nc = NME.coordinates[2] - NME.coordinates[0]
bisector = np.linalg.norm(nh)*nc + np.linalg.norm(nc)*nh
NME.nattach = -amidebond*bisector/np.linalg.norm(bisector) + NME.coordinates[0]
NME.charges = np.array([
-0.415700,
 0.271900,
-0.149000,
 0.097600,
 0.097600,
 0.097600])

PHE = Template()
PHE.name = "PHE"
PHE.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["C", "CD1"],
["H", "HD1"],
["C", "CE1"],
["H", "HE1"],
["C", "CZ" ],
["H", "HZ" ],
["C", "CE2"],
["H", "HE2"],
["C", "CD2"],
["H", "HD2"],
["C", "C"  ],
["O", "O"  ]
])
PHE.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.200813, 5.026064,  1.321087    ],
 [3.911613, 5.857250,  2.409890    ],
 [3.236123, 5.513843,  3.193398    ],
 [4.490014, 7.129513,  2.492354    ],
 [4.264853, 7.776651,  3.340066    ],
 [5.357616, 7.570591,  1.486016    ],
 [5.807943, 8.561138,  1.550220    ],
 [5.646818, 6.739407,  0.397211    ],
 [6.322309, 7.082817, -0.386295    ],
 [5.068419, 5.467143,  0.314744    ],
 [5.293584, 4.820007, -0.532968    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
PHE.backbone = np.array([0,2,-2])
nh = PHE.coordinates[1] - PHE.coordinates[0]
can = PHE.coordinates[2] - PHE.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
PHE.nattach = -amidebond*bisector/np.linalg.norm(bisector) + PHE.coordinates[0]
PHE.cattach = np.array([6.204455,3.702003,-2.420137E-06])
PHE.chi1 = np.array([
    PHE.find("N"), PHE.find("CA"), PHE.find ("CB"), PHE.find("CG") ])
PHE.chi2 = np.array([
    PHE.find("CA"), PHE.find("CB"), PHE.find ("CG"), PHE.find("CD1") ])
PHE.charges = np.array([
-0.415700,
 0.271900,
-0.002400,
 0.097800,
-0.034300,
 0.029500,
 0.029500,
 0.011800,
-0.125600,
 0.133000,
-0.170400,
 0.143000,
-0.107200,
 0.129700,
-0.170400,
 0.143000,
-0.125600,
 0.133000,
 0.597300,
-0.567900])

PRO = Template()
PRO.name = "PRO"
PRO.elements = np.array([
["N", "N"  ],
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
["H", "HA" ],
["C", "C"  ],
["O", "O"  ]
])
PRO.coordinates = np.array([
 [3.326834,  1.557389, -1.603945E-06],
 [4.302147,  0.476598,  0.080119    ],
 [4.419998,  0.019283, -0.902263    ],
 [3.955888, -0.274040,  0.790574    ],
 [5.547126,  1.172441,  0.544693    ],
 [6.413549,  0.741636,  0.042879    ],
 [5.652950,  1.047934,  1.622376    ],
 [5.369091,  2.628184,  0.185227    ],
 [5.969289,  2.861861, -0.694123    ],
 [5.690642,  3.251038,  1.019947    ],
 [3.933610,  2.871277, -0.104508    ],
 [3.611470,  3.488570,  0.734106    ],
 [3.505164,  3.526392, -1.409783    ],
 [2.754240,  2.939065, -2.185412    ]
])
PRO.backbone = np.array([0,-4,-2])
PRO.cattach = np.array([3.904907,4.650696,-1.704043])
cdn = PRO.coordinates[1] - PRO.coordinates[0]
can = PRO.coordinates[10] - PRO.coordinates[0]
bisector = np.linalg.norm(cdn)*can + np.linalg.norm(can)*cdn
PRO.nattach = -amidebond*bisector/np.linalg.norm(bisector) + PRO.coordinates[0]
PRO.charges = np.array([
-0.254800,
 0.019200,
 0.039100,
 0.039100,
 0.018900,
 0.021300,
 0.021300,
-0.007000,
 0.025300,
 0.025300,
-0.026600,
 0.064100,
 0.589600,
-0.574800])

SER = Template()
SER.name = "SER"
SER.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["O", "OG" ],
["H", "HG" ],
["C", "C"  ],
["O", "O"  ]
])
SER.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.230753, 4.925145,  1.196917    ],
 [3.983305, 5.433814,  1.972562    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
SER.backbone = np.array([0,2,-2])
nh = SER.coordinates[1] - SER.coordinates[0]
can = SER.coordinates[2] - SER.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
SER.nattach = -amidebond*bisector/np.linalg.norm(bisector) + SER.coordinates[0]
SER.cattach = np.array([6.204455,3.702003,-2.420137E-06])
SER.chi1 = np.array([
    SER.find("N"), SER.find("CA"), SER.find ("CB"), SER.find("OG") ])
SER.chi2 = None
SER.charges = np.array([
-0.415700,
 0.271900,
-0.024900,
 0.084300,
 0.211700,
 0.035200,
 0.035200,
-0.654600,
 0.427500,
 0.597300,
-0.567900])

THR = Template()
THR.name = "THR"
THR.elements = np.array([
["N", "N"   ],
["H", "H"   ],
["C", "CA"  ],
["H", "HA"  ],
["C", "CB"  ],
["H", "HB"  ],
["C", "CG2" ],
["H", "HG21"],
["H", "HG22"],
["H", "HG23"],
["O", "OG1" ],
["H", "HG1" ],
["C", "C"   ],
["O", "O"   ]
])
THR.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [4.075059, 4.623017,  1.205786    ],
 [2.065936, 3.859425,  1.244383    ],
 [1.567127, 2.890627,  1.271209    ],
 [1.784431, 4.436953,  2.124903    ],
 [1.764699, 4.397847,  0.345796    ],
 [3.971501, 2.947413,  2.411212    ],
 [3.724052, 3.456082,  3.186857    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
THR.backbone = np.array([0,2,-2])
nh = THR.coordinates[1] - THR.coordinates[0]
can = THR.coordinates[2] - THR.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
THR.nattach = -amidebond*bisector/np.linalg.norm(bisector) + THR.coordinates[0]
THR.cattach = np.array([6.204455,3.702003,-2.420137E-06])
THR.chi1 = np.array([
    THR.find("N"), THR.find("CA"), THR.find ("CB"), THR.find("CG2") ])
THR.chi2 = None
THR.charges = np.array([
-0.415700,
 0.271900,
-0.038900,
 0.100700,
 0.365400,
 0.004300,
-0.243800,
 0.064200,
 0.064200,
 0.064200,
-0.676100,
 0.410200,
 0.597300,
-0.567900])

TRP = Template()
TRP.name = "TRP"
TRP.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["C", "CD1"],
["H", "HD1"],
["N", "NE1"],
["H", "HE1"],
["C", "CE2"],
["C", "CZ2"],
["H", "HZ2"],
["C", "CH2"],
["H", "HH2"],
["C", "CZ3"],
["H", "HZ3"],
["C", "CE3"],
["H", "HE3"],
["C", "CD2"],
["C", "C"  ],
["O", "O"  ]
])
TRP.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.200813, 5.026064,  1.321087    ],
 [4.023453, 5.931084,  2.293240    ],
 [3.368841, 5.705466,  3.135071    ],
 [4.811943, 7.073555,  1.949808    ],
 [4.882921, 7.922010,  2.493118    ],
 [5.427347, 6.842060,  0.816764    ],
 [6.297161, 7.689052,  0.119605    ],
 [6.531230, 8.676649,  0.517050    ],
 [6.814091, 7.187011, -1.069023    ],
 [7.498074, 7.791857, -1.664362    ],
 [6.482659, 5.953119, -1.505101    ],
 [6.897660, 5.575648, -2.439654    ],
 [5.604041, 5.117355, -0.785636    ],
 [5.358720, 4.126570, -1.168080    ],
 [5.083390, 5.623004,  0.411545    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
TRP.backbone = np.array([0,2,-2])
nh = TRP.coordinates[1] - TRP.coordinates[0]
can = TRP.coordinates[2] - TRP.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
TRP.nattach = -amidebond*bisector/np.linalg.norm(bisector) + TRP.coordinates[0]
TRP.cattach = np.array([6.204455,3.702003,-2.420137E-06])
TRP.chi1 = np.array([
    TRP.find("N"), TRP.find("CA"), TRP.find ("CB"), TRP.find("CG") ])
TRP.chi2 = np.array([
    TRP.find("CA"), TRP.find("CB"), TRP.find ("CG"), TRP.find("CD1") ])
TRP.charges = np.array([
-0.415700,
 0.271900,
-0.027500,
 0.112300,
-0.005000,
 0.033900,
 0.033900,
-0.141500,
-0.163800,
 0.206200,
-0.341800,
 0.341200,
 0.138000,
-0.260100,
 0.157200,
-0.113400,
 0.141700,
-0.197200,
 0.144700,
-0.238700,
 0.170000,
 0.124300,
 0.597300,
-0.567900])


TYR = Template()
TYR.name = "TYR"
TYR.elements = np.array([
["N", "N"  ],
["H", "H"  ],
["C", "CA" ],
["H", "HA" ],
["C", "CB" ],
["H", "HB2"],
["H", "HB3"],
["C", "CG" ],
["C", "CD1"],
["H", "HD1"],
["C", "CE1"],
["H", "HE1"],
["C", "CZ" ],
["O", "OH" ],
["H", "HH" ],
["C", "CE2"],
["H", "HE2"],
["C", "CD2"],
["H", "HD2"],
["C", "C"  ],
["O", "O"  ]
])
TYR.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.877484, 3.115795,  2.131197    ],
 [4.267328, 4.996267,  1.194946    ],
 [4.059927, 5.918911,  2.227280    ],
 [3.400108, 5.668218,  3.057877    ],
 [4.699998, 7.163547,  2.192791    ],
 [4.538522, 7.881891,  2.996538    ],
 [5.547471, 7.485542,  1.125970    ],
 [6.169255, 8.694617,  1.092468    ],
 [5.956327, 9.246984,  1.848214    ],
 [5.754875, 6.562900,  0.093635    ],
 [6.414694, 6.813595, -0.736962    ],
 [5.114806, 5.318263,  0.128119    ],
 [5.276286, 4.599920, -0.675627    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
TYR.backbone = np.array([0,2,-2])
nh = TYR.coordinates[1] - TYR.coordinates[0]
can = TYR.coordinates[2] - TYR.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
TYR.nattach = -amidebond*bisector/np.linalg.norm(bisector) + TYR.coordinates[0]
TYR.cattach = np.array([6.204455,3.702003,-2.420137E-06])
TYR.chi1 = np.array([
    TYR.find("N"), TYR.find("CA"), TYR.find ("CB"), TYR.find("CG") ])
TYR.chi2 = np.array([
    TYR.find("CA"), TYR.find("CB"), TYR.find ("CG"), TYR.find("CD1") ])
TYR.charges = np.array([
-0.415700,
 0.271900,
-0.001400,
 0.087600,
-0.015200,
 0.029500,
 0.029500,
-0.001100,
-0.190600,
 0.169900,
-0.234100,
 0.165600,
 0.322600,
-0.557900,
 0.399200,
-0.234100,
 0.165600,
-0.190600,
 0.169900,
 0.597300,
-0.567900])


VAL = Template()
VAL.name = "VAL"
VAL.elements = np.array([
["N", "N"   ],
["H", "H"   ],
["C", "CA"  ],
["H", "HA"  ],
["C", "CB"  ],
["H", "HB"  ],
["C", "CG1" ],
["H", "HG11"],
["H", "HG12"],
["H", "HG13"],
["C", "CG2" ],
["H", "HG21"],
["H", "HG22"],
["H", "HG23"],
["C", "C"   ],
["O", "O"   ]
])
VAL.coordinates = np.array([
 [3.325770, 1.547909, -1.607204E-06],
 [3.909407, 0.723611, -2.739882E-06],
 [3.970048, 2.845795, -1.311163E-07],
 [3.671663, 3.400129, -0.889820    ],
 [3.576965, 3.653838,  1.232143    ],
 [2.496995, 3.801075,  1.241379    ],
 [3.997712, 2.900483,  2.489542    ],
 [5.077693, 2.753265,  2.481244    ],
 [3.716972, 3.477628,  3.370558    ],
 [3.499630, 1.931323,  2.516834    ],
 [4.274186, 5.009602,  1.194577    ],
 [3.973781, 5.548460,  0.295972    ],
 [3.993559, 5.587585,  2.075079    ],
 [5.354271, 4.863178,  1.185788    ],
 [5.485541, 2.705207, -4.398755E-06],
 [6.008824, 1.593175, -8.449768E-06]
])
VAL.backbone = np.array([0,2,-2])
nh = VAL.coordinates[1] - VAL.coordinates[0]
can = VAL.coordinates[2] - VAL.coordinates[0]
bisector = np.linalg.norm(nh)*can + np.linalg.norm(can)*nh
VAL.nattach = -amidebond*bisector/np.linalg.norm(bisector) + VAL.coordinates[0]
VAL.cattach = np.array([6.204455,3.702003,-2.420137E-06])
VAL.chi1 = np.array([
    VAL.find("N"), VAL.find("CA"), VAL.find ("CB"), VAL.find("CG1") ])
VAL.chi2 = None
VAL.charges = np.array([
-0.415700,
 0.271900,
-0.087500,
 0.096900,
 0.298500,
-0.029700,
-0.319200,
 0.079100,
 0.079100,
 0.079100,
-0.319200,
 0.079100,
 0.079100,
 0.079100,
 0.597300,
-0.567900])
