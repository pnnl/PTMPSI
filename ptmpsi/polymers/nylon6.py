from ptmpsi.polymers import UnitCell

models = {}
models["Malta1979"] = UnitCell(
  9.71, 8.19, 17.24, 90.0, 90.0, 112.5, 
   np.array([
   ["C", -0.052900,-0.003500,0.000000],
   ["C",  0.047000, 0.003200,0.069700],
   ["C", -0.036600,-0.002500,0.145600],
   ["C",  0.059200, 0.004000,0.214900],
   ["C",  0.069200, 0.004600,0.354400],
   ["C", -0.030700,-0.002000,0.424200],
   ["N ", -0.008200,-0.000600,0.281300],
   ["O ",  0.189300, 0.012700,0.208700],
   ["C",  0.447100,-0.003500,0.490000],
   ["C",  0.547000, 0.003200,0.420200],
   ["C",  0.463400,-0.002500,0.344400],
   ["C",  0.559200, 0.004000,0.275100],
   ["C",  0.569200, 0.004600,0.135500],
   ["C",  0.469300,-0.002000,0.065800],
   ["N ",  0.491800,-0.000600,0.208700],
   ["O ",  0.689300, 0.012700,0.281300],
   ["C", -0.052900, 0.496500,0.214400],
   ["C",  0.047000, 0.503100,0.284100],
   ["C", -0.036600, 0.497600,0.360000],
   ["C",  0.059200, 0.504000,0.429300],
   ["C",  0.069200, 0.504600,0.568800],
   ["C", -0.030700, 0.497900,0.638500],
   ["N ", -0.008200, 0.499400,0.495600],
   ["O ",  0.189300, 0.512700,0.423100],
   ["C",  0.447100, 0.496500,0.704300],
   ["C",  0.547000, 0.503100,0.634600],
   ["C",  0.463400, 0.497600,0.558700],
   ["C",  0.559200, 0.504000,0.489400],
   ["C",  0.569200, 0.504600,0.349900],
   ["C",  0.469300, 0.497900,0.280200],
   ["N ",  0.491800, 0.499400,0.423100],
   ["O ",  0.689300, 0.512700,0.4956]
   ]),
   ref="EPJ 15, 765 (1979)")

models["Quarti"] = UnitCell(
  9.36, 6.92, 17.41, 90.0, 90.0, 111.71,
np.array([
[N,    0.0094,    0.2858,   -0.0110],
[C,   -0.0587,    0.2173,   -0.0018],
[O,   -0.1957,    0.2112,    0.0174],
[C,    0.0424,    0.1470,   -0.0196],
[C,   -0.0456,    0.0713,    0.0135],
[C,    0.0598,    0.0018,   -0.0147],
[C,   -0.0309,   -0.0732,    0.0213],
[C,    0.0758,   -0.1423,   -0.0152],
[H,    0.1181,    0.2872,   -0.0097],
[H,    0.1000,    0.1531,    0.0913],
[H,    0.1348,    0.1482,   -0.1737],
[H,   -0.1045,    0.0670,   -0.0962],
[H,   -0.1357,    0.0698,    0.1690],
[H,    0.1206,    0.0050,    0.0934],
[H,    0.1477,    0.0013,   -0.1724],
[H,   -0.0943,   -0.0758,   -0.0835],
[H,   -0.1176,   -0.0737,    0.1796],
[H,    0.1362,   -0.1413,    0.0934],
[H,    0.1635,   -0.1393,   -0.1729],
[N,   -0.4861,    0.2175,   -0.0070],
[C,    0.4392,    0.2837,   -0.0027],
[O,    0.3017,    0.2850,    0.0153],
[C,   -0.4673,    0.3566,   -0.0238],
[C,    0.4376,    0.4302,    0.0119],
[C,   -0.4583,   -0.4998,   -0.0170],
[C,    0.4536,   -0.4239,    0.0174],
[C,   -0.4358,   -0.3563,   -0.0189],
[H,   -0.3775,    0.2190,   -0.0054],
[H,   -0.4071,    0.3539,    0.0850],
[H,   -0.3768,    0.3582,   -0.1796],
[H,    0.3777,    0.4327,   -0.0971],
[H,    0.3492,    0.4308,    0.1689],
[H,   -0.3976,    0.4960,    0.0910],
[H,   -0.3690,    0.4982,   -0.1735],
[H,    0.3907,   -0.4198,   -0.0877],
[H,    0.3677,   -0.4224,    0.1760],
[H,   -0.3745,   -0.3591,    0.0882],
[H,   -0.3478,   -0.3600,   -0.1763],
[N,    0.0141,    0.4965,    0.4927],
[C,   -0.0608,    0.4305,    0.4970],
[O,   -0.1983,    0.4293,   -0.4851],
[C,    0.0323,    0.3574,    0.4760],
[C,   -0.0631,    0.2839,   -0.4881],
[C,    0.0409,    0.2139,    0.4829],
[C,   -0.0472,    0.1380,   -0.4825],
[C,    0.0636,    0.0705,    0.4813],
[H,    0.1227,    0.4949,    0.4945],
[H,    0.0925,    0.3600,   -0.4153],
[H,    0.1228,    0.3558,    0.3202],
[H,   -0.1230,    0.2815,    0.4029],
[H,   -0.1514,    0.2833,   -0.3311],
[H,    0.1016,    0.2181,   -0.4090],
[H,    0.1303,    0.2158,    0.3265],
[H,   -0.1101,    0.1338,    0.4124],
[H,   -0.1331,    0.1365,   -0.3239],
[H,    0.1249,    0.0733,   -0.4116],
[H,    0.1516,    0.0741,    0.3239],
[N,   -0.4908,    0.4287,    0.4895],
[C,    0.4414,    0.4972,    0.4984],
[O,    0.3045,   -0.4964,   -0.4825],
[C,   -0.4571,   -0.4326,    0.4804],
[C,    0.4552,   -0.3568,   -0.4868],
[C,   -0.4393,   -0.2873,    0.4848],
[C,    0.4699,   -0.2123,   -0.4795],
[C,   -0.4236,   -0.1431,    0.4841],
[H,   -0.3821,    0.4271,    0.4906],
[H,   -0.3994,   -0.4388,   -0.4089],
[H,   -0.3647,   -0.4339,    0.3262],
[H,    0.3962,   -0.3525,    0.4036],
[H,    0.3652,   -0.3552,   -0.3312],
[H,   -0.3785,   -0.2904,   -0.4070],
[H,   -0.3514,   -0.2868,    0.3271],
[H,    0.4067,   -0.2099,    0.4154],
[H,    0.3831,   -0.2119,   -0.3213],
[H,   -0.3632,   -0.1441,   -0.4073],
[H,   -0.3360,   -0.1461,    0.3264]
])
)