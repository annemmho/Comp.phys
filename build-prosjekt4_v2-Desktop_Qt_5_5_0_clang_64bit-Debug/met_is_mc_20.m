mc = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000, 210000, 220000, 230000, 240000, 250000, 260000, 270000, 280000, 290000, 300000, 310000, 320000, 330000, 340000, 350000, 360000, 370000, 380000, 390000, 400000, 410000, 420000, 430000, 440000, 450000, 460000, 470000, 480000, 490000, 500000, 510000, 520000, 530000, 540000, 550000, 560000, 570000, 580000, 590000, 600000, 610000, 620000, 630000, 640000, 650000, 660000, 670000, 680000, 690000, 700000, 710000, 720000, 730000, 740000, 750000, 760000, 770000, 780000, 790000, 800000, 810000, 820000, 830000, 840000, 850000, 860000, 870000, 880000, 890000, 900000, 910000, 920000, 930000, 940000, 950000, 960000, 970000, 980000, 990000, 1e+06, ];
M_var = [5.3878, 2.72738, 1.82607, 1.37266, 1.09975, 0.917447, 0.787051, 0.689161, 0.612969, 0.551985, 0.502057, 0.460429, 0.425198, 0.394989, 0.368803, 0.345894, 0.325668, 0.307689, 0.291602, 0.27712, 0.264023, 0.252107, 0.241226, 0.231252, 0.222072, 0.213601, 0.205755, 0.19847, 0.191687, 0.185352, 0.179429, 0.173877, 0.168659, 0.163746, 0.159113, 0.154742, 0.150607, 0.146687, 0.142972, 0.139435, 0.136076, 0.132873, 0.129822, 0.126906, 0.124119, 0.121457, 0.118906, 0.116459, 0.114116, 0.111867, 0.109707, 0.107631, 0.10563, 0.103704, 0.101846, 0.100056, 0.0983289, 0.0966621, 0.0950504, 0.0934943, 0.0919876, 0.0905299, 0.089117, 0.0877493, 0.0864245, 0.0851388, 0.0838914, 0.0826804, 0.0815061, 0.0803642, 0.079256, 0.0781772, 0.0771278, 0.0761082, 0.0751152, 0.0741482, 0.0732055, 0.0722889, 0.0713948, 0.0705219, 0.0696704, 0.0688411, 0.0680311, 0.0672418, 0.0664684, 0.0657132, 0.0649767, 0.0642564, 0.0635532, 0.0628639, 0.0621898, 0.0615314, 0.0608879, 0.0602562, 0.0596397, 0.0590356, 0.0584442, 0.0578651, 0.0572971, 0.0567391, ];
E_var = [0.81232, 0.422723, 0.290245, 0.224004, 0.184086, 0.157307, 0.138244, 0.124287, 0.113229, 0.104497, 0.0970852, 0.09088, 0.0856508, 0.0811319, 0.077285, 0.0740337, 0.0710432, 0.0684338, 0.0661314, 0.0640046, 0.0622141, 0.0604973, 0.0588951, 0.0574532, 0.0560424, 0.0548202, 0.0536534, 0.0526049, 0.0516036, 0.050628, 0.0497452, 0.0489293, 0.0481647, 0.0474185, 0.0466972, 0.046085, 0.0454973, 0.0449195, 0.0444154, 0.0438281, 0.0433513, 0.0428434, 0.0424052, 0.0419202, 0.0414546, 0.0410672, 0.0406689, 0.0402575, 0.0399217, 0.0395923, 0.039296, 0.0390314, 0.0387349, 0.0384383, 0.0381431, 0.0378673, 0.037611, 0.0373867, 0.0371434, 0.0369304, 0.0367048, 0.0364904, 0.0362619, 0.036061, 0.0358754, 0.0356777, 0.0354863, 0.0353, 0.0351397, 0.0349678, 0.0348161, 0.0346562, 0.0345019, 0.0343763, 0.0342395, 0.0341055, 0.033956, 0.0338393, 0.0337152, 0.0335777, 0.0334441, 0.0333394, 0.0332214, 0.0331281, 0.0330003, 0.0328798, 0.0327832, 0.0326792, 0.0325829, 0.0324717, 0.0323606, 0.0322777, 0.0321996, 0.0320973, 0.0320134, 0.0319318, 0.0318548, 0.0317833, 0.0317034, 0.0316095, ];
MeanM = [0.9812, 0.990232, 0.993253, 0.994756, 0.995662, 0.996268, 0.996697, 0.997008, 0.997257, 0.997452, 0.997621, 0.997761, 0.997877, 0.99798, 0.998067, 0.998141, 0.998209, 0.998266, 0.998318, 0.998366, 0.998406, 0.998444, 0.99848, 0.998511, 0.998544, 0.998571, 0.998598, 0.998621, 0.998644, 0.998666, 0.998686, 0.998705, 0.998722, 0.998738, 0.998755, 0.998768, 0.998782, 0.998795, 0.998806, 0.99882, 0.99883, 0.998842, 0.998852, 0.998862, 0.998873, 0.998881, 0.99889, 0.9989, 0.998907, 0.998914, 0.998921, 0.998927, 0.998933, 0.99894, 0.998947, 0.998953, 0.998958, 0.998963, 0.998969, 0.998973, 0.998979, 0.998983, 0.998988, 0.998993, 0.998997, 0.999001, 0.999005, 0.99901, 0.999013, 0.999017, 0.999021, 0.999024, 0.999028, 0.99903, 0.999033, 0.999036, 0.99904, 0.999042, 0.999045, 0.999048, 0.999051, 0.999053, 0.999056, 0.999058, 0.999061, 0.999064, 0.999066, 0.999068, 0.999071, 0.999073, 0.999076, 0.999078, 0.999079, 0.999081, 0.999083, 0.999085, 0.999086, 0.999088, 0.999089, 0.999092, ];
Probability = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -792, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -792, -792, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -792, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -792, -792, -800, -792, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -800, -792, -800, -800, -800, -800, -800, -800, ];
E_avg = [-1.99052, -1.99383, -1.99497, -1.9955, -1.99584, -1.99607, -1.99623, -1.9963, -1.99638, -1.99644, -1.99652, -1.99657, -1.99662, -1.99667, -1.99671, -1.99673, -1.99676, -1.99677, -1.99679, -1.99681, -1.99681, -1.99682, -1.99684, -1.99684, -1.99687, -1.99687, -1.99688, -1.99689, -1.9969, -1.99691, -1.99692, -1.99693, -1.99694, -1.99694, -1.99695, -1.99696, -1.99696, -1.99697, -1.99697, -1.99698, -1.99698, -1.99699, -1.99699, -1.997, -1.99701, -1.99701, -1.99702, -1.99702, -1.99703, -1.99703, -1.99703, -1.99703, -1.99703, -1.99703, -1.99704, -1.99704, -1.99704, -1.99704, -1.99705, -1.99705, -1.99705, -1.99705, -1.99705, -1.99706, -1.99705, -1.99706, -1.99706, -1.99706, -1.99706, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99707, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, -1.99708, ];
accepted_tot = [0.8071, 0.14575, 0.0933, 0.07235, 0.05916, 0.0479667, 0.0420857, 0.03805, 0.0318667, 0.03002, 0.0260364, 0.0238833, 0.0217538, 0.0197571, 0.01916, 0.0178625, 0.0166235, 0.0171667, 0.0147368, 0.01451, 0.0146571, 0.0134773, 0.0128217, 0.0126917, 0.010872, 0.0112692, 0.010763, 0.0102071, 0.00994483, 0.00920667, 0.00929032, 0.009275, 0.0086, 0.00837059, 0.00802286, 0.00831389, 0.00762162, 0.00746053, 0.00755897, 0.00644, 0.00720488, 0.00657619, 0.00686512, 0.00630455, 0.00593333, 0.00623043, 0.006, 0.00585833, 0.00580816, 0.005672, 0.00566667, 0.00577692, 0.00544151, 0.00533704, 0.00546909, 0.00525, 0.00506667, 0.00508103, 0.00494407, 0.00493, 0.00467869, 0.00468065, 0.00450794, 0.00461875, 0.00455385, 0.00440909, 0.00437313, 0.00412059, 0.00426087, 0.00398714, 0.0042, 0.00394861, 0.00376301, 0.00407703, 0.00394667, 0.00398158, 0.00368831, 0.00384872, 0.00373924, 0.0035175, 0.00359259, 0.00356341, 0.00346024, 0.00341429, 0.00336706, 0.00326744, 0.00332644, 0.00350455, 0.00322022, 0.00323111, 0.00308132, 0.00320652, 0.00317849, 0.00294362, 0.00320316, 0.00312083, 0.00309897, 0.00310816, 0.00294141, 0.002796, ];

plot(mc, E_avg)
