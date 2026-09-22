# All Parametric Survival Models
Ben Geisler
2026-09-22

- [Overview](#overview)
- [Information Criteria](#information-criteria)
- [Coefficient Tables](#coefficient-tables)
- [Ordering-Constrained Selection
  Summary](#ordering-constrained-selection-summary)

# Overview

This technical appendix lists every fitted parametric distribution for
the single joint economic survival model. Overall survival (OS) and
progression-free survival (PFS) are not selected independently. The
base-case pair is the pair with the lowest combined AIC among candidates
that preserve OS $\geq$ PFS for control and every economic biomarker
subgroup at every weekly time point over the modeled horizon.

With the current data and 10-year horizon, this ordering-constrained
procedure selects **gamma for OS and gamma for PFS**.

| Endpoint | Formula |
|:---|:---|
| Overall survival | Surv(OSwk, Death) ~ Age + sex + Rx + crp \* Rx + tmb_braf \* Rx |
| Progression-free survival | Surv(PFSwk, Progression) ~ Age + sex + Rx + crp \* Rx + tmb_braf \* Rx |

Fitted economic survival formulas

# Information Criteria

| Distribution |    AIC |    BIC |                  Endpoint | Delta_AIC | Delta_BIC |
|-------------:|-------:|-------:|--------------------------:|----------:|----------:|
|      weibull | 666.98 | 686.96 |          Overall survival |      0.00 |      0.00 |
|    weibullph | 666.98 | 686.96 |          Overall survival |      0.00 |      0.00 |
|        gamma | 667.31 | 687.28 |          Overall survival |      0.32 |      0.32 |
|       llogis | 668.86 | 688.83 |          Overall survival |      1.87 |      1.87 |
|     gengamma | 668.97 | 691.17 |          Overall survival |      1.99 |      4.21 |
|    lognormal | 670.97 | 690.95 |          Overall survival |      3.99 |      3.99 |
|         genf | 670.98 | 695.39 |          Overall survival |      3.99 |      8.43 |
|     gompertz | 673.66 | 693.63 |          Overall survival |      6.67 |      6.67 |
|  exponential | 686.14 | 703.89 |          Overall survival |     19.15 |     16.93 |
|    lognormal | 634.96 | 654.93 | Progression-free survival |      0.00 |      0.00 |
|        gamma | 636.86 | 656.83 | Progression-free survival |      1.90 |      1.90 |
|     gengamma | 636.96 | 659.15 | Progression-free survival |      2.00 |      4.22 |
|       llogis | 637.39 | 657.37 | Progression-free survival |      2.44 |      2.44 |
|      weibull | 638.25 | 658.23 | Progression-free survival |      3.29 |      3.29 |
|    weibullph | 638.25 | 658.23 | Progression-free survival |      3.29 |      3.29 |
|         genf | 638.96 | 663.37 | Progression-free survival |      4.00 |      8.44 |
|     gompertz | 640.68 | 660.66 | Progression-free survival |      5.73 |      5.73 |
|  exponential | 645.75 | 663.51 | Progression-free survival |     10.80 |      8.58 |

Information criteria for all fitted distributions

# Coefficient Tables

| Endpoint | Distribution | Parameter | Estimate | SE | 95% CI |
|:---|:---|:---|---:|---:|---:|
| OS | exponential | rate | 0.0057 | 0.0058 | \[0.0008, 0.0415\] |
| OS | exponential | Age | 0.0073 | 0.0149 | \[-0.0218, 0.0364\] |
| OS | exponential | sex1 | 0.2141 | 0.2677 | \[-0.3105, 0.7387\] |
| OS | exponential | RxExperimental arm | 0.2869 | 0.3701 | \[-0.4385, 1.0123\] |
| OS | exponential | crp1 | -0.5340 | 0.5502 | \[-1.6123, 0.5443\] |
| OS | exponential | tmb_braf1 | -0.1355 | 0.3975 | \[-0.9145, 0.6436\] |
| OS | exponential | RxExperimental arm:crp1 | -0.1318 | 0.6794 | \[-1.4634, 1.1997\] |
| OS | exponential | RxExperimental arm:tmb_braf1 | -0.1209 | 0.5541 | \[-1.2069, 0.9651\] |
| OS | weibull | shape | 1.7474 | 0.1883 | \[1.4147, 2.1583\] |
| OS | weibull | scale | 158.0336 | 94.1040 | \[49.1910, 507.7067\] |
| OS | weibull | Age | -0.0043 | 0.0088 | \[-0.0216, 0.0130\] |
| OS | weibull | sex1 | -0.1969 | 0.1556 | \[-0.5019, 0.1082\] |
| OS | weibull | RxExperimental arm | -0.2573 | 0.2113 | \[-0.6715, 0.1569\] |
| OS | weibull | crp1 | 0.4732 | 0.3132 | \[-0.1406, 1.0871\] |
| OS | weibull | tmb_braf1 | 0.0589 | 0.2291 | \[-0.3901, 0.5078\] |
| OS | weibull | RxExperimental arm:crp1 | 0.0566 | 0.3898 | \[-0.7074, 0.8206\] |
| OS | weibull | RxExperimental arm:tmb_braf1 | 0.0729 | 0.3223 | \[-0.5587, 0.7046\] |
| OS | weibullph | shape | 1.7474 | 0.1883 | \[1.4147, 2.1583\] |
| OS | weibullph | scale | 0.0001 | 0.0002 | \[0.0000, 0.0023\] |
| OS | weibullph | Age | 0.0075 | 0.0154 | \[-0.0227, 0.0378\] |
| OS | weibullph | sex1 | 0.3440 | 0.2744 | \[-0.1938, 0.8818\] |
| OS | weibullph | RxExperimental arm | 0.4496 | 0.3718 | \[-0.2791, 1.1783\] |
| OS | weibullph | crp1 | -0.8269 | 0.5524 | \[-1.9097, 0.2559\] |
| OS | weibullph | tmb_braf1 | -0.1029 | 0.4000 | \[-0.8868, 0.6810\] |
| OS | weibullph | RxExperimental arm:crp1 | -0.0989 | 0.6809 | \[-1.4333, 1.2356\] |
| OS | weibullph | RxExperimental arm:tmb_braf1 | -0.1275 | 0.5628 | \[-1.2305, 0.9756\] |
| OS | llogis | shape | 2.4310 | 0.2657 | \[1.9623, 3.0117\] |
| OS | llogis | scale | 92.4453 | 63.5361 | \[24.0363, 355.5507\] |
| OS | llogis | Age | -0.0012 | 0.0099 | \[-0.0207, 0.0183\] |
| OS | llogis | sex1 | -0.1538 | 0.1777 | \[-0.5021, 0.1946\] |
| OS | llogis | RxExperimental arm | -0.1625 | 0.2737 | \[-0.6989, 0.3739\] |
| OS | llogis | crp1 | 0.2228 | 0.4011 | \[-0.5634, 1.0089\] |
| OS | llogis | tmb_braf1 | -0.0308 | 0.2912 | \[-0.6016, 0.5399\] |
| OS | llogis | RxExperimental arm:crp1 | 0.3317 | 0.4766 | \[-0.6025, 1.2658\] |
| OS | llogis | RxExperimental arm:tmb_braf1 | 0.2056 | 0.3730 | \[-0.5254, 0.9366\] |
| OS | lognormal | meanlog | 4.8783 | 0.6945 | \[3.5170, 6.2395\] |
| OS | lognormal | sdlog | 0.7411 | 0.0701 | \[0.6157, 0.8920\] |
| OS | lognormal | Age | -0.0076 | 0.0098 | \[-0.0269, 0.0117\] |
| OS | lognormal | sex1 | -0.1874 | 0.1846 | \[-0.5491, 0.1743\] |
| OS | lognormal | RxExperimental arm | -0.1423 | 0.2827 | \[-0.6965, 0.4118\] |
| OS | lognormal | crp1 | 0.0737 | 0.3650 | \[-0.6418, 0.7891\] |
| OS | lognormal | tmb_braf1 | 0.1219 | 0.2800 | \[-0.4268, 0.6707\] |
| OS | lognormal | RxExperimental arm:crp1 | 0.5759 | 0.4598 | \[-0.3253, 1.4770\] |
| OS | lognormal | RxExperimental arm:tmb_braf1 | 0.0858 | 0.3845 | \[-0.6677, 0.8393\] |
| OS | gamma | shape | 2.4303 | 0.4236 | \[1.7270, 3.4198\] |
| OS | gamma | rate | 0.0185 | 0.0123 | \[0.0051, 0.0678\] |
| OS | gamma | Age | 0.0036 | 0.0092 | \[-0.0145, 0.0217\] |
| OS | gamma | sex1 | 0.1796 | 0.1675 | \[-0.1486, 0.5079\] |
| OS | gamma | RxExperimental arm | 0.2215 | 0.2365 | \[-0.2420, 0.6850\] |
| OS | gamma | crp1 | -0.3896 | 0.3306 | \[-1.0375, 0.2583\] |
| OS | gamma | tmb_braf1 | -0.0663 | 0.2474 | \[-0.5513, 0.4187\] |
| OS | gamma | RxExperimental arm:crp1 | -0.1711 | 0.4137 | \[-0.9820, 0.6397\] |
| OS | gamma | RxExperimental arm:tmb_braf1 | -0.0982 | 0.3451 | \[-0.7746, 0.5782\] |
| OS | gompertz | shape | 0.0100 | 0.0024 | \[0.0053, 0.0147\] |
| OS | gompertz | rate | 0.0022 | 0.0024 | \[0.0003, 0.0192\] |
| OS | gompertz | Age | 0.0095 | 0.0156 | \[-0.0210, 0.0401\] |
| OS | gompertz | sex1 | 0.3436 | 0.2753 | \[-0.1960, 0.8832\] |
| OS | gompertz | RxExperimental arm | 0.4435 | 0.3739 | \[-0.2894, 1.1764\] |
| OS | gompertz | crp1 | -0.9710 | 0.5695 | \[-2.0872, 0.1452\] |
| OS | gompertz | tmb_braf1 | -0.1364 | 0.4005 | \[-0.9214, 0.6486\] |
| OS | gompertz | RxExperimental arm:crp1 | 0.0537 | 0.6864 | \[-1.2916, 1.3990\] |
| OS | gompertz | RxExperimental arm:tmb_braf1 | -0.1193 | 0.5633 | \[-1.2235, 0.9848\] |
| OS | gengamma | mu | 4.9547 | 0.6275 | \[3.7249, 6.1846\] |
| OS | gengamma | sigma | 0.5897 | 0.1105 | \[0.4085, 0.8513\] |
| OS | gengamma | Q | 0.9047 | 0.5035 | \[-0.0822, 1.8916\] |
| OS | gengamma | Age | -0.0032 | 0.0089 | \[-0.0206, 0.0142\] |
| OS | gengamma | sex1 | -0.1927 | 0.1615 | \[-0.5092, 0.1237\] |
| OS | gengamma | RxExperimental arm | -0.2487 | 0.2217 | \[-0.6833, 0.1860\] |
| OS | gengamma | crp1 | 0.4569 | 0.3345 | \[-0.1987, 1.1124\] |
| OS | gengamma | tmb_braf1 | 0.0587 | 0.2336 | \[-0.3991, 0.5165\] |
| OS | gengamma | RxExperimental arm:crp1 | 0.0772 | 0.4234 | \[-0.7526, 0.9069\] |
| OS | gengamma | RxExperimental arm:tmb_braf1 | 0.0831 | 0.3298 | \[-0.5633, 0.7296\] |
| OS | genf | mu | 4.9592 | 0.6282 | \[3.7280, 6.1904\] |
| OS | genf | sigma | 0.5874 | 0.1114 | \[0.4050, 0.8519\] |
| OS | genf | Q | 0.9142 | 0.5117 | \[-0.0888, 1.9171\] |
| OS | genf | P | 0.0025 | 0.0733 | \[0.0000, 12109512653441758921864.0000\] |
| OS | genf | Age | -0.0032 | 0.0088 | \[-0.0204, 0.0140\] |
| OS | genf | sex1 | -0.1937 | 0.1612 | \[-0.5096, 0.1222\] |
| OS | genf | RxExperimental arm | -0.2488 | 0.2212 | \[-0.6823, 0.1846\] |
| OS | genf | crp1 | 0.4567 | 0.3339 | \[-0.1977, 1.1111\] |
| OS | genf | tmb_braf1 | 0.0585 | 0.2330 | \[-0.3981, 0.5151\] |
| OS | genf | RxExperimental arm:crp1 | 0.0752 | 0.4225 | \[-0.7530, 0.9034\] |
| OS | genf | RxExperimental arm:tmb_braf1 | 0.0822 | 0.3292 | \[-0.5629, 0.7274\] |
| PFS | exponential | rate | 0.0258 | 0.0249 | \[0.0039, 0.1714\] |
| PFS | exponential | Age | -0.0037 | 0.0144 | \[-0.0320, 0.0245\] |
| PFS | exponential | sex1 | 0.0718 | 0.2674 | \[-0.4522, 0.5959\] |
| PFS | exponential | RxExperimental arm | 0.4714 | 0.3673 | \[-0.2485, 1.1913\] |
| PFS | exponential | crp1 | -0.3244 | 0.5481 | \[-1.3986, 0.7498\] |
| PFS | exponential | tmb_braf1 | -0.2357 | 0.3960 | \[-1.0118, 0.5404\] |
| PFS | exponential | RxExperimental arm:crp1 | -0.6180 | 0.6791 | \[-1.9489, 0.7130\] |
| PFS | exponential | RxExperimental arm:tmb_braf1 | -0.2985 | 0.5660 | \[-1.4078, 0.8109\] |
| PFS | weibull | shape | 1.3913 | 0.1383 | \[1.1450, 1.6906\] |
| PFS | weibull | scale | 40.3187 | 28.7575 | \[9.9628, 163.1662\] |
| PFS | weibull | Age | 0.0040 | 0.0108 | \[-0.0171, 0.0251\] |
| PFS | weibull | sex1 | -0.0211 | 0.1976 | \[-0.4083, 0.3661\] |
| PFS | weibull | RxExperimental arm | -0.4314 | 0.2643 | \[-0.9494, 0.0865\] |
| PFS | weibull | crp1 | 0.2860 | 0.4025 | \[-0.5028, 1.0749\] |
| PFS | weibull | tmb_braf1 | 0.2279 | 0.2874 | \[-0.3355, 0.7913\] |
| PFS | weibull | RxExperimental arm:crp1 | 0.5269 | 0.4991 | \[-0.4512, 1.5050\] |
| PFS | weibull | RxExperimental arm:tmb_braf1 | 0.2938 | 0.4185 | \[-0.5265, 1.1142\] |
| PFS | weibullph | shape | 1.3913 | 0.1383 | \[1.1450, 1.6906\] |
| PFS | weibullph | scale | 0.0058 | 0.0066 | \[0.0006, 0.0531\] |
| PFS | weibullph | Age | -0.0056 | 0.0150 | \[-0.0350, 0.0238\] |
| PFS | weibullph | sex1 | 0.0293 | 0.2747 | \[-0.5090, 0.5677\] |
| PFS | weibullph | RxExperimental arm | 0.6002 | 0.3696 | \[-0.1241, 1.3246\] |
| PFS | weibullph | crp1 | -0.3980 | 0.5603 | \[-1.4962, 0.7003\] |
| PFS | weibullph | tmb_braf1 | -0.3171 | 0.4013 | \[-1.1037, 0.4695\] |
| PFS | weibullph | RxExperimental arm:crp1 | -0.7331 | 0.6940 | \[-2.0933, 0.6272\] |
| PFS | weibullph | RxExperimental arm:tmb_braf1 | -0.4088 | 0.5834 | \[-1.5522, 0.7346\] |
| PFS | llogis | shape | 2.0830 | 0.2182 | \[1.6964, 2.5578\] |
| PFS | llogis | scale | 49.9798 | 40.5110 | \[10.2060, 244.7566\] |
| PFS | llogis | Age | -0.0019 | 0.0118 | \[-0.0250, 0.0211\] |
| PFS | llogis | sex1 | -0.2055 | 0.2112 | \[-0.6195, 0.2084\] |
| PFS | llogis | RxExperimental arm | -0.7299 | 0.3328 | \[-1.3821, -0.0777\] |
| PFS | llogis | crp1 | 0.3212 | 0.4396 | \[-0.5404, 1.1827\] |
| PFS | llogis | tmb_braf1 | 0.0041 | 0.3311 | \[-0.6448, 0.6530\] |
| PFS | llogis | RxExperimental arm:crp1 | 0.9115 | 0.5562 | \[-0.1787, 2.0017\] |
| PFS | llogis | RxExperimental arm:tmb_braf1 | 0.4037 | 0.4481 | \[-0.4745, 1.2819\] |
| PFS | lognormal | meanlog | 3.8454 | 0.7703 | \[2.3357, 5.3552\] |
| PFS | lognormal | sdlog | 0.8258 | 0.0741 | \[0.6926, 0.9847\] |
| PFS | lognormal | Age | -0.0022 | 0.0109 | \[-0.0235, 0.0191\] |
| PFS | lognormal | sex1 | -0.1827 | 0.2054 | \[-0.5852, 0.2199\] |
| PFS | lognormal | RxExperimental arm | -0.5874 | 0.3175 | \[-1.2098, 0.0350\] |
| PFS | lognormal | crp1 | 0.1876 | 0.4336 | \[-0.6622, 1.0374\] |
| PFS | lognormal | tmb_braf1 | 0.1246 | 0.3204 | \[-0.5033, 0.7525\] |
| PFS | lognormal | RxExperimental arm:crp1 | 1.0230 | 0.5358 | \[-0.0270, 2.0731\] |
| PFS | lognormal | RxExperimental arm:tmb_braf1 | 0.2211 | 0.4332 | \[-0.6279, 1.0702\] |
| PFS | gamma | shape | 1.7902 | 0.2925 | \[1.2996, 2.4659\] |
| PFS | gamma | rate | 0.0466 | 0.0343 | \[0.0110, 0.1973\] |
| PFS | gamma | Age | -0.0038 | 0.0107 | \[-0.0248, 0.0171\] |
| PFS | gamma | sex1 | 0.0554 | 0.1987 | \[-0.3341, 0.4448\] |
| PFS | gamma | RxExperimental arm | 0.4732 | 0.2760 | \[-0.0676, 1.0141\] |
| PFS | gamma | crp1 | -0.2964 | 0.4125 | \[-1.1050, 0.5121\] |
| PFS | gamma | tmb_braf1 | -0.1892 | 0.2963 | \[-0.7698, 0.3915\] |
| PFS | gamma | RxExperimental arm:crp1 | -0.6076 | 0.5076 | \[-1.6026, 0.3873\] |
| PFS | gamma | RxExperimental arm:tmb_braf1 | -0.2977 | 0.4200 | \[-1.1209, 0.5254\] |
| PFS | gompertz | shape | 0.0094 | 0.0034 | \[0.0027, 0.0160\] |
| PFS | gompertz | rate | 0.0182 | 0.0183 | \[0.0025, 0.1304\] |
| PFS | gompertz | Age | -0.0035 | 0.0150 | \[-0.0328, 0.0258\] |
| PFS | gompertz | sex1 | 0.0247 | 0.2734 | \[-0.5111, 0.5605\] |
| PFS | gompertz | RxExperimental arm | 0.5788 | 0.3688 | \[-0.1440, 1.3015\] |
| PFS | gompertz | crp1 | -0.3481 | 0.5575 | \[-1.4407, 0.7446\] |
| PFS | gompertz | tmb_braf1 | -0.3710 | 0.4045 | \[-1.1639, 0.4218\] |
| PFS | gompertz | RxExperimental arm:crp1 | -0.7779 | 0.7000 | \[-2.1499, 0.5940\] |
| PFS | gompertz | RxExperimental arm:tmb_braf1 | -0.4210 | 0.5876 | \[-1.5727, 0.7306\] |
| PFS | gengamma | mu | 3.8407 | 0.8179 | \[2.2376, 5.4438\] |
| PFS | gengamma | sigma | 0.8256 | 0.0757 | \[0.6898, 0.9882\] |
| PFS | gengamma | Q | 0.0087 | 0.5810 | \[-1.1301, 1.1475\] |
| PFS | gengamma | Age | -0.0021 | 0.0127 | \[-0.0270, 0.0228\] |
| PFS | gengamma | sex1 | -0.1809 | 0.2353 | \[-0.6421, 0.2802\] |
| PFS | gengamma | RxExperimental arm | -0.5862 | 0.3257 | \[-1.2247, 0.0522\] |
| PFS | gengamma | crp1 | 0.1900 | 0.4586 | \[-0.7087, 1.0888\] |
| PFS | gengamma | tmb_braf1 | 0.1247 | 0.3204 | \[-0.5031, 0.7526\] |
| PFS | gengamma | RxExperimental arm:crp1 | 1.0172 | 0.6584 | \[-0.2733, 2.3076\] |
| PFS | gengamma | RxExperimental arm:tmb_braf1 | 0.2228 | 0.4497 | \[-0.6585, 1.1042\] |
| PFS | genf | mu | 3.8415 | 0.8179 | \[2.2384, 5.4446\] |
| PFS | genf | sigma | 0.8254 | 0.0757 | \[0.6896, 0.9879\] |
| PFS | genf | Q | 0.0079 | 0.5817 | \[-1.1322, 1.1480\] |
| PFS | genf | P | 0.0004 | 0.0143 | \[0.0000, 3767584399327780893062886026626.0000\] |
| PFS | genf | Age | -0.0021 | 0.0126 | \[-0.0269, 0.0226\] |
| PFS | genf | sex1 | -0.1810 | 0.2354 | \[-0.6423, 0.2803\] |
| PFS | genf | RxExperimental arm | -0.5864 | 0.3257 | \[-1.2249, 0.0520\] |
| PFS | genf | crp1 | 0.1901 | 0.4586 | \[-0.7087, 1.0889\] |
| PFS | genf | tmb_braf1 | 0.1243 | 0.3203 | \[-0.5034, 0.7521\] |
| PFS | genf | RxExperimental arm:crp1 | 1.0173 | 0.6587 | \[-0.2736, 2.3083\] |
| PFS | genf | RxExperimental arm:tmb_braf1 | 0.2235 | 0.4497 | \[-0.6579, 1.1049\] |

Coefficient estimates for all fitted parametric models

# Ordering-Constrained Selection Summary

| Endpoint                  | Distribution |    AIC |
|:--------------------------|:-------------|-------:|
| Overall survival          | gamma        | 667.31 |
| Progression-free survival | gamma        | 636.86 |

Selected minimum-combined-AIC pair satisfying OS \>= PFS

The pair is selected jointly. Consequently, an endpoint’s selected
family need not be its independently lowest-AIC family. The full
candidate-pair audit and the combined-AIC penalty relative to
independent endpoint selection are stored in `models$ordered_selection`
and `models$best_fit$ordering_aic_penalty`, respectively.

------------------------------------------------------------------------

**Report completed on:** 2026-09-22  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 5.0
