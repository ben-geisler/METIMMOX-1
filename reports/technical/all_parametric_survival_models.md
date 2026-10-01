# All Parametric Survival Models
Ben Geisler
2026-10-01

- [Overview](#overview)
- [Information Criteria](#information-criteria)
- [Coefficient Tables](#coefficient-tables)
- [Ordering-Constrained Selection
  Summary](#ordering-constrained-selection-summary)
  - [Validation warnings](#validation-warnings)

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
|     gengamma | 496.05 | 518.24 | Progression-free survival |      0.00 |      0.23 |
|    lognormal | 498.04 | 518.02 | Progression-free survival |      1.99 |      0.00 |
|         genf | 498.07 | 522.48 | Progression-free survival |      2.02 |      4.46 |
|        gamma | 500.16 | 520.13 | Progression-free survival |      4.11 |      2.11 |
|      weibull | 501.29 | 521.26 | Progression-free survival |      5.24 |      3.24 |
|    weibullph | 501.29 | 521.26 | Progression-free survival |      5.24 |      3.24 |
|       llogis | 501.56 | 521.53 | Progression-free survival |      5.51 |      3.52 |
|     gompertz | 505.31 | 525.29 | Progression-free survival |      9.26 |      7.27 |
|  exponential | 509.48 | 527.24 | Progression-free survival |     13.43 |      9.22 |

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
| PFS | exponential | rate | 0.0330 | 0.0359 | \[0.0039, 0.2789\] |
| PFS | exponential | Age | -0.0073 | 0.0158 | \[-0.0383, 0.0237\] |
| PFS | exponential | sex1 | 0.1380 | 0.2973 | \[-0.4447, 0.7208\] |
| PFS | exponential | RxExperimental arm | 0.5189 | 0.4102 | \[-0.2850, 1.3228\] |
| PFS | exponential | crp1 | -0.0651 | 0.6570 | \[-1.3529, 1.2227\] |
| PFS | exponential | tmb_braf1 | -0.3421 | 0.4728 | \[-1.2689, 0.5846\] |
| PFS | exponential | RxExperimental arm:crp1 | -0.8879 | 0.7900 | \[-2.4364, 0.6605\] |
| PFS | exponential | RxExperimental arm:tmb_braf1 | -0.3330 | 0.6409 | \[-1.5892, 0.9232\] |
| PFS | weibull | shape | 1.4630 | 0.1602 | \[1.1804, 1.8132\] |
| PFS | weibull | scale | 31.0771 | 23.6241 | \[7.0045, 137.8805\] |
| PFS | weibull | Age | 0.0072 | 0.0110 | \[-0.0143, 0.0288\] |
| PFS | weibull | sex1 | -0.0851 | 0.2071 | \[-0.4910, 0.3208\] |
| PFS | weibull | RxExperimental arm | -0.4777 | 0.2815 | \[-1.0294, 0.0740\] |
| PFS | weibull | crp1 | -0.0283 | 0.4553 | \[-0.9207, 0.8641\] |
| PFS | weibull | tmb_braf1 | 0.2816 | 0.3312 | \[-0.3676, 0.9308\] |
| PFS | weibull | RxExperimental arm:crp1 | 0.8651 | 0.5516 | \[-0.2160, 1.9463\] |
| PFS | weibull | RxExperimental arm:tmb_braf1 | 0.3664 | 0.4548 | \[-0.5250, 1.2578\] |
| PFS | weibullph | shape | 1.4630 | 0.1602 | \[1.1804, 1.8132\] |
| PFS | weibullph | scale | 0.0066 | 0.0082 | \[0.0006, 0.0758\] |
| PFS | weibullph | Age | -0.0106 | 0.0161 | \[-0.0422, 0.0210\] |
| PFS | weibullph | sex1 | 0.1245 | 0.3025 | \[-0.4683, 0.7173\] |
| PFS | weibullph | RxExperimental arm | 0.6989 | 0.4160 | \[-0.1164, 1.5141\] |
| PFS | weibullph | crp1 | 0.0415 | 0.6663 | \[-1.2645, 1.3475\] |
| PFS | weibullph | tmb_braf1 | -0.4120 | 0.4852 | \[-1.3629, 0.5390\] |
| PFS | weibullph | RxExperimental arm:crp1 | -1.2657 | 0.8155 | \[-2.8641, 0.3328\] |
| PFS | weibullph | RxExperimental arm:tmb_braf1 | -0.5361 | 0.6688 | \[-1.8468, 0.7747\] |
| PFS | llogis | shape | 2.0380 | 0.2305 | \[1.6329, 2.5437\] |
| PFS | llogis | scale | 33.7290 | 29.6553 | \[6.0202, 188.9727\] |
| PFS | llogis | Age | 0.0021 | 0.0127 | \[-0.0228, 0.0270\] |
| PFS | llogis | sex1 | -0.2072 | 0.2361 | \[-0.6699, 0.2555\] |
| PFS | llogis | RxExperimental arm | -0.5991 | 0.3538 | \[-1.2925, 0.0944\] |
| PFS | llogis | crp1 | 0.0566 | 0.5385 | \[-0.9989, 1.1121\] |
| PFS | llogis | tmb_braf1 | 0.1411 | 0.3842 | \[-0.6119, 0.8941\] |
| PFS | llogis | RxExperimental arm:crp1 | 1.1401 | 0.6431 | \[-0.1204, 2.4007\] |
| PFS | llogis | RxExperimental arm:tmb_braf1 | 0.2831 | 0.5028 | \[-0.7023, 1.2685\] |
| PFS | lognormal | meanlog | 3.5067 | 0.8264 | \[1.8869, 5.1264\] |
| PFS | lognormal | sdlog | 0.8171 | 0.0817 | \[0.6717, 0.9939\] |
| PFS | lognormal | Age | 0.0014 | 0.0115 | \[-0.0213, 0.0240\] |
| PFS | lognormal | sex1 | -0.1929 | 0.2203 | \[-0.6247, 0.2390\] |
| PFS | lognormal | RxExperimental arm | -0.5186 | 0.3322 | \[-1.1697, 0.1325\] |
| PFS | lognormal | crp1 | -0.0861 | 0.4850 | \[-1.0367, 0.8645\] |
| PFS | lognormal | tmb_braf1 | 0.2160 | 0.3520 | \[-0.4738, 0.9058\] |
| PFS | lognormal | RxExperimental arm:crp1 | 1.2897 | 0.5813 | \[0.1503, 2.4291\] |
| PFS | lognormal | RxExperimental arm:tmb_braf1 | 0.1801 | 0.4617 | \[-0.7249, 1.0851\] |
| PFS | gamma | shape | 1.8980 | 0.3348 | \[1.3432, 2.6818\] |
| PFS | gamma | rate | 0.0644 | 0.0511 | \[0.0136, 0.3045\] |
| PFS | gamma | Age | -0.0069 | 0.0112 | \[-0.0288, 0.0150\] |
| PFS | gamma | sex1 | 0.1142 | 0.2114 | \[-0.3002, 0.5287\] |
| PFS | gamma | RxExperimental arm | 0.4970 | 0.2963 | \[-0.0837, 1.0777\] |
| PFS | gamma | crp1 | 0.0190 | 0.4703 | \[-0.9028, 0.9408\] |
| PFS | gamma | tmb_braf1 | -0.2555 | 0.3378 | \[-0.9175, 0.4066\] |
| PFS | gamma | RxExperimental arm:crp1 | -0.9574 | 0.5642 | \[-2.0633, 0.1485\] |
| PFS | gamma | RxExperimental arm:tmb_braf1 | -0.3255 | 0.4553 | \[-1.2178, 0.5669\] |
| PFS | gompertz | shape | 0.0115 | 0.0045 | \[0.0027, 0.0204\] |
| PFS | gompertz | rate | 0.0240 | 0.0269 | \[0.0027, 0.2163\] |
| PFS | gompertz | Age | -0.0085 | 0.0161 | \[-0.0399, 0.0230\] |
| PFS | gompertz | sex1 | 0.1544 | 0.3009 | \[-0.4354, 0.7442\] |
| PFS | gompertz | RxExperimental arm | 0.6906 | 0.4189 | \[-0.1304, 1.5116\] |
| PFS | gompertz | crp1 | 0.0460 | 0.6650 | \[-1.2574, 1.3494\] |
| PFS | gompertz | tmb_braf1 | -0.3986 | 0.4795 | \[-1.3383, 0.5412\] |
| PFS | gompertz | RxExperimental arm:crp1 | -1.2884 | 0.8330 | \[-2.9211, 0.3443\] |
| PFS | gompertz | RxExperimental arm:tmb_braf1 | -0.6333 | 0.6826 | \[-1.9712, 0.7047\] |
| PFS | gengamma | mu | 4.0298 | 0.3623 | \[3.3196, 4.7399\] |
| PFS | gengamma | sigma | 0.3510 | 0.0468 | \[0.2703, 0.4558\] |
| PFS | gengamma | Q | -3.6054 | 0.1448 | \[-3.8892, -3.3216\] |
| PFS | gengamma | Age | -0.0213 | 0.0041 | \[-0.0293, -0.0133\] |
| PFS | gengamma | sex1 | -0.3585 | 0.1837 | \[-0.7186, 0.0016\] |
| PFS | gengamma | RxExperimental arm | -0.2495 | 0.1824 | \[-0.6070, 0.1080\] |
| PFS | gengamma | crp1 | -0.2335 | 0.2586 | \[-0.7403, 0.2733\] |
| PFS | gengamma | tmb_braf1 | 0.2753 | 0.2476 | \[-0.2100, 0.7605\] |
| PFS | gengamma | RxExperimental arm:crp1 | 1.6668 | 0.2419 | \[1.1926, 2.1411\] |
| PFS | gengamma | RxExperimental arm:tmb_braf1 | -0.5251 | 0.2970 | \[-1.1073, 0.0571\] |
| PFS | genf | mu | 4.0358 | 0.3487 | \[3.3524, 4.7193\] |
| PFS | genf | sigma | 0.3368 | 0.0444 | \[0.2601, 0.4360\] |
| PFS | genf | Q | -3.6666 | 0.1637 | \[-3.9875, -3.3458\] |
| PFS | genf | P | 0.0000 | 0.0004 | \[0.0000, 28112.4609\] |
| PFS | genf | Age | -0.0216 | 0.0038 | \[-0.0291, -0.0141\] |
| PFS | genf | sex1 | -0.3871 | 0.1890 | \[-0.7576, -0.0167\] |
| PFS | genf | RxExperimental arm | -0.2332 | 0.1776 | \[-0.5814, 0.1150\] |
| PFS | genf | crp1 | -0.2121 | 0.2646 | \[-0.7307, 0.3065\] |
| PFS | genf | tmb_braf1 | 0.3062 | 0.2475 | \[-0.1788, 0.7913\] |
| PFS | genf | RxExperimental arm:crp1 | 1.6559 | 0.2417 | \[1.1822, 2.1295\] |
| PFS | genf | RxExperimental arm:tmb_braf1 | -0.5689 | 0.2967 | \[-1.1503, 0.0126\] |

Coefficient estimates for all fitted parametric models

# Ordering-Constrained Selection Summary

| Endpoint                  | Distribution |    AIC |
|:--------------------------|:-------------|-------:|
| Overall survival          | gamma        | 667.31 |
| Progression-free survival | gamma        | 500.16 |

Selected minimum-combined-AIC pair satisfying OS \>= PFS

The pair is selected jointly. Consequently, an endpoint’s selected
family need not be its independently lowest-AIC family. The full
candidate-pair audit and the combined-AIC penalty relative to
independent endpoint selection are stored in `models$ordered_selection`
and `models$best_fit$ordering_aic_penalty`, respectively.

------------------------------------------------------------------------

**Report completed on:** 2026-10-01  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 5.0

## Validation warnings

No warnings recorded during rendering.
