## Microbial communities of sequential consumers self-organizeinto complementary resource preference (Draft title, to be decided later)

This repository contains the code used to simulate the community assembly of diauxic microbes  under serial dilution. 

#### modules.ipynb
The library that has functions used to simulate the community assembly process (Without considering the lag of diauxic behaviors)

#### modulesWithLags.ipynb
The patch that contains functions that adds the feature of lags. 

#### (Other notebooks that shows how to make figures)

| Parameter  | Meaning                                                      | Value | Notes                                            |
| ---------- | ------------------------------------------------------------ | ----- | ------------------------------------------------ |
| $S$        | Number of species                                            | 11520 | $M!\cdot S_0$                                    |
| $M$        | Number of resources                                          | 4     |                                                  |
| $S_0$      | Number of species that shares the same preference order      | 480   |                                                  |
| $\mu_g$    | Mean of growth rates' distribution before rectified          | 6     |                                                  |
| $\sigma_g$ | Standard error of growth rates' distribution before rectified | 1     |                                                  |
| $\theta_g$ | Threshold of rectification                                   | 0.1   |                                                  |
| $Y$        | Yield of all species                                         | 0.5   |                                                  |
| $\tau$     | Lag time of all species on all resources                     | 0     | 0.1 and 1 in supp. materials                     |
| $D$        | Dilution factor                                              | 100   | 20 and 5 in supp.materials **NEEDS MODIF LATER** |
| $T$        | Period of dilution                                           | 10    |                                                  |
| $c_{0}$    | Concentration of every resource at the beginning of the first cycle | 1e0   |                                                  |
| $b_0$      | Microbe concentration at invasion                            | 1e-6  |                                                  |
| $b_{thrs}$ | Threshold under which a species is seen as eliminated        | 1e-6  |                                                  |
| $c_{thrs}$ | Threshold under which a resource is seen as depleted         | 1e-9  |                                                  |
| $N_{sys}$  | Number of simulated independent systems                      | 1600  | **NEEDS MODIF. LATER**                           |

