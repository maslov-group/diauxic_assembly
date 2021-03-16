## Complementary resource preferences spontaneously emerge in diauxic microbial communities

This repository contains the code used to simulate the community assembly of diauxic microbes under serial dilution, under conditions that are considered in this paper.  

#### Parameters of the main simulation

| Parameter                                                    | Value | Notes                                                        |
| ------------------------------------------------------------ | ----- | ------------------------------------------------------------ |
| Number of species                                            | 11520 | Equals to #(possible preference orders) * #(species of each preference order). |
| Number of resources                                          | 4     | #(possible preference orders) = 4!.                          |
| Number of species that shares the same preference order      | 480   |                                                              |
| Mean of growth rates' distribution before rectified          | 0.25  | Growth rate unit is hr^-1^.                                  |
| Standard error of growth rates' distribution before rectified | 0.05  |                                                              |
| Threshold of rectification                                   | 0.1   |                                                              |
| Yield of all species                                         | 0.5   |                                                              |
| Lag time of all species on all resources                     | 0     | Time unit is hr. Lags of 1 hr and 5 hr are also considered.  |
| Dilution factor                                              | 100   | In supp. a dilution factor of 2 is also considered.          |
| Period of dilution                                           | 24    |                                                              |
| Concentration of every resource at the beginning of the first cycle | 1e0   | All concentration unit is a.u. In main text all resources are provided at equal amount; in supp. a ratio of 100:1:1:1 is considered with the total resources concentration being constant. |
| Microbe concentration at invasion                            | 1e-6  |                                                              |
| Threshold under which a species is seen as eliminated        | 1e-6  |                                                              |
| Threshold under which a resource is seen as depleted         | 1e-9  |                                                              |

ModulesWithLags.py contains the functions used in the simulation. 

Parameters.py contains the value of parameters. 

RunsWithLags.py and Example_run.ipynb run the simulation. 