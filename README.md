## Complementary resource preferences spontaneously emerge in diauxic microbial communities

Many microbes grow diauxically, utilizing the available resources one at a time rather than simultaneously. The properties of communities of microbes growing diauxically remain poorly understood, largely due to a lack of theory and models of such communities.  Here, we develop and study a minimal model of diauxic microbial communities assembling in a serially diluted culture. We find that unlike co-utilizing communities, diauxic community assembly repeatably and spontaneously leads to communities with complementary resource preferences, i.e., communities where species prefer different resources as their top choice. Simulations and theory explain that the emergence of complementarity is driven by the disproportionate contribution of the top choice resource to the growth of a diauxic species. Additionally, we develop a geometric approach for analyzing serially diluted communities, with or without diauxie, which intuitively explains several additional emergent community properties, such as the apparent lack of species which grow fastest on a resource other than their most preferred resource. Overall, our work provides testable predictions for the assembly of natural as well as synthetic communities of diauxically shifting microbes.

#### Description of code:

2 parts of the code are present. 

- In folder community_assembly there is code simulating the community assembly of diauxic species under serial dilutions, which corresponds to our main result in figure xxx. 
  - ModulesWithLags.py contains the functions used in the simulation. 
  - Parameters.py contains the value of parameters (for details see the table below).
  - RunsWithLags.py runs the simulation and generates data. The data's size we used for our main text results has exceeded 200MB, so we here have a small subset separated from the main data. This example data can be generated by running the simulation with certain parameters. 
  - The jupyter notebook processes the data and make plots in fig. 2 and 3.
- In folder pairwise_invasions there is code simulating pairwise invasion experiments of laboratory E. coli species, which corresponds to our newly added figure 6. 
  - The barthe_strains.xlsx summarizes data from Barthe et al. (see reference in main text) of the laboratory strains we chose to simulate. 
  - The jupyter notebook does simulations and displays the outcomes, as in fig. 6. 

----

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
| Lag time of all species on all resources                     | 0     | Time unit is hr. Can me changed to any other non-zero number in simulations. |
| Dilution factor                                              | 100   | In supp. a dilution factor of 2 is also considered.          |
| Period of dilution                                           | 24    | In unit of hours                                             |
| Concentration of every resource at the beginning of the first cycle | 1e0   | All concentration unit is a.u. In main text all resources are provided at equal amount; in supp. a ratio of 100:1:1:1 is considered with the total resources concentration being constant. |
| Microbe concentration at invasion                            | 1e-6  |                                                              |
| Threshold under which a species is seen as eliminated        | 1e-6  |                                                              |
| Threshold under which a resource is seen as depleted         | 1e-9  |                                                              |

