# Diauxie

### Without considering the lags:

The first idea is inspired by the C and N model. Generate a species pool with random growth rates, yields and preference lists, and find uninvadable communities under serial dilutions. As we find these states, we may try to find there corresponding volume in the nutirent supply space, and see if there is multistability. 

However, a major problem is that the rules (such as the 2 basic rules in C and N model) are not clear in this setup. What we only know is that one species survives if for one nutrient it prevails. Without effcient reduction of states, it is hard to find all the uninvadable states that are possible under at least one set of nutrient supply. 

Another problem is that the current algorithm can only find one uninvadable state per pool. Finding more requires significantly more time.

--------

The second idea is of game theory. The first step is to only have pure strategies of using preference list of nutrient. For each speices there is a set of fixed growth rates, and the preference order can be changed to reach its best payoff. We try to find the Nash Equilibrium of the system. Some trivial NEs include the states where certain species never survive no matter what strategy it uses. This is also the most usual case. It is also possible to check the volume in growth rate space each NE occupies in this case. 

The second step is to have mixed strategies, in which we might be able to map these mixed strategies to certain coutilization proportions.

Another thing to do is this: discretize the growth rates and set up a budget. eg: choose between 1, 2, 3 and set budget 6. there are 7 strategies: permutations of 1, 2, 3 and 2, 2, 2. 

It is not yet clear which evolves faster: the pref order or the growth rates. We may assume that they differ by at least an order so that we could do it one after one.


