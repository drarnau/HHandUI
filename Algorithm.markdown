# Algorithm for _On Households and Unemployment Insurance_ by Sekyu Choi and Arnau Valladares-Esteban
<!-- ### To-do/keep in mind
- Codes are for computers to execute and for humans to read.
- The code for couples should:
  - Have separate utility for the male and the female within the same maximization problem -->

### Functions
##### Need to be created
- Instantaneous utility
  - Inputs: parameters, consumption (real), employed [0,1], search [0,1], search cost (element of a vector/array)
  - Output: utility
  - Given log utility for consumption, if consumption is negative, the function should return a very negative number instead of evaluating
- UI Benefits
  - Input: productivity shock, average productivity (equilibrium object), wage (equilibrium object)
  - Output: UI benefit

##### Already out there
- Linear grid
- Curved grid
- Tauchen discretization for AR(1)
- Rowenhorst discretization for AR(1)
- Linear interpolation
  - *KMRS* do not allow extrapolation
- Golden-section search. Depending on implementation:
  - Can we use the already-out-there codes?
  - Is it efficient to use the already-out-there codes in terms of computation and design time vs. creating our application-specific implementation?

### Steps
1. Load/Input parameters
  - Keep in mind difference between *model* parameters vs. *numerical* parameters (e.g. grid points).
2. Create grids and transition matrices for:
  - Assets (value function iteration): *KMRS* use 48 grid points with more points close to 0. Range: [0,1440].
  - Assets (simulation): *KMRS* use 1000 grid points evenly distributed. Same range as before.
  - Productivity shock: *KMRS* do it through Tauchen. 20 grid points that cover 2 standard deviations.
  - Match quality shock: *KMRS* do it through Tauchen. 7 grid points that cover 2 standard deviations.
3. Define functions for the expected values of next period of the Bellman equations
  - The expected values have different structure across value functions. Looks better to create a different function associated to each value function.
  - The expected values are functions of assets, productivity, and UI eligibility. The other shocks are unconditional on state variables.
  - *KMRS* use linear interpolation to have an approximation of the expected value tomorrow given a choice for assets tomorrow. They treat the whole expectation operator as one function and then they interpolate.
  - An alternative is to respect that each value function *within* the expectation operator can be interpolated and *after* compute the expected value. The latter seems more accurate and not particularly costly. The problem is that the *max* operators might create *jumps* which imply a *violation* of the assumptions needed for the golden search algorithm.
4. Define functions for each value function as a function of assets tomorrow
  - Each function needs to be properly link to the expected value functions.
  - Assets tomorrow affect:
    - Instantaneous utility today through the budget constraint (consumption).
    - Expected value tomorrow.
5. Initialize value functions
  - There are 3 value functions: W, U, and N. J and W are straight forward to compute, once we know the values of W, U, and N.
  - Start with a guess
6. Use Golden search to find policy functions
  - Update guess
  - Keep iterating until guess and outcome of golden search are very similar
