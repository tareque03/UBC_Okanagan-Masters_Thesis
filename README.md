# Perturbation Theory of Linear Systems with Applications to Structured Epidemic Models

**Md Tareque Hossain**
MSc Thesis in Mathematics, The University of British Columbia (Okanagan), December 2023
Supervisor: Dr. Eric Foxall

📄 **[Read the full thesis (PDF)](thesis.pdf)**

## Abstract

We explore the impact of parameter changes on infection spread within structured populations. We study perturbations of linear models in both discrete and continuous time. By interpreting the resulting formulas in the context of an associated Markov chain, new insights are gained. Specifically, we explain how the rate of change of key quantities, that we represent with the so-called canonical flow, is determined by expected hitting times of corresponding states in the Markov chain, from different initial distributions. We give intuitive direct proofs of known results and obtain novel formulas in the case of a general primitive or irreducible Metzler matrix. We apply our results to two models from the literature: the fast pathogen movement limit of a structured model of Cholera spread, and the slow-fast decomposition of critical structured SIR models near the disease-free set. We explain how the latter context can be used to model the effects of control measures with uncontrolled response variables.

## Thesis outline

| Chapter | Topic |
|---|---|
| 1 | Introduction |
| 2 | Literature review |
| 3 | Linear models in discrete and continuous time (Perron-Frobenius theory, canonical flow of Perron vectors, group inverse and hitting-time formula) |
| 4 | Slow-fast decomposition of a structured SIR model (numerical algorithm and simulations) |
| 5 | Conclusion |

## Figures from the thesis

### Chapter 3: Linear models in discrete and continuous time

**Figure 3.1.** Transition diagram of the n-state life-cycle model.

![Figure 3.1](images/fig3_1_life_cycle_transition_diagram.png)

**Figure 3.2.** The path taken from leaf nodes ℓ and k to leaf node i, through the hub node 1 in a star network.

![Figure 3.2](images/fig3_2_star_network_path.png)

### Chapter 4: Slow-fast decomposition

**Figure 4.1.** The x variable changes quickly to approach the slow manifold and changes slowly along the slow manifold as λ is varied.

![Figure 4.1](images/fig4_1_slow_manifold.png)

**Figure 4.2.** Structured SIR model with small I and a restriction on B (taking η₂ = 0.2) such that α(Q) = 0.

![Figure 4.2](images/fig4_2_structured_sir_small_I.png)

**Figure 4.3.** Structured SIR model with different η₂ values on the given time intervals: 0–15: η₂ = 0.1, 15–30: η₂ = 0.25, 30–45: η₂ = 0.1, 45–60: η₂ = 0.25.

![Figure 4.3](images/fig4_3_structured_sir_varying_eta2.png)

## Code

The simulations in Chapter 4 (Section 4.1.2) are in the `code/` folder. In the code, `rho2` corresponds to η₂ in the thesis.

| Script | Reproduces |
|---|---|
| [`code/fig4_2_fixed_rho2.py`](code/fig4_2_fixed_rho2.py) | Figure 4.2, with ρ₂ fixed at 0.2 for 60 days |
| [`code/fig4_3_varying_rho2.py`](code/fig4_3_varying_rho2.py) | Figure 4.3, with ρ₂ switched every 15 days (0.1, 0.25, 0.1, 0.25) |

Both scripts use the same parameters: ε = 0.01, γ = (0.4, 0.2), initial susceptible proportions (0.6, 0.4), time step h = 0.05, and T = 60.

### Run

```bash
pip install -r requirements.txt
cd code
python fig4_2_fixed_rho2.py
python fig4_3_varying_rho2.py
```

Each script displays the plot and saves a PNG (`fig4_2_output.png`, `fig4_3_output.png`) in the current folder.

## Citation

```
Hossain, M. T. (2023). Perturbation theory of linear systems with applications to
structured epidemic models. MSc thesis, The University of British Columbia (Okanagan).
```

© Md Tareque Hossain, 2023
