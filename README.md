# Hematopoietic Stem Cell Dynamics – PNAS Modeling Code

This repository contains MATLAB code, data, and sample figures from the project on **Dynamically adjusted cell fate decisions and resilience to mutant invasion during steady-state hematopoiesis revealed by an experimentally parameterized mathematical model**, as published in *Proceedings of the National Academy of Sciences (PNAS)*.  

The code implements mathematical models of HSC population dynamics, equilibria, and the effects of mutant cell introduction. It is intended for reproducibility, reference, and further exploration by researchers interested in computational biology, mathematical modeling, and systems biology.  

---

## Repository Structure

```
hematopoietic-modeling/
├─ src/
│  ├─ MendozaCode_MathExPLR2022.m    # Old exploratory script (math refresher)
│  ├─ SystemModelTrial1.m            # First trial at system modeling
│  ├─ Trial2.m                       # Second trial script
│  ├─ EquilibriaAnalytics.m          # Analysis of equilibrium states
│  ├─ MutantCode1.m                  # Mutant population simulation
│  ├─ MutantCode2.m                  # Extended mutant simulation
│  ├─ m1.mat                         # Example dataset
│  ├─ m2.mat                         # Example dataset
│  ├─ m3.mat                         # Example dataset
│  ├─ m4.mat                         # Example dataset
│  ├─ m5.mat                         # Example dataset
├─ figures/
│  ├─ AllMutant[S;0.05].jpg          # Sample simulation output
│  ├─ AllMutant[S;0.01].jpg          # Sample simulation output
│  ├─ AllMutants[W;0.01].jpg         # Sample simulation output
│  ├─ AllMutants[W;0.05].jpg         # Sample simulation output
│  ├─ FinalMutant[S;0.01].jpg        # Final version of simulation
│  ├─ FinalMutant[S;0.05].jpg        # Final version of simulation
│  ├─ FinalMutant[W;0.01].jpg        # Final version of simulation
│  ├─ FinalMutant[W;0.05].jpg        # Final version of simulation
└─ README.md
```

## Requirements

- **MATLAB R2022b or later**  
- Toolboxes:  
  - Symbolic Math Toolbox  
  - Statistics and Machine Learning Toolbox (recommended)  

---

## How to Run

1. Clone or download this repository.  
2. Open MATLAB and navigate to the `src/` folder.  
3. Run the desired script, for example:  
