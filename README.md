# Three-Body Recombination Model

This repository contains the Python code used for calculations and plotting graphs presented in the paper on three-body recombination of y-charged particles in the early Universe. The code models the evolution of the relative density of these particles as a function of the temperature of photons from the cosmic microwave background (CMB).

## Overview

The code performs the following tasks:

1. **Computation of Relative Density**: It calculates the evolution of the relative density of y-charged particles for different recombination mechanisms (three-body, radiative, and classical recombination).
   
2. **Limitations of Three-Body Recombination**: The code also evaluates the applicability of the three-body recombination formula within certain model parameters and presents a restriction graph showing the boundaries of its validity.

3. **Visualization**: The script generates multiple plots that compare different recombination mechanisms and highlight the constraints on the use of the three-body recombination formula.

## Main Features

- **Temperature-Dependent Calculations**: The code computes the relative density of y-charged particles as a function of the temperature of photons.
- **Multiple Recombination Models**: It includes calculations for three types of recombination:
  - Three-body recombination.
  - Radiative recombination.
  - Classical recombination.
- **Validity Range Calculation**: The code calculates the limits of applicability for the three-body recombination formula and plots these regions on the graphs.
- **Graphical Output**: The resulting graphs are saved and can be used for further analysis or directly in publications.

## Requirements

To run the code, you need the following dependencies:
- Python 3.x
- `numpy` for numerical calculations.
- `matplotlib` for plotting.
- `scipy` for scientific computing functions.

To install the required dependencies, you can run:
```bash
pip install numpy matplotlib scipy
```

## Usage

1. **Run the Script**: Simply execute the Python script `rec.py` to perform all calculations and generate the graphs.
   ```bash
   python rec.py
   ```

2. **Adjust Parameters**: You can modify the model parameters (e.g., particle masses, temperature range) directly in the script by changing the respective variables at the top of the file.

3. **Output**: The script will output the following plots:
   - **Relative density of y-charged particles** vs. **Temperature** for various recombination mechanisms.
   - **Limitations of three-body recombination** with shaded regions indicating the valid parameter space.

   The graphs will be saved as `.png` files in the working directory.

## Key Sections of the Code

- **Computation of Relative Density**: This section solves the equations that describe the change in relative density ($r/r_0$) based on the recombination mechanism.
- **Limitations of Applicability**: This part of the code calculates the boundaries where the three-body recombination formula is valid.
- **Plotting Functions**: The final part of the code generates the graphs, formatting the curves and adding labels and legends.

## Graphs

- The script produces the following graphs:
  - **Three-Body Recombination vs. Other Mechanisms**: A comparison of different recombination processes.
  - **Applicability of the Three-Body Formula**: Displays the region where the formula is valid, with an overlay on the recombination curves.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

The code and models are based on the methods and formulas presented in the research paper: IN PROCESS.

---
