# vitraPro: simulation of viral transduction and propagation for biomanufacturing

Simulation of viral transduction and propagation for up to two viral species in suspension cell cultures with model and numerics presented in:

Destro, F. and Braatz, R. D. (2024). Efficient simulation of viral transduction and propagation for biomanufacturing, biorXiv. DOI: 10.1101/2024.03.30.587435

MATLAB and Python implementations are available. A (faster) model for systems with only one viral species is also implemented. 

## Repository structure
- `simulator`: folder containing code for running simulations
  - `matlab`: matlab implementation of the simulator
    - `one virus`: simulate systems with one viral species
    - `two viruses`: simulate systems with two viral species
  - `python`: python implementation of the simulator
    - `two viruses`: simulate systems with two viral species
- `case studies`: MATLAB code used to generate the results reported in the paper
  - `CS_1_S1`: case studies 1 and S1 (simulate by running `Run_CS1.m` or `Run_CS_S1.m`)
  - `CS_2`: case study 2 (simulate by running `Run_CS2.m`)
  - `CS_3`: case study 3 (simulate by running `Run_CS3.m`)

## Simulation
For starting a simulation, go to a simulator folder based on the desired language (MATLAB/Python) and system type (one/two viral species).
When only one virus is present in the system, the simulators for one and two viral species gave the same results, but the simulator for one viral species is faster.

Files in simulator folder `simulator/matlab/two viruses`:
- `RunSimulation.m`: external wrapper to be used for setting the simulation inputs and starting a simulation
- `def_parameters.m`: function to be used for setting the desired model parameters 
- `main.m`: function called by external wrapper to solve the model equations
- `inf_model.m`: function containing model equations
- `sample_figures.m`: script to generate sample plots

To start a simulation, specify the model parameters in `def_parameters.m` and the simulation inputs in `RunSimulation.m`. Then, execute the script `RunSimulation.m`.

The MATLAB simulator for systems with one virus and the Python simulator contain analogous files.

## Simulation input/output structure
`RunSimulation` contains detailed comments on simulation inputs and outputs and on how to specify the simulation inputs.

**Simulation inputs**
- Simulation duration
- Viable cell density at process onset
- Nonviable cell density at process onset
- Substrate concentration
- Presence of DIPs (Y/N)
- Viral inoculation conditions (two options):
  - Inoculation with virions: specify multiplicity of infection (inoculated plaque forming units per viable cell)
  - Inoculation with infected cells: specify concentration of inoculated infected cells. In this case, viral infection is initiated by the progeny released by the inoculated infected cells
- System type (batch, perfusion, or continuous)
  - Specified by appropriately setting dilution rate `D` and bleeding ratio `r`
  - For continuous and perfusion: specify uninfected cell and substrate concentration in feed
- Numerics
  - Infection age mesh: bin size, maximum infection age, infection age beyond which a cell can't get reinfected (if any)
  - Scheme for integrating ODEs after PDE discretization: RK23 or RK45
  - Expert users can also modify the parameter used for normalizing the states (`scaling`) and the RK convergence tolerance
- Generate sample figures (Y/N)

**Simulation outputs (systems with two viral species)**

      Scalar states: the output is a vector that represents a time profile
          T:              Uninfected cell concentration time profile - [cell/mL]        
          V1:             Virion 1 concentration time profile - [PFU/mL]            
          V2:             Virion 2 concentration time profile - [PFU/mL]   
          NV:             Nonviable cell concentration time profile  - [cell/mL]
          S:              Substrate concentration time profile - [nmol/mL]
          N1_pc_avg:      Time profile of avg concentration of virus 1 genome in 
                          infected cells - [vg/cell] (per cell basis)
          N2_pc_avg:      Time profile of avg concentration of virus 2 genome in 
                          infected cells - [vg/cell] (per cell basis)
          N1_Co_pc_avg:   Time profile of avg concentration of virus 1 genome in 
                          coinfected cells - [vg/cell] (per cell basis)
          N2_Co_pc_avg:   Time profile of avg concentration of virus 2 genome in 
                          coinfected cells - [vg/cell] (per cell basis)

      States distributed with respect to one infection age: the output is a matrix. Each row of the matrix represents the state distribution at a time instant.
          I1:         Time profile of concentration of cells infected by 
                      virus 1 - [cell/mL] 
          I2:         Time profile of concentration of cells infected by 
                      virus 2 - [cell/mL] 
          B1:         Time profile of concentration of virus 1 genome 
                      bound to infected cells - [vg/mL] (total conc. in system)
          B2:         Time profile of concentration of virus 2 genome 
                      bound to infected cells - [vg/mL] (total conc. in system)
          B1_pc:      Time profile of concentration of virus 1 genome 
                      bound to infected cells - [vg/mL] (total conc. in system)
          B2_pc:      Time profile of concentration of virus 2 genome 
                      bound to infected cells - [vg/mL] (total conc. in system)
          N1:         Time profile of concentration of virus 1 genome in 
                      nucleus of infected cells - [vg/mL] (total conc. in system)
          N2:         Time profile of concentration of virus 2 genome in 
                      nucleus of infected cells - [vg/mL] (total conc. in system)
          N1_pc:      Time profile of concentration of virus 1 genome in 
                      nucleus of infected cells - [vg/mL] (total conc. in system)
          N2_pc:      Time profile of concentration of virus 2 genome in 
                      nucleus of infected cells - [vg/mL] (total conc. in system)

      States distributed with respect to two infection ages: the output is a matrix. Each row of the matrix represents the state distribution at a time instant.
          Co:         Time profile of concentration of coinfected cells - [cell/mL]
          B1_Co:      Time profile of concentration of virus 1 genome 
                      bound to coinfected cells - [vg/cell] (per cell basis)
          B2_Co:      Time profile of concentration of virus 2 genome 
                      bound to coinfected cells - [vg/cell] (per cell basis)
          B1_Co_pc:   Time profile of concentration of virus 1 genome 
                      bound to coinfected cells - [vg/cell] (per cell basis)
          B2_Co_pc:   Time profile of concentration of virus 2 genome 
                      bound to coinfected cells - [vg/cell] (per cell basis)
          N1_Co:      Time profile of concentration of virus 1 genome in 
                      nucleus of coinfected cells - [vg/cell] (per cell basis)
          N2_Co:      Time profile of concentration of virus 2 genome in 
                      nucleus of coinfected cells - [vg/cell] (per cell basis)
          N1_Co_pc:   Time profile of concentration of virus 1 genome in 
                      nucleus of coinfected cells - [vg/cell] (per cell basis)
          N2_Co_pc:   Time profile of concentration of virus 2 genome in 
                      nucleus of coinfected cells - [vg/cell] (per cell basis)



Additional information on the model input/output structure and on the model parameters to be specified in `def_parameters` is reported in the paper. The names of parameters and variables in the code reflect the symbols used in the papers.

## Expand the simulator
New parameters, variables, and settings can be added to the simulator by editing the simulator functions and scripts. 

For instance, simulation of a continuous feed that contains more than one nutrient or suspended virions/infected cells can be set up by editing the relevant equations in `inf_model` and adding these additional inputs in `RunSimulation` and `main`.

## How to cite
Destro, F. and Braatz, R. D. (2024). Efficient simulation of viral transduction and propagation for biomanufacturing, biorXiv. DOI: 10.1101/2024.03.30.587435_Mol

## License
The code in this repository is provided under a CC BY-NC-ND 4.0 license, as detailed in the `LICENSE` file.

