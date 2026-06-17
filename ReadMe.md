# Inverse Design of Hard-Magnetoelastic Plates for Shape-morphing: Analytical Formulas and Experimental Validation

This repository provides supplementary ABAQUS files for the manuscript:

> **Inverse Design of Hard-Magnetoelastic Plates for Shape-morphing: Analytical Formulas and Experimental Validation**.
> Zhanfeng Li, Yafei Wang, Jiong Wang, Xiaohu Yao, Mokarram Hossain, Junjie Lan, Zuodong Wang, and Jianbin Wu.  

The files include ABAQUS input files, user-element subroutines, command examples, post-processing scripts, and selected output data for reproducing the forward-model validations, analytical inverse-design examples, benchmark studies, branch-selection checks, and experimental comparisons reported in the manuscript.

## Repository Contents

- `ForwardProblem/`: forward simulations and model-verification cases.
- `Verification-branch/`: numerical checks for branch selection, branch reachability, and loading-path effects.
- `InverseCase1-Bending/`: inverse-design examples for monotonic large-angle bending.
- `InverseCase2-Wave/`: inverse-design examples for distributed wavy morphing.
- `InverseCase3-Letter/`: inverse-design examples for S-C-U-T letter-like targets.
- `Benchmark-Epsilon-Th004/`, `Benchmark-Th-Epsilon003/`, `Benchmark-LinearM0/`: sensitivity and non-uniform magnetization-magnitude benchmark cases.
- `Experiment/`: FE files, scripts, outputs, and snapshots related to the experimental S-C-U-T demonstrations.

Most case folders follow the same structure:

- `*.inp`: ABAQUS input files.
- `*.for`: ABAQUS UEL or related user subroutines.
- `InputFiles*/`: mesh and auxiliary input fragments included by the main `.inp` files.
- `Command.txt`: example commands for submitting the jobs from an ABAQUS Command window.
- `PythonScript/`: scripts for extracting deformed coordinates from `.odb` result files.
- `Output/`: exported coordinate data used for plotting or comparison.

## Simulation

Running these examples requires ABAQUS/Standard with a working Fortran user-subroutine toolchain. On Windows, make sure Visual Studio and Intel oneAPI or Intel Parallel Studio are correctly linked to ABAQUS before submitting jobs. One possible linking guide is available [here](https://www.researchgate.net/publication/349991987_Linking_ABAQUS_20192020_and_Intel_oneAPI_Base_Toolkit_FORTRAN_Compiler).

Submit jobs from the ABAQUS Command window:

1. Open an ABAQUS Command window.
2. Change to the target case folder. The path shown in each `Command.txt` is the author's local working path, so replace it with the path to your local clone.
3. Copy the corresponding `abaqus job=... user=...` command from `Command.txt`.
4. Run the command in the case folder. If ABAQUS asks whether to overwrite an existing job, answer `y`.

Example:

```bat
cd /d D:\path\to\inverse-design-HME-plate\InverseCase1-Bending
abaqus job=Beam-Th004-1Pi cpus=8 user=UEL-HardMagneto-FM-Q2.for
```

The `job` name should match the target `.inp` file without the extension, and the `user` argument should point to the corresponding `.for` subroutine in the same case folder.

## Exporting Coordinates

The post-processing scripts in `PythonScript/` export deformed bottom-surface coordinates from ABAQUS `.odb` files. They use ABAQUS/CAE session APIs, so run them with `abaqus cae noGUI`, not with a standard Python interpreter.

Example:

```bat
cd /d D:\path\to\inverse-design-HME-plate\InverseCase1-Bending
abaqus cae noGUI=PythonScript\script-TIMEvsCOORD-Th004.py
```

Before running a script, check these fields inside the script:

- `odb_path`: folder containing the generated `.odb` files.
- `odb_files`: list of `.odb` files to process.
- `folder_name`: output folder for the exported `.csv` files.
- `session.Path(...)`: the node list defining the extraction path. Update it if the mesh or model changes.

The scripts extract the deformed `COORD` components along the selected bottom path and write paired coordinate data to `Output/*.csv`.

## Manuscript Figures

The manuscript figure files referenced by the LaTeX source were regenerated as PNG files and placed under `figures/`.

### Kinematics

| Forward problem | Inverse problem |
| --- | --- |
| <img src="figures/Kinematics-Forward.png" width="430" alt="Forward problem"> | <img src="figures/Kinematics-Inverse.png" width="430" alt="Inverse problem"> |

**Kinematics and forward/inverse problem setting.** (a) A plate with prescribed remanent magnetization $\mathbf{M}$ deforms from the reference configuration $\mathcal{K}_r$ to the current configuration $\mathcal{K}_t$ under the applied field $\mathbf{H}_a$. (b) The target shape and applied field are prescribed, and the required magnetization distribution $\mathbf{M}$ is solved.

### Programmable Domain

| $\bar{h}=0.01$ | $\bar{h}=0.02$ | $\bar{h}=0.04$ |
| --- | --- | --- |
| <img src="figures/ProgrammableRegion-1.png" width="300" alt="Programmable domain for hbar 0.01"> | <img src="figures/ProgrammableRegion-2.png" width="300" alt="Programmable domain for hbar 0.02"> | <img src="figures/ProgrammableRegion-3.png" width="300" alt="Programmable domain for hbar 0.04"> |

<p align="center">
  <img src="figures/ProgrammableEndpoint.png" width="520" alt="Endpoint workspace for hbar 0.02 and epsilon 0.1">
</p>

**Programmable domain for the constant-magnetization case.** (a)--(c) The green-tinted lower boundary surface denotes the critical mathematical boundary $\epsilon_c$, the red plane denotes the physical upper bound $\epsilon_{\max}=0.3$, and the green shaded volume indicates the region that is both mathematically solvable and physically reachable. (d) Programmable endpoint workspace.

### Branch Selection

| Distribution of $\theta_M$ | Distribution of $\psi_m$ |
| --- | --- |
| <img src="figures/OptimalSolution_ThetaM.png" width="340" alt="Distribution of theta_M"> | <img src="figures/OptimalSolution_Energy.png" width="340" alt="Distribution of psi_m"> |

| Profile of branch 1 | Profile of branch 2 |
| --- | --- |
| <img src="figures/OptimalSolution_Sample1.png" width="420" alt="Profile of branch 1"> | <img src="figures/OptimalSolution_Sample2.png" width="420" alt="Profile of branch 2"> |

**Representative branch-selection example based on magnetic potential-energy ranking.** (a) Two analytical branches of $\theta_M(\xi)$, where the red curve denotes branch 1 and the blue curve denotes branch 2. (b) Corresponding magnetic potential-energy density distributions. (c) Magnetization vectors for branch 1 in the reference and current configurations. (d) Magnetization vectors for branch 2 in the reference and current configurations.

### Inverse-Design Workflow

<p align="center">
  <img src="figures/Framework.png" width="760" alt="Inverse-design workflow">
</p>

**Inverse-design workflow for shape-morphing of hard-magnetoelastic cantilever plates.** For a target shape $g(\xi)$ and inputs $(L,2h,C_0,H_0,\theta_H,M_0)$, the workflow checks feasibility ($\Delta(\xi)\ge 0$ and loading/magnetization bounds), solves candidate $\theta_M$ branches, selects the preferred branch, and verifies the programmed shape through 3D FE simulations.

### Case 1: Different Bending Angles

| Case 1-1 with $\Theta_{\mathrm{end}}=\pi/2$. | Case 1-2 with $\Theta_{\mathrm{end}}=\pi$. |
| --- | --- |
| <img src="figures/Example1-1.png" width="420" alt="Case 1-1"> | <img src="figures/Example1-2.png" width="420" alt="Case 1-2"> |

| Case 1-3 with $\Theta_{\mathrm{end}}=3\pi/2$. | Case 1-4 with $\Theta_{\mathrm{end}}=2\pi$. |
| --- | --- |
| <img src="figures/Example1-3.png" width="420" alt="Case 1-3"> | <img src="figures/Example1-4.png" width="420" alt="Case 1-4"> |

**Case 1 inverse-design results for increasing end rotation.** Gray profiles denote the initial shapes, blue profiles denote the prescribed target shapes, and green profiles denote the FE-simulated deformation sequence.

### Case 2: Different Waviness Levels

| Case 2-1 with $A=\pi/2$ and $n=1/2$. | Case 2-2 with $A=\pi/4$ and $n=3/2$. |
| --- | --- |
| <img src="figures/Example2-1.png" width="420" alt="Case 2-1"> | <img src="figures/Example2-2.png" width="420" alt="Case 2-2"> |

| Case 2-3 with $A=\pi/3$ and $n=5/2$. | Case 2-4 with $A=\pi/6$ and $n=7/2$. |
| --- | --- |
| <img src="figures/Example2-3.png" width="420" alt="Case 2-3"> | <img src="figures/Example2-4.png" width="420" alt="Case 2-4"> |

**Case 2 inverse-design results for different waviness levels.** Gray profiles denote the initial shapes, blue profiles denote the prescribed target shapes, and green profiles denote the FE-simulated deformation sequence.

### Case 3: Letter-Like Target Shapes

<p align="center">
  <img src="figures/SCUT.png" width="560" alt="Assembled SCUT layout">
</p>

| Case 3-1 for letter S. | Case 3-2 for letter C. |
| --- | --- |
| <img src="figures/Example3-1.png" width="420" alt="Case 3-1 for letter S"> | <img src="figures/Example3-2.png" width="420" alt="Case 3-2 for letter C"> |

| Case 3-3 for letter U. | Case 3-4 for letter T. |
| --- | --- |
| <img src="figures/Example3-3.png" width="420" alt="Case 3-3 for letter U"> | <img src="figures/Example3-4.png" width="420" alt="Case 3-4 for letter T"> |

**Case 3 inverse-design and simulation results for the assembled "SCUT" pattern.** Gray profiles denote the initial shapes, blue profiles denote the prescribed target shapes, and green profiles denote the FE-simulated deformation sequence.

### Sensitivity Benchmarks

| $\epsilon=0.02$ | $\epsilon=0.05$ | $\epsilon=0.10$ | $\epsilon=0.15$ |
| --- | --- | --- | --- |
| <img src="figures/Benchmark-1-1.png" width="230" alt="Benchmark epsilon 0.02"> | <img src="figures/Benchmark-1-2.png" width="230" alt="Benchmark epsilon 0.05"> | <img src="figures/Benchmark-1-3.png" width="230" alt="Benchmark epsilon 0.10"> | <img src="figures/Benchmark-1-4.png" width="230" alt="Benchmark epsilon 0.15"> |

<p align="center">
  <img src="figures/Benchmark-1-Error.png" width="380" alt="Root-mean-square error for benchmark B_epsilon">
</p>

**Shape-morphing results for benchmark $\mathcal{B}_{\epsilon}$ evaluated at different magnetic-loading levels $\epsilon$.** In (a)--(d), the dashed curves denote the deformed bottom-surface profiles obtained from 3D FE simulations for $\epsilon=0.02$, $0.05$, $0.10$, and $0.15$, respectively; (e) the root-mean-square error.

| $2\bar{h}=0.01$ | $2\bar{h}=0.03$ | $2\bar{h}=0.05$ | $2\bar{h}=0.07$ |
| --- | --- | --- | --- |
| <img src="figures/Benchmark-2-1.png" width="230" alt="Benchmark hbar 0.01"> | <img src="figures/Benchmark-2-2.png" width="230" alt="Benchmark hbar 0.03"> | <img src="figures/Benchmark-2-3.png" width="230" alt="Benchmark hbar 0.05"> | <img src="figures/Benchmark-2-4.png" width="230" alt="Benchmark hbar 0.07"> |

<p align="center">
  <img src="figures/Benchmark-2-Error.png" width="380" alt="Root-mean-square error for benchmark B_h">
</p>

**Shape-morphing results for benchmark $\mathcal{B}_{\bar{h}}$ evaluated across different total normalized thicknesses $2\bar{h}$.** In (a)--(d), the dashed curves denote the deformed bottom-surface profiles obtained from 3D FE simulations for $2\bar{h}=0.01$, $0.03$, $0.05$, and $0.07$, respectively; (e) the root-mean-square error.

| $a=1,\ b=0$ (constant distribution) | $a=0.5,\ b=0.5$ (linearly increasing) |
| --- | --- |
| <img src="figures/Benchmark-3-1.png" width="420" alt="Constant magnetization magnitude distribution"> | <img src="figures/Benchmark-3-2.png" width="420" alt="Linearly increasing magnetization magnitude distribution"> |

| $a=0.3,\ b=0.04$ (near-minimal dosage) | Root-mean-square error |
| --- | --- |
| <img src="figures/Benchmark-3-3.png" width="420" alt="Near-minimal dosage magnetization magnitude distribution"> | <img src="figures/Benchmark-3-Error.png" width="420" alt="Root-mean-square error for non-uniform magnetization magnitude"> |

**Illustrative shape-morphing results for non-uniform magnetization magnitude.** In (a)--(c), the dashed curves denote the deformed bottom-surface profiles obtained from 3D FE simulations for constant, linearly increasing, and near-minimal-dosage magnetization-magnitude distributions, respectively; (d) the root-mean-square error.

### Experimental Comparison

<p align="center">
  <img src="figures/Experiment.png" width="760" alt="Experimental workflow">
</p>

**Experimental workflow for material characterization, theoretical design, fixture-assisted magnetization, and magnetic loading.** The workflow includes VSM and uniaxial tests, inverse-design and FE verification, fabrication and magnetization, and comparison between the loaded profile and the target shape.

<p align="center">
  <img src="figures/Experiment-SCUT.png" width="760" alt="Experimental realization of SCUT shapes">
</p>

**Experimental realization of S-C-U-T shapes.** Each row corresponds to one target profile, from S to T. The left column shows the inverse-design branches, selected magnetization branch, auxiliary magnetizing fixture, and resulting flat magnetized sample. The middle column shows the corresponding finite element predictions under the prescribed upward actuation field, with the color map denoting the von Mises stress $\sigma_M$. The right column shows the experimental loading results under the uniform external field $\mathbf{H}_a$.

### Verification Figures

| Case 1 | Case 2 | Case 3 |
| --- | --- | --- |
| <img src="figures/Forward-Case1.png" width="300" alt="Forward verification case 1"> | <img src="figures/Forward-Case2.png" width="300" alt="Forward verification case 2"> | <img src="figures/Forward-Case3.png" width="300" alt="Forward verification case 3"> |

| Case 4 | Case 5 |
| --- | --- |
| <img src="figures/Forward-Case4.png" width="300" alt="Forward verification case 4"> | <img src="figures/Forward-Case5.png" width="300" alt="Forward verification case 5"> |

**Comparison between the plate model and 3D finite element simulations for representative cases.** Solid lines denote the plate-model predictions, and scatter markers denote the 3D FE results extracted from the bottom surface of the plate.

| $\mathbb{F}$-based and rotation-based profiles | Model discrepancy |
| --- | --- |
| <img src="figures/Forward-Case6.png" width="420" alt="F-based and rotation-based model comparison"> | <img src="figures/FR_Model_RMSE_BarChart.png" width="420" alt="Difference between F-based and rotation-based model simulations"> |

**Constitutive-model comparison for the magnetization update.** (a) Comparison of the deformed profiles computed from the $\mathbb{F}$-based formulation (solid lines) and the rotation-based 3D FE implementation (scattered points) for varying $C_0$. (b) Difference between the $\mathbb{F}$- and rotation-based model simulations.

<p align="center">
  <img src="figures/OptimalSolution_Simulation.png" width="560" alt="Numerical simulations of the two branches under different loading paths">
</p>

**Numerical simulations of the two branches under different loading paths.** The colored 3D FE configurations show the magnetic-potential density $\psi_m$. Branch 1 does not reach the prescribed target under the default path, but reaches it under the specified path; branch 2 reaches the target under the default path.

<p align="center">
  <img src="figures/ExpSetupDetail.png" width="760" alt="Experimental implementation details">
</p>

**Experimental implementation details.** Fabrication by curing a flat silicone plate, fixture-assisted magnetization in a 3D-printed auxiliary mold under the pulsed field $\mathbf{B}_m$, and magnetic loading of the clamped sample under the applied field $\mathbf{H}_a$.

## License

This repository is released under the Apache License 2.0. See [LICENSE](LICENSE) for details.
