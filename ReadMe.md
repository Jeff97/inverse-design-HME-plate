# Inverse Design of Hard-Magnetoelastic Plates for Shape Morphing

This repository provides supplementary ABAQUS files for the manuscript:

> **Inverse design of hard-magnetoelastic plates for shape morphing under plane-strain conditions**.
> Zhanfeng Li, Yafei Wang, Junjie Lan, Zuodong Wang, Jianbin Wu, Xiaohu Yao, Mokarram Hossain, and Jiong Wang.  

The files include ABAQUS input files, user-element subroutines, command examples, post-processing scripts, and selected output data for the forward-model validations, inverse-design examples, benchmark studies, and experimental comparisons reported in the manuscript.

## Repository Contents

- `ForwardProblem/`: forward simulations and model-verification cases.
- `Verification-branch/`: numerical checks for branch selection and loading-path effects.
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

## Representative Figures

The PDF figures used below were converted from the manuscript figure files to PNG and placed in this repository root.

### Problem Setting

![Kinematics and forward/inverse problem setting for the plane-strain plate.](Kinematics.png)

**Kinematics and forward/inverse problem setting for the plane-strain plate.** In the forward problem, a plate of length $L$, width $W$, and thickness $2h$ with programmed remanent magnetization $\mathbf{M}$ deforms from the reference configuration to the current configuration under the applied field $\mathbf{H}_a$. In the inverse problem, the target shape and applied field are prescribed, and the required magnetization distribution $\mathbf{M}$ is solved.

### Inverse-Design Workflow

![Inverse-design workflow for shape programming of hard-magnetoelastic cantilever plates.](Framework.png)

**Inverse-design workflow for shape programming of hard-magnetoelastic cantilever plates.** For a target shape $g(\xi)$ and inputs $(L,2h,C_0,H_0,\theta_H,M_0)$, the workflow checks feasibility, solves candidate magnetization-direction branches, selects the preferred branch by magnetic-potential ranking, and verifies the result by 3D FE simulations.

### Programmable Domain

**Programmable domain for the constant-magnetization case under fabrication/loading limits.** The mathematical limit defines the local solvability boundary, the physical limit gives the upper bound of the available magnetic loading, and the green region indicates the portion that is both mathematically solvable and physically reachable.

<table align="center">
  <tr>
    <th>h = 0.01</th>
    <th>h = 0.02</th>
    <th>h = 0.04</th>
  </tr>
  <tr>
    <td align="center"><img src="ProgrammableRegion-1.png" width="300" alt="Programmable domain for hbar 0.01"></td>
    <td align="center"><img src="ProgrammableRegion-2.png" width="300" alt="Programmable domain for hbar 0.02"></td>
    <td align="center"><img src="ProgrammableRegion-3.png" width="300" alt="Programmable domain for hbar 0.04"></td>
  </tr>
</table>

<p align="center">
  <img src="ProgrammableEndpoint.png" width="420" alt="Endpoint workspace for hbar 0.02 and epsilon 0.1">
</p>

<p align="center"><strong>Endpoint workspace for h = 0.02 and epsilon = 0.1.</strong></p>

The gray/blue surface denotes the critical mathematical boundary $\epsilon_c$, the red plane denotes the physical upper bound $\epsilon_{\max}=0.3$, and the green shaded volume indicates the feasible programmable region. The endpoint workspace is obtained by stochastic inverse-design sampling in the upper half-plane.

### Forward-Model Verification

**Comparison between the asymptotic plate model and 3D finite element simulations for representative cases.** Solid lines denote the plate-model predictions, and scatters denote the 3D FE results. Subfigures correspond to Cases 1-5.

<table align="center">
  <tr>
    <th>Case 1</th>
    <th>Case 2</th>
    <th>Case 3</th>
  </tr>
  <tr>
    <td align="center"><img src="Forward-Case1.png" width="300" alt="Forward verification case 1"></td>
    <td align="center"><img src="Forward-Case2.png" width="300" alt="Forward verification case 2"></td>
    <td align="center"><img src="Forward-Case3.png" width="300" alt="Forward verification case 3"></td>
  </tr>
</table>

<table align="center">
  <tr>
    <th>Case 4</th>
    <th>Case 5</th>
  </tr>
  <tr>
    <td align="center"><img src="Forward-Case4.png" width="300" alt="Forward verification case 4"></td>
    <td align="center"><img src="Forward-Case5.png" width="300" alt="Forward verification case 5"></td>
  </tr>
</table>

### Case 1: Different Bending Angles

**Case-1 illustrative examples arranged in a $2\times2$ layout.**

| Case 1-1 with $\Theta_{\mathrm{end}}=\pi/2$ | Case 1-2 with $\Theta_{\mathrm{end}}=\pi$ |
| --- | --- |
| <img src="Example1-1.png" width="420"> | <img src="Example1-2.png" width="420"> |

| Case 1-3 with $\Theta_{\mathrm{end}}=3\pi/2$ | Case 1-4 with $\Theta_{\mathrm{end}}=2\pi$ |
| --- | --- |
| <img src="Example1-3.png" width="420"> | <img src="Example1-4.png" width="420"> |

### Case 2: Different Waviness Levels

**Case-2 illustrative examples arranged in a $2\times2$ layout.**

| Case 2-1 with $A=\pi/2$ and $n=1/2$ | Case 2-2 with $A=\pi/4$ and $n=3/2$ |
| --- | --- |
| <img src="Example2-1.png" width="420"> | <img src="Example2-2.png" width="420"> |

| Case 2-3 with $A=\pi/3$ and $n=5/2$ | Case 2-4 with $A=\pi/6$ and $n=7/2$ |
| --- | --- |
| <img src="Example2-3.png" width="420"> | <img src="Example2-4.png" width="420"> |

### Case 3: Letter-Like Target Shapes

**Case-3 results showing an assembled "SCUT" overview and the four corresponding single-letter targets.**

<p align="center">
  <img src="SCUT.png" width="520" alt="Assembled SCUT layout">
</p>

**Assembled "SCUT" layout from the four deformed samples using rigid-body placement.**

| Case 3-1 for letter S | Case 3-2 for letter C |
| --- | --- |
| <img src="Example3-1.png" width="420"> | <img src="Example3-2.png" width="420"> |

| Case 3-3 for letter U | Case 3-4 for letter T |
| --- | --- |
| <img src="Example3-3.png" width="420"> | <img src="Example3-4.png" width="420"> |

### Sensitivity Benchmarks

**Shape programming results for benchmark $\mathcal{B}_{\epsilon}$ evaluated at different magneto-mechanical parameters $\epsilon$.**

| $\epsilon=0.02$ | $\epsilon=0.05$ | $\epsilon=0.10$ | $\epsilon=0.15$ |
| --- | --- | --- | --- |
| <img src="Benchmark-1-1.png" width="230"> | <img src="Benchmark-1-2.png" width="230"> | <img src="Benchmark-1-3.png" width="230"> | <img src="Benchmark-1-4.png" width="230"> |

<p align="center">
  <img src="Benchmark-1-Error.png" width="380" alt="Root-mean-square error for benchmark B_epsilon">
</p>


**Shape programming results for benchmark $\mathcal{B}_{\bar{h}}$ evaluated across different total normalized thicknesses $2\bar{h}$.**

| $2\bar{h}=0.01$ | $2\bar{h}=0.03$ | $2\bar{h}=0.05$ | $2\bar{h}=0.07$ |
| --- | --- | --- | --- |
| <img src="Benchmark-2-1.png" width="230"> | <img src="Benchmark-2-2.png" width="230"> | <img src="Benchmark-2-3.png" width="230"> | <img src="Benchmark-2-4.png" width="230"> |

<p align="center">
  <img src="Benchmark-2-Error.png" width="380" alt="Root-mean-square error for benchmark B_h">
</p>


**Illustrative shape programming results for non-uniform magnetization magnitude.** Panels compare the target profile with the simulated profiles and show the corresponding position errors.

| $a=1,\ b=0$ | $a=0.5,\ b=0.5$ |
| --- | --- |
| <img src="Benchmark-3-1.png" width="420"> | <img src="Benchmark-3-2.png" width="420"> |

| $a=0.3,\ b=0.04$ | Position error |
| --- | --- |
| <img src="Benchmark-3-3.png" width="420"> | <img src="Benchmark-3-Error.png" width="420"> |

### Experimental Comparison

![Experimental workflow for material characterization, theoretical design, fixture-assisted magnetization, and magnetic loading.](Experiment.png)

**Experimental workflow for material characterization, theoretical design, fixture-assisted magnetization, and magnetic loading.** The workflow includes VSM and uniaxial tests, inverse-design and FE verification, fabrication and magnetization using an auxiliary mold, and comparison between the loaded experimental profile and the target shape.

![Experimental realization of S-C-U-T shape programming.](Experiment-SCUT.png)

**Experimental realization of S-C-U-T shape programming.** Each row corresponds to one target profile, from S to T. The left column shows the inverse-design branches, selected magnetization branch, auxiliary magnetizing fixture, and resulting flat magnetized sample. The middle column shows finite element predictions under the prescribed upward actuation field. The right column shows the experimental loading results under the uniform external field.

## License

This repository is released under the Apache License 2.0. See [LICENSE](LICENSE) for details.
