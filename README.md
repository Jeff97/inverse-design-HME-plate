# Inverse design of hard-magnetoelastic plates for shape morphing

_Analytical formulas, branch selection, Abaqus UEL simulations, post-processing scripts, benchmarks, and experiments for inverse-designed hard-magnetoelastic plates._

[![GitHub stars](https://img.shields.io/github/stars/Jeff97/inverse-design-HME-plate?style=social)](https://github.com/Jeff97/inverse-design-HME-plate/stargazers)
[![License](https://img.shields.io/github/license/Jeff97/inverse-design-HME-plate?color=blue)](LICENSE)
[![Status](https://img.shields.io/badge/status-companion%20manuscript-6f42c1)](#-manuscript)

---

## 📋 Overview

This repository provides the supplementary finite-element and experimental resources for the manuscript _Inverse Design of Hard-Magnetoelastic Plates for Shape-morphing: Analytical Formulas and Experimental Validation_.

<p align="center">
  <img src="figures/Framework.png" width="780" alt="Inverse-design workflow for hard-magnetoelastic cantilever plates">
</p>

_Figure 1: From a prescribed target shape and magnetic loading to feasible magnetization branches, branch selection, and 3D finite-element verification._

The repository supports five connected tasks:

- Verify the forward plate model against three-dimensional finite elements
- Determine whether a target shape is mathematically and physically reachable
- Solve and rank candidate magnetization branches
- Reproduce bending, wavy, letter-like, and sensitivity examples
- Compare analytical design, 3D FE predictions, and experiments

## ⚡ Quick start

### Prerequisites

- Abaqus/Standard with a working Fortran user-subroutine toolchain
- Visual Studio plus Intel oneAPI/Parallel Studio when using the common Windows Abaqus setup
- Abaqus/CAE for the included `.odb` post-processing scripts

### Run a representative inverse-design case

```bat
git clone https://github.com/Jeff97/inverse-design-HME-plate.git
cd inverse-design-HME-plate\InverseCase1-Bending
abaqus job=Beam-Th004-1Pi cpus=8 user=UEL-HardMagneto-FM-Q2.for
```

The `job` value must match an `.inp` filename without its extension, and `user` must point to the corresponding `.for` file in the case directory. The supplied `Command.txt` files contain additional case-specific examples; replace any author-machine paths before use.

## 📚 Repository map

| Path | Purpose |
| --- | --- |
| [`ForwardProblem/`](ForwardProblem/) | Forward simulations and model verification |
| [`Verification-branch/`](Verification-branch/) | Branch selection, reachability, and loading-path checks |
| [`InverseCase1-Bending/`](InverseCase1-Bending/) | Monotonic large-angle bending targets |
| [`InverseCase2-Wave/`](InverseCase2-Wave/) | Distributed wavy targets |
| [`InverseCase3-Letter/`](InverseCase3-Letter/) | S-C-U-T letter-like targets |
| [`Benchmark-Epsilon-Th004/`](Benchmark-Epsilon-Th004/) | Magnetic-loading sensitivity |
| [`Benchmark-Th-Epsilon003/`](Benchmark-Th-Epsilon003/) | Thickness sensitivity |
| [`Benchmark-LinearM0/`](Benchmark-LinearM0/) | Non-uniform magnetization-magnitude study |
| [`Experiment/`](Experiment/) | FE files, scripts, outputs, and experimental comparisons |
| [`figures/`](figures/) | Regenerated manuscript figures |

Most case directories use the following structure:

| Item | Role |
| --- | --- |
| `*.inp` | Abaqus model and analysis definition |
| `*.for` | UEL or related Fortran user subroutine |
| `InputFiles*/` | Mesh and auxiliary fragments included by the main input file |
| `Command.txt` | Example Abaqus submission commands |
| `PythonScript/` | Abaqus/CAE coordinate-extraction scripts |
| `Output/` | Exported coordinates used for comparison and plotting |

## 🔧 Exporting coordinates

The post-processing scripts export deformed bottom-surface coordinates from Abaqus `.odb` files. They use Abaqus/CAE session APIs and therefore must run through Abaqus rather than a standard Python interpreter.

```bat
cd inverse-design-HME-plate\InverseCase1-Bending
abaqus cae noGUI=PythonScript\script-TIMEvsCOORD-Th004.py
```

Before running a script, review these fields:

| Field | What to verify |
| --- | --- |
| `odb_path` | Directory containing the generated `.odb` files |
| `odb_files` | Files selected for extraction |
| `folder_name` | Destination for exported `.csv` data |
| `session.Path(...)` | Bottom-surface node path, especially after mesh changes |

The scripts extract deformed `COORD` components along the selected path and write paired coordinate data under the corresponding `Output/` directory.

## 📊 Design cases and figures

<details>
<summary><strong>📊 Core theory and inverse-design cases</strong></summary>

### Forward and inverse kinematics

| Forward problem | Inverse problem |
| --- | --- |
| <img src="figures/Kinematics-Forward.png" width="430" alt="Forward hard-magnetoelastic plate problem"> | <img src="figures/Kinematics-Inverse.png" width="430" alt="Inverse hard-magnetoelastic plate problem"> |

_Figures 2–3: Forward prediction with prescribed magnetization and inverse recovery of the magnetization distribution from a target shape._

### Programmable domain and branch selection

| Programmable endpoint workspace | Candidate branch energy |
| --- | --- |
| <img src="figures/ProgrammableEndpoint.png" width="430" alt="Programmable endpoint workspace"> | <img src="figures/OptimalSolution_Energy.png" width="430" alt="Magnetic potential energy of candidate branches"> |

_Figures 4–5: Feasibility boundaries and energy-based ranking of analytical branches._

### Bending and wave targets

| Large-angle bending | Distributed wave |
| --- | --- |
| <img src="figures/Example1-4.png" width="430" alt="Inverse-designed two-pi bending case"> | <img src="figures/Example2-3.png" width="430" alt="Inverse-designed wavy plate case"> |

_Figures 6–7: Representative inverse-design results; gray, blue, and green profiles denote initial, target, and FE-simulated configurations._

### Letter-like targets

<p align="center">
  <img src="figures/SCUT.png" width="600" alt="Assembled SCUT target layout">
</p>

_Figure 8: Assembled S-C-U-T target layout. Individual FE comparisons are available as `figures/Example3-1.png` through `figures/Example3-4.png`._

</details>

---

<details>
<summary><strong>📊 Benchmarks, verification, and experiments</strong></summary>

### Sensitivity benchmarks

| Magnetic loading | Plate thickness | Magnetization magnitude |
| --- | --- | --- |
| <img src="figures/Benchmark-1-Error.png" width="280" alt="RMSE across magnetic-loading benchmark cases"> | <img src="figures/Benchmark-2-Error.png" width="280" alt="RMSE across thickness benchmark cases"> | <img src="figures/Benchmark-3-Error.png" width="280" alt="RMSE for non-uniform magnetization benchmark cases"> |

_Figures 9–11: Root-mean-square errors for loading, thickness, and magnetization-magnitude benchmarks. Individual configurations are stored alongside each summary figure._

### Forward-model and loading-path verification

| Forward verification | Constitutive comparison | Branch reachability |
| --- | --- | --- |
| <img src="figures/Forward-Case3.png" width="280" alt="Representative forward-model verification"> | <img src="figures/FR_Model_RMSE_BarChart.png" width="280" alt="Difference between F-based and rotation-based simulations"> | <img src="figures/OptimalSolution_Simulation.png" width="280" alt="Branch reachability under different loading paths"> |

_Figures 12–14: Plate-model verification, constitutive-model comparison, and numerical reachability of the analytical branches._

### Experimental comparison

<p align="center">
  <img src="figures/Experiment.png" width="760" alt="Experimental workflow for hard-magnetoelastic plates">
</p>

_Figure 15: Material characterization, inverse design, FE verification, fabrication, fixture-assisted magnetization, and magnetic loading._

<p align="center">
  <img src="figures/Experiment-SCUT.png" width="760" alt="Experimental realization of inverse-designed SCUT shapes">
</p>

_Figure 16: Analytical branches, FE predictions, and experimental realizations for the S-C-U-T targets._

</details>

---

## 🔍 Reproduction notes

- Use the input/subroutine pair supplied within one case directory
- Replace absolute paths in `Command.txt` and Python scripts with paths to your clone
- Keep the intended loading path when evaluating analytical branch reachability
- Update `session.Path(...)` if the mesh or node labels change
- Compare exported bottom-surface coordinates with the matching files under `Output/`
- Treat the supplied models as research reproduction files and review units, constraints, material parameters, magnetic fields, and state variables before reuse

## ✍️ Manuscript

> Z. Li, Y. Wang, J. Wang, X. Yao, M. Hossain, J. Lan, Z. Wang, and J. Wu, “Inverse Design of Hard-Magnetoelastic Plates for Shape-morphing: Analytical Formulas and Experimental Validation,” manuscript.

No DOI, journal, acceptance, or formal publication status is asserted in this repository. Replace this entry with the published citation when one becomes available.

```bibtex
@unpublished{LiInverseDesignHMEPlate,
  title  = {Inverse Design of Hard-Magnetoelastic Plates for Shape-morphing: Analytical Formulas and Experimental Validation},
  author = {Li, Zhanfeng and Wang, Yafei and Wang, Jiong and Yao, Xiaohu and Hossain, Mokarram and Lan, Junjie and Wang, Zuodong and Wu, Jianbin},
  note   = {Manuscript and supplementary files}
}
```

If this repository supports your work, consider starring it and citing the manuscript or its future published version.

## 🔐 License

The source code and documentation are available under the [Apache License 2.0](LICENSE). Manuscript text, publication figures, experimental photographs, and third-party material may be subject to separate terms.
