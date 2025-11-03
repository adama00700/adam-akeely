# ESP-INFM: Reversed Phonon Nonlinearity in Engineered SiC Lattices
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by/4.0/)

**© 2025 Adam Akeely. All rights reserved.**  
**Independent Researcher**  
**November 2025**

---

## Abstract
A **nonlinear field-theoretic model** couples a **pre-engineered 6×6×6 SiC lattice with Ni impurity** to an **effective vacuum-like field** via **terahertz Gaussian excitation**.

**Key result**:  
> **Higher drive amplitude → smaller atomic displacement**

This is **not numerical instability** — it is **field self-stabilization** via:
$$
E_{\rm eff} \approx \alpha \left( \ddot{u} + \beta u^3 \right)
$$

---

## Key Results
| Amplitude | $\Delta z_{\max}$ (Å) |
|-----------|-----------------------|
| 0.10      | 0.591                 |
| 0.50      | 0.164                 |
| 1.00      | 0.055                 |

→ **Negative log-log slope** = **inverse nonlinearity** = **lattice-mediated field saturation**

---

## Repository Contents
- `lammps.py` – Full simulation script (Google Colab ready)
- `SiC.tersoff` – Erhart/Albe Tersoff potential
- `POSCAR_SiC_Ni_clean` – VASP structure file
- `disp_ni_*.dat` – Raw displacement data (3 amplitudes)
- `plot_loglog.py` – Generate Fig. 1 (optional)

---

## How to Run
1. Open in [Google Colab](https://colab.research.google.com/github/adama0700/adam-akeely/blob/main/lammps.py)
2. Run all cells
3. Output: `ni_displacement_plot.png` + `disp_ni.dat`

---

## Citation
```bibtex
@article{akeely2025espinfm,
  title={ESP-INFM: Reversed Phonon Nonlinearity in Engineered SiC Lattices},
  author={Akeely, Adam},
  journal={arXiv preprint arXiv:2501.xxxxx},
  year={2025},
  note={\url{https://github.com/adama0700/adam-akeely}}
}
