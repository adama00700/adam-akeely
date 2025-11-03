# adam-akeely
ESP-INFM: Reversed Phonon Nonlinearity in Engineered SiC Lattices
# ESP-INFM: Reversed Phonon Nonlinearity in Engineered SiC Lattices

**Author**: Adam Akeely  
**Independent Researcher**  
**November 2025**

---

## Abstract
A nonlinear field-theoretic model coupling a **pre-engineered SiC lattice** (with Ni impurity) to an effective vacuum-like field via **terahertz Gaussian excitation**.  
**Key result**: *Higher drive amplitude → smaller atomic displacement* — **inverse nonlinearity**.

$$
E_{\rm eff} \approx \alpha \left( \ddot{u} + \beta u^3 \right)
$$

---

## Key Results
| Amplitude | $\Delta z_{\max}$ (Å) |
|---------|----------------|
| 0.10    | 0.591          |
| 0.50    | 0.164          |
| 1.00    | 0.055          |

→ **Negative log-log slope** = field self-stabilization.

---

## Repository Contents
- `lammps.py` – Full simulation script (Colab-ready)
- `SiC.tersoff` – Erhart/Albe potential
- `POSCAR_SiC_Ni_clean` – 6×6×6 supercell
- `disp_ni_*.dat` – Displacement data for A = 0.1, 0.5, 1.0
- `plot_loglog.py` – Generate Fig. 1

---

## How to Run
1. Open in [Google Colab](https://colab.research.google.com/github/adama0700/adam-akeely/blob/main/lammps.py)
2. Run all cells
3. Output: `ni_displacement_plot.png` + `disp_ni.dat`

---

## Paper
Preprint: [arXiv:2501.xxxxx](https://arxiv.org/abs/2501.xxxxx) *(coming soon)*

---

## License
[MIT License](LICENSE) – Free to use, cite, modify.

> *"The lattice doesn't break under pressure. It learns to contain it."* – Adam Akeely
