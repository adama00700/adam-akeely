import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from google.colab import files
import sys

print("="*80)
print("Project Harmony/Master - Phonon Hammer with Tersoff + NVT")
print("FIX v55: Final stable version — correct fix print syntax, remove unsupported fixes.")
print("="*80)

# === 1. تثبيت Conda (Miniforge) في Colab ===
print("\n[1/9] Installing Conda (Miniforge)...")
try:
    import condacolab
    print("CondaColab is already installed and active.")
    if 'google.colab' in sys.modules and not os.environ.get("CONDA_PREFIX"):
        print("... CondaColab imported but not active. Re-installing to force kernel switch.")
        condacolab.install()
        print("Kernel will restart. RERUN THIS CELL after the kernel restarts.")
        exit()
except ImportError:
    print("... CondaColab not found. Installing now.")
    import condacolab
    condacolab.install()
    print("CondaColab installed. Kernel will restart.")
    print("RERUN THIS CELL after the kernel restarts.")
    exit()
except Exception as e:
    print(f"Unexpected Conda error: {e}")
    import condacolab
    condacolab.install()
    print("Kernel will restart. RERUN THIS CELL after the kernel restarts.")
    exit()

try:
    import condacolab
    if not os.environ.get("CONDA_PREFIX"):
        print("CondaColab failed to activate after restart. Please re-run the cell manually.")
        exit()
    print("CondaColab check passed (post-restart).")
except:
    print("CondaColab failed to import. Please re-run the cell.")
    exit()

# === 2. تثبيت LAMMPS و OpenMPI ===
print("\n[2/9] Installing modern LAMMPS (MPI build) and 'openmpi'...")
try:
    subprocess.run(["conda", "install", "-c", "conda-forge", "lammps", "openmpi", "-y"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    print("Modern 'lammps' (MPI) and 'openmpi' (mpirun) installed.")
except Exception as e:
    print(f"Conda installation failed: {e}")
    exit()

# === 3. تحديد مسارات lmp و mpirun ===
print("\n[3/9] Locating LAMMPS and mpirun executables...")
lmp_path = subprocess.run(['which', 'lmp'], capture_output=True, text=True).stdout.strip()
mpirun_path = subprocess.run(['which', 'mpirun'], capture_output=True, text=True).stdout.strip()

if not lmp_path or not mpirun_path:
    print(f"Critical executable missing. lmp: '{lmp_path}', mpirun: '{mpirun_path}'")
    exit()

lammps_cmd = [mpirun_path, "-np", "1", "--allow-run-as-root", lmp_path]
print(f"LAMMPS executable found: {lmp_path}")
print(f"mpirun executable found: {mpirun_path}")
print(f"Using command: mpirun -np 1 --allow-run-as-root {lmp_path}")

# === 4. إنشاء ملف SiC.tersoff ===
print("\n[4/9] Creating SiC.tersoff...")
tersoff_content = """# SiC Tersoff Potential - Erhart & Albe 2005
Si  Si  Si  1.0  1.0  0.0  3.8049e4  1.6217e1  -0.59825  0.78734  1.099e-6  1.7322e1  4.7118e2  2.85  0.15  2.4799e1  1.8308e3
C   C   C   1.0  1.0  0.0  1.0139e5  1.6146e1  -0.48489  0.72751  1.1000e-6  1.7322e1  4.3844e2  1.95  0.15  3.4647e1  1.3936e3
Si  Si  C   1.0  1.0  0.0  3.8049e4  1.6217e1  -0.48189  0.72751  1.0999e-6  1.7322e1  4.7118e2  2.85  0.15  2.4799e1  1.8308e3
Si  C   Si  1.0  1.0  0.0  3.8049e4  1.6217e1  -0.59825  0.78734  1.0999e-6  1.7322e1  4.7118e2  2.85  0.15  2.4799e1  1.8308e3
C   Si  Si  1.0  1.0  0.0  1.0139e5  1.6146e1  -0.48489  0.72751  1.1000e-6  1.7322e1  4.3844e2  1.95  0.15  3.4647e1  1.3936e3
C   Si  C   1.0  1.0  0.0  1.0139e5  1.6146e1  -0.48489  0.72751  1.1000e-6  1.7322e1  4.3844e2  1.95  0.15  3.4647e1  1.3936e3
C   C   Si  1.0  1.0  0.0  1.0139e5  1.6146e1  -0.48489  0.72751  1.1000e-6  1.7322e1  4.3844e2  1.95  0.15  3.4647e1  1.3936e3
Si  C   C   1.0  1.0  0.0  3.8049e4  1.6217e1  -0.59825  0.78734  1.0999e-6  1.7322e1  4.7118e2  2.85  0.15  2.4799e1  1.8308e3
"""
with open("SiC.tersoff", "w") as f:
    f.write(tersoff_content)
print("SiC.tersoff created.")

# === 5. إنشاء data.sic_ni (تم تصحيح تحديد موقع Ni) ===
print("\n[5/9] Creating data.sic_ni...")
a = 4.36
n = 6      # <-- عدل من 4 إلى 6 للحصول على 1728 ذرة
L = n * a
si_basis = np.array([[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) * a
c_basis = si_basis + np.array([0.25,0.25,0.25]) * a
atoms = []
atom_id = 1
center = np.array([L/2, L/2, L/2])
ni_id = None
min_dist = float('inf')  # تم التصحيح: استخدام inf بدلاً من 999

for i in range(n):
    for j in range(n):
        for k in range(n):
            shift = np.array([i,j,k]) * a
            # ذرات السيليكون (Si)
            for pos in si_basis:
                r = shift + pos
                dist = np.linalg.norm(r - center)
                atoms.append((atom_id, 1, *r))  # type 1 = Si
                if dist < min_dist:
                    min_dist = dist
                    ni_id = atom_id
                atom_id += 1
            # ذرات الكربون (C)
            for pos in c_basis:
                r = shift + pos
                atoms.append((atom_id, 2, *r))  # type 2 = C
                atom_id += 1

# استبدال الذرة الأقرب للمركز بـ Ni (type 3)
for idx, atom in enumerate(atoms):
    if atom[0] == ni_id:
        atoms[idx] = (atom[0], 3, *atom[2:])
        break

with open("data.sic_ni", "w") as f:
    f.write("# SiC-3C 4x4x4 with Ni at center\n\n")
    f.write(f"{len(atoms)} atoms\n3 atom types\n\n")
    f.write(f"0.0 {L:.6f} xlo xhi\n0.0 {L:.6f} ylo yhi\n0.0 {L:.6f} zlo zhi\n\n")
    f.write("Masses\n\n1 28.0855\n2 12.0107\n3 58.6934\n\nAtoms # atomic\n\n")
    for atom in atoms:
        f.write(f"{atom[0]} {atom[1]} {atom[2]:.6f} {atom[3]:.6f} {atom[4]:.6f}\n")
print(f"data.sic_ni created – Ni at atom #{ni_id}")

# === 6. إنشاء ملف in.thz_hammer ===
print("\n[6/9] Creating in.thz_hammer (HAMMER RUN v55 - Amp=0.01)...")
AMP_TEST = 0.01
L_val = L
dt = 0.00025  # ps (مطابق للـ timestep في السكربت)

lammps_script = f"""
units           metal
atom_style      atomic
boundary        p p p
read_data       data.sic_ni
group           gNi type 3

fix             pinNi gNi spring/self 100.0

pair_style      hybrid/overlay tersoff lj/cut 8.0
pair_coeff      * * tersoff SiC.tersoff Si C NULL
pair_coeff      1 3 lj/cut 0.01 2.75
pair_coeff      2 3 lj/cut 0.02 2.40
pair_coeff      3 3 lj/cut 0.01 2.75

neighbor        2.0 bin
neigh_modify    delay 10 check yes

minimize        1.0e-12 1.0e-12 10000 100000
velocity        all zero linear

reset_timestep  0
timestep        {dt}

variable        freq equal 7.5
variable        amp equal {AMP_TEST}
variable        pulse_z equal "v_amp*sin(2*PI*v_freq*time)"
variable        L equal {L_val}
variable        L_minus_2 equal ${{L}}-2.0
region          bottom block INF INF INF INF 0 2.0
region          top block INF INF INF INF ${{L_minus_2}} ${{L}}
region          box_boundary union 2 bottom top

group           boundary_atoms region box_boundary
group           inner_atoms subtract all boundary_atoms

fix             hammer inner_atoms addforce 0.0 0.0 v_pulse_z
fix             1 all nvt temp 50 300 0.1

compute         disp gNi displace/atom
compute         zdisp gNi reduce ave c_disp[3]

variable        istep equal step
variable        zval equal c_zdisp

thermo          200
thermo_style    custom step temp pe ke etotal press v_zval

fix             print_z all print 100 "${{istep}} ${{zval}}" file disp_ni.dat title "# Step Z_disp_Ni" screen no
dump            1 all custom 500 dump.soliton.lammpstrj id type x y z
dump_modify     1 sort id

run             20000

unfix           1
unfix           print_z




"""

with open("in.thz_hammer", "w") as f:
    f.write(lammps_script)
print(f"in.thz_hammer created (HAMMER Amp={AMP_TEST})")

# === 7. تشغيل LAMMPS ===
print("\n[7/9] Running LAMMPS...")
for f in ["disp_ni.dat","dump.soliton.lammpstrj","log.lammps","ni_displacement_plot.png"]:
    if os.path.exists(f):
        os.remove(f)
        print(f"   ... cleaned {f}")

# تحسينات OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

required_files = ["data.sic_ni", "SiC.tersoff", "in.thz_hammer"]
for f in required_files:
    if not os.path.exists(f):
        print(f"Missing required file: {f}")
        exit()

try:
    result = subprocess.run(
        lammps_cmd + ["-in", "in.thz_hammer"],
        capture_output=True, text=True, timeout=1800
    )
    with open("log.lammps", "w") as f:
        f.write(result.stdout + "\n--- STDERR ---\n" + result.stderr)
    if result.returncode == 0:
        print("Simulation completed successfully!")
    else:
        print("Simulation failed. Check log.lammps.")
except subprocess.TimeoutExpired:
    print("ERROR: Simulation timeout (30 minutes).")
    exit()

# === 8. تحليل النتائج ===
print("\n[8/9] Analyzing results...")
if os.path.exists("disp_ni.dat") and os.path.getsize("disp_ni.dat") > 20:
    data = np.loadtxt("disp_ni.dat", comments="#")
    t = data[:, 0] * dt  # ps
    z = data[:, 1]
    plt.figure(figsize=(10,5))
    plt.plot(t, z, "m-", lw=2)
    plt.xlabel("Time (ps)")
    plt.ylabel("Z-displacement (Å)")
    plt.title(f"Phonon Hammer v55 (Amp={AMP_TEST})")
    plt.grid(True, ls="--", alpha=0.6)
    plt.savefig("ni_displacement_plot.png")
    plt.show()
    print(f"Max displacement: {np.max(np.abs(z)):.4f} Å")
else:
    print("disp_ni.dat missing or empty. Check log.lammps.")

# === 9. تحميل الملفات ===
print("\n[9/9] Preparing files for download...")

output_files = [
    "disp_ni.dat",
    "dump.soliton.lammpstrj",
    "log.lammps",
    "ni_displacement_plot.png",
    "thermo.log",
    "out_energy_profile.txt",
]

# for f in output_files:
#     if os.path.exists(f):
#         os.remove(f)
#         print(f"   ... cleaned {f}")

print("\n" + "="*80)
print(f"HAMMER SCRIPT (v55, Amp={AMP_TEST}) COMPLETE!")
print("="*80)
