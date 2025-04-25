# 🚀 HopIt – Hopper Iterative System Analysis Tool

**HopIt** is a modular Python tool for performing an iterative system-level mass estimation of a hopper (small rocket vehicle). It evaluates propellant mass, tank sizing, propulsion, batteries, structure, and more – until mass convergence is reached.

---

## 📁 Project Structure

hopit/
├── main.py              # Entry point: runs iterative system analysis
├── system.json          # Example system configuration (optional)
└── core/
    ├── functions.py     # Utility functions for mass estimation
    └── casadi_core.py   # CasADi-based tank sizing core

---

## ▶️ Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/flxbrt/hopit.git
cd hopit
```
### 2. Install dependencies
Make sure you have Python ≥ 3.8 installed. Then install required packages:

```bash
pip install -r requirements.txt
```
### 3. Run the tool
```bash
python main.py
```

By default, it runs with a hardcoded config defined inside main.py. The tool will:

>>Load system configuration

>>Estimate subsystem masses

>>Iterate until total mass convergence

>>Print the result (if sys_print = True)

>>Optionally store output (if store_system = True)


⚙️ Configuration

```bash
config = {
    'Flight Time [s]': 80,
    'Oxidizer': 'O2',
    'Fuel': 'C2H5OH',
    'Oxidizer Pump': False,
    'Fuel Pump': False,
    'Fuel Self Pressurised': False,
    'Oxidizer Self Pressurised': False
}
```

Later versions may load this from a .json or .yaml file.

💡 Features

Iterative loop for consistent subsystem mass convergence

Tank sizing via CasADi-based nonlinear solving

Modular code structure (easy to extend with additional subsystems)

Designed for conceptual rocket vehicle studies

self pressurized vs non self pressurized

chamber sizing

pump fed or pressure fed




📌 Planned Enhancements
CLI interface with argparse

JSON-based config loading

Plotting of key subsystem masses

Web-based visualization (long-term)

Couling to THERMAT

addind performance mode --> allows to copmute system performance given a system --> helps to understand the performance of the systm for different available COTS components
