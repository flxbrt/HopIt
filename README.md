# 🚀 HopIt – Hopper Iterative System Analysis Tool

**HopIt** is a modular Python tool for performing an iterative system-level mass estimation of a hopper (small rocket vehicle). It evaluates propellant mass, tank sizing, propulsion, batteries, structure, and more – until mass convergence is reached.

---

## 📁 Project Structure
```plaintext
hopit/
├── main.py              # Entry point: runs iterative system analysis
├── system.json          # Example system configuration (optional)
└── core/
    ├── functions.py     # Utility functions for mass estimation
    └── casadi_core.py   # CasADi-based tank sizing core
```
---

## 🛫 Getting Started

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
- Load system configuration
- Estimate subsystem masses and compute engine performance parameter
- Iterate until total mass convergence
- Print the result (if sys_print = True)
- Optionally store output (if store_system = True)


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
- Designed for conceptual Vertical take off vertical landing vehicle studies propelled by a bi-liquid rocket engine
- Modular code structure (easy to adapt individual subsystems or extend with additional subsystems)
- Tank sizing equations solved via root finding with CasADi
- Allows self pressurized and non self pressurized propellants
- Allows pressure fed and epump fed cycles
- Computes main combustion chamber performance parameters

🚀 Planned Enhancements
- Plotting of system sketch with subsystems to scale
- Adding 'Performance Mode' to compute the flight time and engine performance given a system (reverse to current mode); this is practical when one once to understand the potential performance of the system for different COTS components
- Coupling to [THERMAT](https://github.com/flxbrt/THERMAT)
