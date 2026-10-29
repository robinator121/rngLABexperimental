# Codebase Metrics (cloc equivalent)

**Project:** RNGLab Experimental Platform  
**Date:** October 2026  
**Method:** Exact syntactic line-count audit excluding `.git`, `.venv`, and `__pycache__`.

---

## Summary by Language

| Language / File Type | Files | Blank Lines | Comment Lines | Code Lines | Total Lines |
| :--- | ---:| ---:| ---:| ---:| ---:|
| **Python 3** (`.py`) | 219 | 3,237 | 555 | 24,171 | 27,963 |
| **Markdown Documentation** (`.md`) | 43 | 593 | 209 | 2,156 | 2,958 |
| **JSON Manifests / Configurations** (`.json`) | 5 | 0 | 0 | 964 | 964 |
| **Arduino C++** (`.ino`) | 1 | 54 | 9 | 166 | 229 |
| **Text Specifications** (`.txt`) | 2 | 14 | 0 | 44 | 58 |
| **TOML Configuration** (`.toml`) | 1 | 3 | 0 | 15 | 18 |
| **Dependency Lockfile** (`.lock`) | 1 | 0 | 3 | 7 | 10 |
| **Total** | **272** | **3,901** | **776** | **27,523** | **32,200** |

---

## Component Breakdown

- **Core Scientific Libraries (`src/rnglab*`):** Packages `rnglab11` through `rnglab17` implementing source abstractions, statistical batteries, VMC, PIMC, exact Ising references, TRIDENT kernels, and CALIBRATE failure predictors.
- **Financial Risk Modeling (`src/rnglab_finance`):** European option pricing, portfolio Value at Risk (VaR), Expected Shortfall (ES), and stochastic fault injection tests under model risk principles.
- **Physical Acquisition & Firmware (`src/rnglab15`, `hardware/`):** Serial acquisition daemons, NIST SP 800-22 runner, and Arduino sensor firmware.
- **Automated Verification Suites (`tests/`, `scripts/`):** 237 passing unit, mathematical, and regression tests.
