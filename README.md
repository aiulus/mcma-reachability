# Comparing Direct and Indirect Reachability Analysis:
## A Review on Conformance-Based and Model-Agnostic Methods
## Overview
TODO
## Repository Structure

```plaintext
.
├── conformance-based/
├── data-driven/             # Experiments on parameter tuning
│   ├── dd_pc/               # Predictive control
│   ├── dd_ra/               # Reachability analysis
├── experiments/
│   ├── scripts/             # Standalone executables that produce comparative analyses
│   │   ├── compLinearDT.m   # Linear discrete-time systems
│   │   ├── compPolyDT.m     # Polynomial discrete-time systems
│   │   ├── compLipDT.m      # Lipschitz-continuous discrete-time systems
├── outputs/                 # Destination for objects meant for long-term storage
├── workspaces/              # Configuration files
└── README.md                # Repository documentation

```
## Usage
TODO
### Installation
TODO
#### Dependencies
* [MPT](https://www.mpt3.org)
* [mosek](https://www.mosek.com/products/academic-licenses/)
