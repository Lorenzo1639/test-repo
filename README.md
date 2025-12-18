# Batch Topology Optimization of Airfoil Sections in Abaqus

This repository contains a Python script for Abaqus/CAE that performs
batch topology optimization of 2D airfoil sections subjected to pressure loads.
The script automatically reads airfoil geometry and pressure distributions
from CSV files, builds the finite element model, applies boundary conditions
and loads, and runs a stiffness-based topology optimization for multiple
cases and slices.

The repository also includes a Visual Studio Code task (`tasks.json`) that
allows running the script in noGUI mode directly from the editor.

The code was developed and used in the context of a semester project and is
provided to support the reproducibility of the numerical results.
 
