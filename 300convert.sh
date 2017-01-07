#!/usr/bin/env bash
jupyter nbconvert Lab/300\ -\ Ligand\ framework.ipynb --to python
mv Lab/300\ -\ Ligand\ framework.py Lab/moleculehelper.py
