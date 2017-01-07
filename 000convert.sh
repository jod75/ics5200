#!/usr/bin/env bash
jupyter nbconvert Lab/000\ -\ ChEMBL\ Helper.ipynb --to python
jupyter nbconvert Lab/001\ -\ Python\ Helper.ipynb --to python
jupyter nbconvert Lab/002\ -\ Elements.ipynb --to python
jupyter nbconvert Lab/003\ -\ HDFS\ Helper.ipynb --to python
echo renaming \(000\ -\ ChEMBL\ Helper.py\) --\> chemblhelper.py
mv Lab/000\ -\ ChEMBL\ Helper.py Lab/chemblhelper.py
echo renaming \(001\ -\ Python\ Helper.py\) --\> pythonhelper.py
mv Lab/001\ -\ Python\ Helper.py Lab/pythonhelper.py
echo renaming \(002\ -\ Elements.py\) --\> elements.py
mv Lab/002\ -\ Elements.py Lab/elements.py
echo renaming \(003\ -\ HDFS\ Helper.py\) --\> hdfshelper.py
mv Lab/003\ -\ HDFS\ Helper.py Lab/hdfshelper.py
