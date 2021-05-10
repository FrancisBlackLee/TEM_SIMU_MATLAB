# TEM_SIMU_MATLAB
Updated on 2021/5/10 by Francis Black Lee
## Introduction
TEM_SIMU_MATLAB is part of a simulation software suite for electron microscopy and spectroscopy. The repo is written in MATLAB, mainly for demonstrating the simulation principles and fast implementation of ideas. Other APIs built in this project will be included in this repo in the future.
## Installation
Download the repo to your load path, and add the folder and its subfolders to your MATLAB path, note that the GUIs provided in this repo were designed using MATLAB appdesigner, so make sure your MATLAB version supports the app designer and recompile them towards to your MATLAB version before using them.
## Samples
The repo provides severals samples (under ``\script_samples`` folder) as demonstration, they will give you a quick start when you start to write your first script. You can also refer to the ``\tests`` folder for more scripting examples.
## Use the Newest APIs
New APIs should always be labelled with ``_X`` as the suffix of their names before they are accepted as formal APIs, ``X`` denotes the experimental version. Note that even as experimental APIs, they are always tested before being commited to the master branch, so feel free to use these APIs and you are always welcome to report any problem of the codes to the author.
## Legacy Codes
Legacy codes are placed under the ``\legacycodes`` folder, any formal or experimental API should be labelled with exact version index on its name before it is put into the ``\legacycodes`` folder.
## Current Development Status
This repo is still under active development, now it is capable of simply handling crystal structures (including a CIF parser), TEM imaging/diffraction and STEM imaging/diffraction.
