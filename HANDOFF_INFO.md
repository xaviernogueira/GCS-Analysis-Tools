# Hand-off info

**Date:** January 5th 2025.
**Scope:** Information that will assist a developer in managing this repository, and moving forward with improvements.

## Setup

In a video I made previously I demonstrate set-up from an IDE (i.e., PyCharm). The easier way to get started is on the Windows 
Powershell command line. Note that using Windows is a must due to the `ArcGIS Pro` dependency.
1. **Setup Powershell w/ conda:** Open powershell and install `wget` with `winget` (built-in), then refresh the shell, then [install miniconda](https://docs.anaconda.com/miniconda/install/) using the command shown in the linked docs,
and following the installation Wizard. Finally reset the shell once more. This can be done as follows:
```bash
winget install wget
powershell

wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" -outfile "./Downloads/Miniconda3-latest-Windows-x86_64.exe"
powershell

# verify conda access
conda --version
```

2. **Create a clone of the default ArcGIS Pro conda env:** To use `arcpy` in the code, we need to be using an ArcPro conda 
environment. Using ArcPro UI (see [here](https://pro.arcgis.com/en/pro-app/latest/arcpy/get-started/clone-an-environment.htm)), one can identify this environment, and "clone" an exact copy to a desired location. We need to 
do this, since we will be installing additional dependencies, which is not allowed in the built-in conda environment. The name you
give it does not matter, nor should the location.
3. **Activate the clone environment:** Next we will want to make sure we can activate the clone ArcPro environment. In Powershell
run the following:
```bash
conda env list
```
This should show the cloned environment with whatever name you give it. Let's say it's `arcpro-env-clone`, we then activate:
```bash
conda activate arcpro-env-clone
```
4. **Install additional dependencies:** Next with the cloned environment active, install the additional software dependencies using conda:
```bash
conda install -c conda-forge scipy
conda install -c conda-forge pillow
conda install -c conda-forge plotly
conda install -c conda-forge seaborn
conda install -c conda-forge openpyxl
conda install -c conda-forge pytest # optional, to run tests

# or in one command
conda install -c conda-forge scipy pillow plotly seaborn openpyxl pytest
```
5. **Run the GUI!:** With the environment still active, navigate to wherever you cloned `GCS-Analysis-Tools` repo. Once there 
you can fire up the GUI with simply:
```bash
python gui_tkinter.py
```

## Repo Terminology
To make the rest of this guide and future communications clean, I will define 
I mean by certain terms.
* "GUI functions": the Python function that takes arguments directly from user input. 
These functions are 1-to-1 with what each "Run" button in the GUI will do.
* "Process functions": this will refer to the underlying geospatial functions that 
are triggered by the GUI functions. For example, making a centerline is a process function, 
as is creating "station lines" from said centerline. However, both may be run sequentially by a GUI function.


## Tests

To make sure the handoff is successfull, and to give you an easier time avoiding breaking anything, 
I decided to write a few tests for the GUI functions specifically. While tests 
for the process functions would be cool to assure behavior in edge cases (i.e., traditional "unit tests"),
that would have been overly time consuming, and considering code will change, it made sense to verify 
the behavior from a user input -> user output way (i.e., GUI functions only).

