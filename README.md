## Geomorphic Covariance Structure (GCS) analysis GUI
The `gcs_gui.py` script in this package when ran will pull up a graphical-user interface / GUI that simplifies the
process of multi-flow stage GCS analysis. 

For more information on GCS analysis please reference my PI Gregory Pasternack's [video](https://www.youtube.com/watch?v=VSMK72FbTfI) which is part of a series introducing the topic. To see real-world outputs of this software within an appropriate scientific context please reference [this recording](https://www.youtube.com/watch?v=fp1Nag4kN7s) of my graduate/MS exit seminar.

### Package contents
- **Graphic User Interface (GUI) executable script (`gcs_gui.py`)**
- Python files storing the underlying methods
- Copy of RapidLasso's [LAStools software](http://lastools.org/) `.exe` command line executables (`LAStools/`).
- Copy of ['Breeze' GUI theme](https://github.com/MaxPerl/ttk-Breeze) (`ttk-Breeze-master/`).

### Pre-requisite packages
- `arcpy` (+ Spatial Analyst licence) -- Comes with the ArcPro Python 3 default environment. 
  Clone environment to enable external package installation. 
- [`pillow`](https://python-pillow.org/) to power pop up image outputs.
- [`seaborn`](https://seaborn.pydata.org/) for publication-ready static plotting.
- [`plotly`](https://plotly.com/) for [Sankey diagram](https://plotly.com/python/sankey-diagram/) plotting.

### **Please read our [documentation](https://gcs-gui-documentation.readthedocs.io) before attempting to use this GUI!**
