# Hocking_Chi_LEMs
Inputs, driver scripts, and analyses used to create the landscape evolution models of (manuscript in progress).

Requirements:
* [Landlab](https://landlab.github.io/#/)
* [QGIS from Conda](https://plugins.qgis.org/planet/tag/conda/)
* [Topotoolbox](https://topotoolbox.wordpress.com/)
* [Topographic Analysis Toolkit](https://github.com/amforte/Topographic-Analysis-Kit)

Each folder inside Hocking_Chi_Lems_master contains the input parameters and scripts used to run and analyze a unique landscape evolution model. For each model directory, the nested sturcture is as follows:

* **Input:** A folder containing two .csv files (Spin_up.csv and Main_run.csv). These files prescribe the parameters used to drive the model spin-up and main run, respectively. We used the .csv format to make the models more accessible to undergraduate students.
* **Model:** A folder containing the python script used to drive the model. The script expects the master directory to have the same nested file structure as the Github repository. The script imports values specified in Spin_up.csv and Main_run.csv, creates an "Output" folder in the model directory, drives the model, and periodically exports data to the output folder. One of the exported products are DEMs of specified timesteps (stored in "Output\Main_Run\topographic__elevation__tif"); this directory and its contents are used by scripts in the analysis step below.
* **Analysis:** A folder containing materials used to analyze the landscape evolution models. Nested structure is as follows:
  - **Input:** A folder containing a .csv file (Parameters.csv) that prescribes parameters used to drive the Hocking_Chi_LEMs_Analysis Matlab files. We used the .csv format to make the models more accessible to undergraduate students.
  - **Scripts:** A folder containing Matlab scripts used to analyze the landscape evolution models. Scripts are as follows:
    - **Hocking_Chi_LEMs_Analysis_UL.m or Hocking_Chi_LEMs_Analysis_VL.m:** Main driver analysis script. The script expects the master directory to have the same nested file structure as the Github repository. The script imports values specified in Parameters.csv, creates an "Output" folder in the analysis directory, loops through the DEMs stored in "Output\Main_Run\topographic__elevation__tif", and exports analysis products.
    - **SegmentPicker_CS.m:** A modified version of the SegmentPicker function from the [Topographic Analysis Toolkit](https://github.com/amforte/Topographic-Analysis-Kit). The function has been modified to run faster within the Hocking_Chi_LEMs_Analysis_UL.m and Hocking_Chi_LEMs_Analysis_VL.m scripts. 
