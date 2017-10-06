# zvi-metagenomics
ProjectTemplate use overview
* data that you will frequently load are kept in the `data` directory, limit use as the data are loaded into memory when the project is loaded. 

Loading the project provides access to data and functions in lib and automatically loads packages defined in the `config/config.dcf` file.  
To limit over loading your environment only use the `lib` and `data` directories for commonly used functions and data. 
The `reports` and `scratch` directories have data subdirectories for storing intermediate results that are used in other analysis but do not want to load every time you load the project. 


To load the project from the project root directory run 
```
library('ProjectTemplate')
load.project()
```