# AGECalibrator
This is the source code for the AGECalibrator app

The files: 
* `functions.R` contains the necessary functions. 
* `main.R` runs the Shiny app
* `ui.R` and `server.R` contain the code for the user interface and server, respectively.  


The `data` folder contains test data, which you can use to test your installation.
The `docs` folder contains the manual. 

## Local installation
### Dependencies
* `R` => 3.6.0
* `shiny`
* `sninydashboard`
* `DT`
* `ggplot2`
* `dplyr`
* `tidyr`
* `remotes`
* `rPeaks`
* `Rcpp`
* `coin`
* `multcomp`

### Installation

The following installation instructions are provided for the case if you (for any reason) prefer to run the app locally. 

Open the `main.R` file in RStudio and run it to install the necessary packages. Then press the Run button to run the app. 

The installation was tested in Linux-based systems only but should be platform independent (provided that you know how to install all necessary packages). 
