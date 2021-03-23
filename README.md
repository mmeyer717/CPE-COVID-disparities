# Center for Policing Equity COVID-19 Modeling Project

We sought to simulate COVID-19 spread in a large, synthetic U.S. city, focusing on contributions of three domains in driving spread of COVID-19 and related racial disparities:
a)	low-wage, essential work, 
b)	police-public contact, and 
c)	jail and prison churn.

We used a Susceptible-Infected-Recovered (SIR) model to illustrate the high level dynamics of the disease spread (rather than projecting precise case numbers in an actual city).

Please refer to the Center for Policing Equity COVID-19 Modeling Project documentation for further details about our modeling choices. 

## Directories

* **model**: Run our model in its entirety using Model.ipynb. Necessary functions for running the model can be found in python files also in this directory.
* **input**: Find the contact matrices and initial subpopulation group sizes as spreadsheets; these are necessary inputs into our model.

## Running The Code
The code can notebook Model.ipynb can be run through jupyter notebooks. It is written in python3 and depends only on pandas, numpy and matplotilib.

