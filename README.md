# Time Series Analysis PM2.5 Delhi

This repository contains files related to the project of Time Series Analysis of pollutant PM2.5 in Delhi over the years of 2015-2021. This project explores the various Staistical Methods of Time Series Analysis fit for the given dataset.

## Files

### Data Store
- [city_hour.csv](city_hour.csv) : Contains hourly recorded concentration of pollutants in various cities over the mentioned period
- [out.csv](out.csv) : Intermediate dataset of the above one, cleaned for the use in R file
- [out2.csb](out2.csv) : Same file as above, with indices

### Code
- [1st Jupyter Notebook](Final_FTSA.ipynb): Code for initial cleaning and analysis along with some simple models
- [2nd Jupyter Notebook](Final_FTSA_2.ipynb): Code for Random Forest, XG Boost and LSTM Models alongside TabPy Integration
- [1st R File](FTSA_.R): Code for Dynamic Harmonic Regression and TBATS Model
- [2nd R File](FTSA2_.R): Code for Outlier Treatment, GAM Model

### Report
- [Word Report](Final_Masters_Project.docx): Report for the entire project
