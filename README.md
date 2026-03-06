# NECOFS_FVCOM_Daily_Data_Request
Automated workflow to request forecasted bottom water temperature data from the NECOFS/FVCOM daily forecast models

This repository hosts an R script which requests forecasted bottom water temperature data from the NECOFS/FVCOM daily forecast models, specifically:
NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST
NECOFS_FVCOM_OCEAN_NORTHEAST_FORECAST

Data is requested via OpenDAP.
The R script is run via a github workflow (.yml file).
The github workflow is automatically triggered by a job hosted on cron-jobs.org.
