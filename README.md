# Wader Scenarios Repository

<img src="Plots/MarkdownPics/lolandwaderscenarios logo.png" id="id" class="class" width=25% height=25% > 

Currently just using this repository to host the code and main plots for this project. The outputs from the code and the raw data are not stored on here was the file sizes are too large.


## Folder Structure

### `Code` folder
- `Guideline Creation`: code to turn the stakeholder rules into graded rasters that can be added together. There is a separate script for each of the four landscapes. 
- `Optimisation`: **WIP** contains 1 script that loops over each field and upgrades it to the quality of a reserve.
- `Paper plots`: contains a single script that creates plots that are needed for the main manuscript but are not needed elsewhere
- `Scenarios`: contains a 6 script workflow to run wader scenarios. It build the canvas and spatially defines targeting strategies then runs the scenarios for each region in `5- Scenario Creation.R`
- `Wader Abundance`: contains a 5 script workflow that cleans the BWWM survey data, adds attribute data to each land parcel and then builds species-specific random forest models. 
- `Xtra Scripts`: contains a series of scripts that are mainly used as references. None of them are used in the main analysis. 
