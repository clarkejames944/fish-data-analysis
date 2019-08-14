# fish-data-analysis
Analysis of Tiger Flathead otolith growth data

This project details the production of an Integral Projection Model (IPM) to describe the otolith size dynamics of fish.
The project began with exploratory analysis on the otolith dataset for the Tiger Flathead from fishing activities off the south
east coast of Australia.

I will detail what each of the folders and files are for.

# Analysis folder

## compare.R 
- this r script details the comparison of the Bayesian models once they have been fitted using Stan.

## explore.R 
-this indicates the initial exploration of the dataset including the fitting of multiple GAMs to understand factors
which impact the otolith growth

## fit_hinge_age_effect_all_fixed_effects.R

Thes files beginning with fit are all used to run the corresponding stan models

## fit_hinge_basemodel.R

## fit_hinge_intercept+pre_thres_slope_age_effect.R

## fit_hinge_intercept+threshold_age_effect.R

## fit_hinge_intercept_age_effect.R

## fit_hinge_sex_diff.R

## fit_hinge_temp_effect.R

## fit_hinge_threshold+pre_thres_slope_age_effect.R

## fit_hinge_threshold_age_effect.R

## fit_hinge_zone_effect.R

## fit_hinge_zone_temp_sex_diff.R

## fit_intercept.R

## fit_intercept_and_threshold.R

## setup.R

Setup script for the packages we want to use in our analyses and the manipulation of the model
that we want.

## summarise.R

Used to summarise the outputs from the regression models, e.g. to create and check the pairs plot produced by the 
Stan model.
###################################################

# Figures folder 

A collection of figures produced for the project.

Including a Map.R file which was used to build the map file

plots.R was a script used for plotting some graphs to describe the dataset.

####################################################

# IPM folder

Folder for the running of the IPM and for the creation of some graphs from this IPM.
When collecting multiple iterations I had a seperate R script

####################################################

# IPM summaries

All the data-frames produced by the IPM to be used to make the graphs that
are found in the IPM folder

####################################################

# maps

This folder was created for when I made the Australia map and I was hoping to incorporate some shape files,
but I didn't get the chance to do this.

The main map script is found in the analysis folder

####################################################

# miscellaneous 

Some of the very early analysis I carried out on the data, which wasn't used in the final analysis.

####################################################

# models

Saved the fitted Stan models here

####################################################

# notation

Some files for building the mathematical notation that I was to use in the model

####################################################

#otoliths(working)

The original datasets and some preliminary work carried out by Federico Giazza

####################################################

# par_names

I stored the parameter names from each regression model to be used for the IPM in here

This was so that I didn't have to create a separate IPM for each but instead intput this section.
There is a call to this folder in the IPM code.

####################################################

# posterior predictive checks

This is where I carried out the posterior predictive checks on the model chosen by the comparison, to check
if the model sufficiently described our data.

####################################################

# stan

This is where all the Stan files are stored for each regression model.

#####################################################

# vitals

This is similar to the 'par_names' folder in that this contains the vital rates functions for each regression model,
to be input into the IPM. This was so that the IPM could remain a general one. There is a call to this file 
in the IPM code.
