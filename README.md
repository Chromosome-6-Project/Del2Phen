# Chromosome6

This is a work-in-progress program written in both Python and R for the purposes of predicting phenotypes as part of the Chromosome 6 Project.

It currently consists of 2 parts:

- chr6.py: This is the main work horse of the project.
  - Compiles local Chromosome 6 Project data
  - Organizes data in an object-oriented way
  - Performs patient-patient comparisons
  - Performs phenotype predictions
  - Exports reports
  - Provides several graphing functionalities

- Chr6Network: This is an RShiny app that serves as a prototype for an interactive website for graphically displaying the Chromosome 6 Project patient genetic network
  - Uses local data exported from chr6.py
  - Plots an interactive network connecting patients with genetic overlaps
  - Allows manipulation and filtering
  - Displays summary data for patients and genes of interest
