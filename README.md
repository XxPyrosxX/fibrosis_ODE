# Fibrosis ODE Documentation:

## Introduction
The Fibrosis ODE code simulates a model of fibrosis, specifically targeting myocardial infarction scenarios.    
It employs ordinary differential equations (ODEs) to describe the dynamics of various species involved in fibrosis, such as reactive oxygen species (ROS) and transforming growth factor-beta (TGF-B).  

## Code Structure:

### The code is organized into several functions and sections:
    run_simulation Function:
        Responsible for running the ODE simulation using the ode solver from SciPy.
        Takes initial conditions, parameters, and simulation time span as inputs.
        Iterates through the simulation time, integrating ODEs at each step.

    plot_results Function:
        Plots the simulation results using Matplotlib.
        Allows for plotting specific nodes and sorting them based on their final values.

    create_results_dict Function:
        Creates a dictionary containing the activity levels of each species at a specific time.

    generate_bar_graph Function:
        Generates a bar graph comparing control and knocked-down values for specified species.
        Provides options to sort and customize the graph.

    run_simulation_and_return_dict Function:
        Runs a specific time in the simulation and returns the results as a dictionary.

    run_graphs Function:
        Utilizes the above functions to run various simulations and generate graphs.
        Takes input parameters from a Tkinter GUI.

## Tkinter GUI:
The Tkinter window serves as a user interface for setting simulation parameters and controlling graph generation. Key features include:

    Entry fields for specifying the knocked-down value, specific time, and node for the graph.
    Checkbutton to show/hide extra graphs for all species.
    Entry fields for up to five nodes for their individual graph.
    Options for selecting the type of graph (high, reg, kd, kdhigh).

## Data Handling:
The regular ODE model manually specifies parameters like tau, ymax, w, n, and EC50 in the NetfluxODE_params module.    
In contrast, the rapid fibrosis ODE model reads parameters from an Excel file (data.xlsx), allowing for easy modification and experimentation without modifying the code.   
 
## Integration with Agent-Based Model:
There is ongoing work to integrate an agent-based model with the existing ODE model.    
The combined model aims to simulate fibroblast activity during myocardial infarction.   
This integration would provide a more comprehensive understanding of fibrosis dynamics by considering both cellular-level interactions (agent-based) and biochemical processes (ODE).   

## Conclusion:
The Fibrosis ODE code, coupled with the Tkinter window, provides a flexible and user-friendly tool for simulating and visualizing fibrosis scenarios.     
The ongoing integration with an agent-based model demonstrates a multidimensional approach to studying complex biological processes.   
This combined model holds the potential to offer more accurate and insightful simulations of fibroblast activity during myocardial infarction.   
    
