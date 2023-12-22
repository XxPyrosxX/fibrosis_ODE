import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import NetfluxODE
import NetfluxODE_params
from collections import OrderedDict

import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as messagebox


# This method runs the simulation where we get the t values and the results
def run_simulation(tspan, y0, params, speciesNames, specific_time):
    t = []
    dt = tspan[1] / ((specific_time * 150) / 10)
    r = ode(NetfluxODE.ODEfunc).set_integrator('vode', method='adams', order=10, rtol=0, atol=1e-6, with_jacobian=False)
    r.set_initial_value(y0, tspan[0]).set_f_params(*params)
    results = np.empty([0, len(speciesNames)])
    while r.successful() and r.t <= tspan[1]:
        r.integrate(r.t + dt)
        results = np.append(results, [r.y], axis=0)
        t.append(r.t)
    return t, results


# This method is in charge of plotting the graph when given the t and results from the previous method
def plot_results(t, results, speciesNames, title, certain_nodes=False, node_entries=[]):
    if certain_nodes:
        # Use node_entries directly if certain_nodes is True
        speciesNames = [node for node in node_entries]

        # Get the indices of nodes in node_entries within speciesNames
        node_indices = [speciesNames.index(node) for node in node_entries]

        # Filter the results array to include only the specified nodes
        results = results[:, node_indices]

    sorting_indices = np.argsort(results[-1, :])[::-1]
    speciesNames_sorted = [speciesNames[i] for i in sorting_indices]

    fig, ax = plt.subplots()
    ax.plot(t, results)
    ax.set(xlabel='Time', ylabel='Fractional activation')
    ax.legend(speciesNames_sorted)
    ax.set_title(title)


# Creates a dictionary with the activity levels of each of the species
def create_results_dict(results_at_specific_time, speciesNames):
    results_dict = {}
    for i in range(len(speciesNames)):
        node = speciesNames[i]
        value = results_at_specific_time[i]
        results_dict[node] = value
    return results_dict


# Generates a bar graph comparing the control to the knocked down value
def generate_bar_graph(regular_results_dict, different_dict, knocked_down_val,
                       context_title, sorted_species=[]):
    # Extract species names and values
    species_names = list(regular_results_dict.keys())
    difference_dict = {key: different_dict[key] - regular_results_dict.get(key, 0) for key in different_dict}
    difference_values = list(difference_dict.values())

    sorted_species_names = []
    sorted_difference_values = []
    if len(sorted_species) != 0:
        ordered_dict = {key: difference_dict[key] for key in sorted_species}
        sorted_difference_values = list(ordered_dict.values())
        sorted_species_names = sorted_species
    else:
        # Combine the two lists into tuples using zip
        combined_data = list(zip(difference_values, species_names))

        # Sort the combined data based on the values in the first element of each tuple (the integers)
        sorted_data = sorted(combined_data, key=lambda x: x[0])

        # Unpack the sorted data into separate lists
        sorted_difference_values, sorted_species_names = zip(*sorted_data)

    delta_symbol = "\u0394"  # Unicode character for delta symbol
    # Create a bar graph
    fig, ax = plt.subplots()
    ax.bar(sorted_species_names, sorted_difference_values, color='blue', alpha=0.7)
    ax.set(ylabel=f"{delta_symbol} Activity (KD - Control)",
           title=f'Knocked Down Node: {knocked_down_val}; {context_title}')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    return sorted_species_names


# This runs a specific time in the simulation to gather results and returns the dictionary
def run_simulation_and_return_dict(tspan, y0, params, specific_time, speciesNames):
    t, results = run_simulation(tspan, y0, params, speciesNames, specific_time)
    index = int((specific_time - tspan[0]) / (tspan[1] / ((specific_time * 150) / 10)))-1
    results_at_specific_time = results[index, :]
    results_dict = create_results_dict(results_at_specific_time, speciesNames)
    return results_dict

# Runs all the graph and deals with all the dictionaries and calculating the change in activity
def run_graphs():
    knocked_down_val = knocked_down_val_entry.get()
    show_extra_graphs = extra_graphs_var.get()
    specific_time = float(time_entry.get())
    display_type = graph_type_entry.get()
    node = node_entry.get()
    if(display_type != "high" and display_type != "reg" and display_type != "kd" and display_type != "kdhigh"):
        display_type = ""
    node_entries = []
    incorrect_entries = []
    node_entries.append(node_entry1.get().strip())
    node_entries.append(node_entry2.get().strip())
    node_entries.append(node_entry3.get().strip())
    node_entries.append(node_entry4.get().strip())
    node_entries.append(node_entry5.get().strip())
    new_node_entries = list(filter(lambda x: x != '', node_entries))
    new_node_entries = list(OrderedDict.fromkeys(new_node_entries))
    node_entries_original = new_node_entries.copy()

    #####################################################
    # REGULAR IMPLEMENTATION (NO VALUES ARE KNOCKED DOWN)
    #####################################################
    [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.loadRegularParams()
    regular_results_dict = run_simulation_and_return_dict([0, specific_time], y0, [tau, ymax, w, n, EC50],
                                                          specific_time, speciesNames)

    if node not in speciesNames:
        node = "ROS"

    for entry in node_entries_original:
        if entry not in speciesNames:
            incorrect_entries.append(entry)
            new_node_entries.remove(entry)

    if incorrect_entries:
        # Display a message box with the list of incorrect entries
        error_message = f"Unfortunately, we weren't able to display the following nodes due to the fact they don't exist:\n{', '.join(incorrect_entries)}\nMake sure there are no typos."
        messagebox.showinfo("Incorrect Node Values", error_message)

    ##############################################################
    # KNOCKED DOWN IMPLEMENTATION (A SINGLE VALUE IS KNOCKED DOWN)
    ##############################################################
    [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.loadKnockedDownParams(knocked_down_val)
    knocked_down_results_dict = run_simulation_and_return_dict([0, specific_time], y0, [tau, ymax, w, n, EC50],
                                                              specific_time, speciesNames)

    ##############################################################
    # KNOCKED DOWN VALUE (A SINGLE VALUE IS KNOCKED DOWN) AND HIGH TGF-B (90% FOR THE WEIGHT PARAMETER FOR THE TGF-B REACTION)
    ##############################################################
    # WHEN COMPARING TO REGULAR_RESULTS IT RESULTS IN POSITIVE VALUES WHEN MOST SHOULD BE NEGATIVE
    # THE VALUES SHOULD BE WAY SMALLER AS A RESULT OF HIGH TGFB CONCENTRATION INSTEAD OF WAY HIGHER AND THI
    # THIS INCLUDES KNOCKING DOWN THE ROS NODE FURTHER THAN THE ORIGINAL -0.14 AND SHOULD BE AROUND -1 INSTEAD
    [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.bothOperationsParams(knocked_down_val)
    both_operation_results_dict = run_simulation_and_return_dict([0, specific_time], y0, [tau, ymax, w, n, EC50],
                                                                  specific_time, speciesNames)

    ##############################################################
    # HIGH TGF-B IMPLEMENTATION (90%)
    ##############################################################
    [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.hightgfbParams()
    high_tgfb_dict = run_simulation_and_return_dict([0, specific_time], y0, [tau, ymax, w, n, EC50],
                                                               specific_time, speciesNames)

    ##############################################################
    # Bar Graph Showing the Change in Activity (Knocked Down Val - Control), Signaling Context: Baseline
    ##############################################################
    sorted_species = generate_bar_graph(regular_results_dict, knocked_down_results_dict, knocked_down_val,
                                        'Signaling Context: Baseline')

    ##############################################################
    # Bar Graph Showing the Change in Activity (Knocked Down Val - Control), Signaling Context: High TGF-B
    ##############################################################
    generate_bar_graph(high_tgfb_dict, both_operation_results_dict, knocked_down_val,
                       'Signaling Context: High TGF-B', sorted_species)

    if len(new_node_entries) > 0:
        graph_text = ""
        if display_type == "kd":
            [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.loadKnockedDownParams(node)
            graph_text = "Knocked Down " + node
        elif display_type == "kdhigh":
            [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.bothOperationsParams(node)
            graph_text = "Knocked Down " + node + " and High TGF-B"
        elif display_type == "high":
            [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.hightgfbParams()
            graph_text = "High TGF-B"
        else:
            [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.loadRegularParams()
            graph_text = "Regular Parameters"
        tspan = [0, specific_time]
        t, results = run_simulation(tspan, y0, [tau, ymax, w, n, EC50], speciesNames, specific_time)
        plot_results(t, results, speciesNames, graph_text, True, new_node_entries)

    if show_extra_graphs:
        [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.loadRegularParams()
        tspan = [0, specific_time]
        t, results = run_simulation(tspan, y0, [tau, ymax, w, n, EC50], speciesNames, specific_time)
        plot_results(t, results, speciesNames, 'Regular Unchanged Parameters', )

        [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.loadKnockedDownParams(knocked_down_val)
        tspan = [0, specific_time]
        t, results = run_simulation(tspan, y0, [tau, ymax, w, n, EC50], speciesNames, specific_time)
        plot_results(t, results, speciesNames, f'Knocked Down {knocked_down_val}')

        [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.bothOperationsParams(knocked_down_val)
        tspan = [0, specific_time]
        t, results = run_simulation(tspan, y0, [tau, ymax, w, n, EC50], speciesNames, specific_time)
        plot_results(t, results, speciesNames, f'High TGF-B + Knocked Down ' + knocked_down_val)

        [speciesNames, tau, ymax, y0, w, n, EC50] = NetfluxODE_params.hightgfbParams()
        tspan = [0, specific_time]
        t, results = run_simulation(tspan, y0, [tau, ymax, w, n, EC50], speciesNames, specific_time)
        plot_results(t, results, speciesNames, f'High TGF-B ' + knocked_down_val)

    plt.show()


# Create the main window
window = tk.Tk()
window.title("Graph Generator")

# Create labels
knocked_down_val_label = ttk.Label(window, text="Knocked Down Value:")
knocked_down_val_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")

time_label = ttk.Label(window, text="Specific Time:")
time_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")

# Create entry widgets
knocked_down_val_entry = ttk.Entry(window, width=20)
knocked_down_val_entry.grid(row=0, column=1, padx=10, pady=5)

time_entry = ttk.Entry(window, width=20)
time_entry.grid(row=1, column=1, padx=10, pady=5)

# Create a checkbutton for extra graphs
extra_graphs_var = tk.BooleanVar()
extra_graphs_checkbutton = ttk.Checkbutton(window, text="Show Extra Graphs", variable=extra_graphs_var)
extra_graphs_checkbutton.grid(row=2, column=0, columnspan=2, pady=5)

# Create entry widgets for selected nodes
node_label1 = ttk.Label(window, text="Node 1:     ")
node_label1.grid(row=3, column=0, padx=10, pady=5, sticky="e")
node_entry1 = ttk.Entry(window, width=20)
node_entry1.grid(row=3, column=1, padx=10, pady=5)

node_label2 = ttk.Label(window, text="Node 2:     ")
node_label2.grid(row=4, column=0, padx=10, pady=5, sticky="e")
node_entry2 = ttk.Entry(window, width=20)
node_entry2.grid(row=4, column=1, padx=10, pady=5)

node_label3 = ttk.Label(window, text="Node 3:     ")
node_label3.grid(row=5, column=0, padx=10, pady=5, sticky="e")
node_entry3 = ttk.Entry(window, width=20)
node_entry3.grid(row=5, column=1, padx=10, pady=5)

node_label4 = ttk.Label(window, text="Node 4:     ")
node_label4.grid(row=6, column=0, padx=10, pady=5, sticky="e")
node_entry4 = ttk.Entry(window, width=20)
node_entry4.grid(row=6, column=1, padx=10, pady=5)

node_label5 = ttk.Label(window, text="Node 5:     ")
node_label5.grid(row=7, column=0, padx=10, pady=5, sticky="e")
node_entry5 = ttk.Entry(window, width=20)
node_entry5.grid(row=7, column=1, padx=10, pady=5)

graph_type = ttk.Label(window, text="Graph? (high, reg, kd, kdhigh)")
graph_type.grid(row=8, column=0, padx=10, pady=5, sticky="e")
graph_type_entry = ttk.Entry(window, width=20)
graph_type_entry.grid(row=8, column=1, padx=10, pady=5)

node = ttk.Label(window, text="Which Node for Graph?")
node.grid(row=9, column=0, padx=10, pady=5, sticky="e")
node_entry = ttk.Entry(window, width=20)
node_entry.grid(row=9, column=1, padx=10, pady=5)


# Create a button to run the script
run_button = ttk.Button(window, text="Run Graphs", command=run_graphs)
run_button.grid(row=10, column=0, columnspan=2, pady=10)

# Run the main loop
window.mainloop()
