import pandas as pd
import numpy as np
from scipy.integrate import ode
import NetfluxODE

SPECIES_NAMES = ['AngII', 'AT1R', 'AGT', 'ACE', 'NOX', 'ROS', 'ET1', 'ETAR', 'DAG', 'PKC', 'TRPC', 'NE', 'BAR',
                    'Forskolin', 'AC', 'cAMP', 'PKA', 'CREB', 'CBP', 'TGFB', 'TGFB1R', 'smad3', 'smad7', 'latentTGFB',
                    'BAMBI', 'PDGF', 'PDGFR', 'NP', 'NPRA', 'cGMP', 'PKG', 'mechanical', 'B1int', 'Rho', 'ROCK', 'Ca',
                    'calcineurin', 'NFAT', 'IL6', 'gp130', 'STAT', 'IL1', 'IL1RI', 'TNFa', 'TNFaR', 'NFKB', 'PI3K',
                    'Akt', 'p38', 'TRAF', 'ASK1', 'MKK3', 'PP1', 'JNK', 'abl', 'Rac1', 'MEKK1', 'MKK4', 'ERK', 'Ras',
                    'Raf', 'MEK1', 'FAK', 'epac', 'Factin', 'FA', 'migration', 'cmyc', 'CTGF', 'proliferation', 'SRF',
                    'EDAFN', 'aSMA', 'AP1', 'TIMP1', 'TIMP2', 'PAI1', 'proMMP14', 'proMMP1', 'proMMP2', 'proMMP9',
                    'MMP1', 'MMP2', 'MMP9', 'MMP14', 'fibronectin', 'periostin', 'CImRNA', 'CIIImRNA', 'CI', 'CIII']

SPECIFIC_TIME = 100

# copied over from Netflux_ODE
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

# copied over from Netflux_ODE
def create_results_dict(results_at_specific_time, speciesNames):
    results_dict = {}
    for i in range(len(speciesNames)):
        node = speciesNames[i]
        value = results_at_specific_time[i]
        results_dict[node] = value
    return results_dict

def run_simulation_and_return_dict(tspan, y0, params, specific_time, speciesNames):
    t, results = run_simulation(tspan, y0, params, speciesNames, specific_time)
    index = int((specific_time - tspan[0]) / (tspan[1] / ((specific_time * 150) / 10)))-1
    results_at_specific_time = results[index, :]
    results_dict = create_results_dict(results_at_specific_time, speciesNames)
    return results_dict

# Read the spreadsheet into a pandas DataFrame
df = pd.read_excel("data.xlsx")  # Replace with the actual file path and name

# Specify the columns to be considered (excluding the first column)
columns_to_read = ['TAU', 'Y Max', 'Y Init', 'W', 'N', 'EC50']

results_list = []

# Iterate through the rows and extract parameter values
for index, row in df[columns_to_read].iterrows():
    # Check if there's data in the row
    if not pd.isnull(row['TAU']):
        # Split the string by commas and convert to a NumPy array
        tau_values = np.array([float(value) for value in row['TAU'].split(',')])
        y_max_values = np.array([float(value) for value in row['Y Max'].split(',')])
        y_init_values = np.array([float(value) for value in row['Y Init'].split(',')])
        w_values = np.array([float(value) for value in row['W'].split(',')])
        n_values = np.array([float(value) for value in row['N'].split(',')])
        ec50_values = np.array([float(value) for value in row['EC50'].split(',')])

        result_dict = run_simulation_and_return_dict([0, SPECIFIC_TIME], y_init_values,
                                       [tau_values, y_max_values, w_values, n_values, ec50_values],
                                       SPECIFIC_TIME, SPECIES_NAMES)

        # Append the result dictionary to the DataFrame
        # Append the result dictionary to the list
        results_list.append(result_dict)

# Convert the list of dictionaries to a DataFrame
results_df = pd.DataFrame(results_list)

# Save the DataFrame to an Excel file
results_df.to_excel("output_results.xlsx", index=False)

print("Results exported to output_results.xlsx")
