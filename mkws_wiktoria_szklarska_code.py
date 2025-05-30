import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product


def generate_title(mix_label, mode, plot_type):
    suffix = {
        ("HP", "isotherms"): "adiabatic temperature at constant pressure",
        ("UV", "isotherms"): "adiabatic temperature at constant volume",
        ("HP", "gradients"): "adiabatic temperature at constant pressure",
        ("UV", "gradients"): "adiabatic temperature at constant volume"
    }[(mode, plot_type)]
    
    base = {
        "H2 + O2 (poor)": "Hydrogen + oxygen, poor mixture, ",
        "H2 + O2 (rich)": "Hydrogen + oxygen, rich mixture, ",
        "H2 + Air (poor)": "Hydrogen + air, poor mixture, ",
        "H2 + Air (rich)": "Hydrogen + air, rich mixture, ",
    }[mix_label]

    return base + suffix

def get_group_title(group_label, mode, plot_type):
    prefix = {
        "O2": "Hydrogen + oxygen, ",
        "Air": "Hydrogen + air, "
    }[group_label]

    suffix = {
        ("HP", "isotherms"): "adiabatic temperature at constant pressure",
        ("UV", "isotherms"): "adiabatic temperature at constant volume",
        ("HP", "gradients"): "adiabatic temperature at constant pressure",
        ("UV", "gradients"): "adiabatic temperature at constant volume"
    }[(mode, plot_type)]

    return prefix + suffix

#INITIAL CONDITION - CANTERA USES SI UNITS
p_initial = 101325 
T_o_range = np.arange(25, 101, 5) 

#CONCENTRATION OF HYDROGEN RANGE FOR EACH MIXTURE AND LFL, UFL
H2_ranges = {
    "H2 + O2 (poor)": np.arange(0.03, 0.0601, 0.0025),
    "H2 + Air (poor)": np.arange(0.03, 0.0601, 0.0025),
    "H2 + Air (rich)": np.arange(0.75, 0.8501, 0.0025),
    "H2 + O2 (rich)": np.arange(0.85, 0.9501, 0.0025),
}

all_data = [] #array for storing data

#CALCULATION LOOPS

for mix_label, H2_fracs in H2_ranges.items():
    for H2_frac, T_C in product(H2_fracs, T_C_range):
        T_K = T_C + 273.15 #conversion to Kelvins
        
        #DEIFNING THE COMPOSITION OF MIXTURES DEPENDING ON HYDROGEN FRACTION
        if "O2" in mix_label:
            O2_frac = 1 - H2_frac
            composition = f"H2:{H2_frac}, O2:{O2_frac}"
        elif "Air" in mix_label:
            O2_frac = (1 - H2_frac) / 4.76
            N2_frac = O2_frac * 3.76
            composition = f"H2:{H2_frac}, O2:{O2_frac}, N2:{N2_frac}"

        try:
            # Create Cantera gas object
            gas = ct.Solution("h2o2.yaml")

            # Set state: constant pressure
            gas.TPX = T_K, P_initial, composition
            gas.equilibrate("HP")
            T_ad_HP = gas.T

            # Set state: constant volume
            gas.TPX = T_K, P_initial, composition
            gas.equilibrate("UV")
            T_ad_UV = gas.T

            # Append result to list
            all_data.append({
                "Mixture": mix_label,
                "H2 Concentration [%]": H2_frac * 100,
                "Initial T [C]": T_C,
                "T_ad_HP [C]": T_ad_HP - 273.15,
                "T_ad_UV [C]": T_ad_UV - 273.15,
            })

        except Exception as e:
            print(f"Error with {mix_label}, H2={H2_frac:.4f}, T={T_C}C: {e}")

#DATA VISUALIZATION 

for mix_label in df["Mixture"].unique():
    df_part = df[df["Mixture"] == mix_label]
    h2_vals = sorted(df_part["H2 Concentration [%]"].unique())
    T_vals = sorted(df_part["Initial T [C]"].unique())
    X_grid, Y_grid = np.meshgrid(h2_vals, T_vals) #return a tuple of coordinate matrices from coordinate vectors -> for isotherms
    
    #RESULTS ARRAYS - DEFINE
    Z_HP = np.empty_like(X_grid)
    Z_UV = np.empty_like(X_grid)
    
    #RESULTS ARRAYS - FILL UP 
    for i, T in enumerate(T_vals):
        for j, H2 in enumerate(h2_vals):
            row = df_part[(df_part["H2 Concentration [%]"] == H2) & (df_part["Initial T [C]"] == T)]
            Z_HP[i, j] = row["T_ad_HP [C]"].values[0]
            Z_UV[i, j] = row["T_ad_UV [C]"].values[0]
    #PLOTS!
    fig1, axs1 = plt.subplots(1, 2, figsize=(14, 6)) #1 row, 2 columns, size in inches
    c1 = axs1[0].contour(X_grid, Y_grid, Z_HP, colors='black') #contour create the isotherm maps, X-grid stores concentration of hydrogen, Y initial temperature (the pair T0-%), Z has corresponding values of adiabatic temperatures
    axs1[0].clabel(c1, inline=True, fontsize=8)
    axs1[0].set_title(get_title(mix_label, "HP", "isotherms"))
    axs1[0].set_xlabel("H2 Concentration [%]")
    axs1[0].set_ylabel("Initial Temperature [C]")

    c2 = axs1[1].contour(X_grid, Y_grid, Z_UV, colors='blue')
    axs1[1].clabel(c2, inline=True, fontsize=8)
    axs1[1].set_title(get_title(mix_label, "UV", "isotherms"))
    axs1[1].set_xlabel("H2 Concentration [%]")
    axs1[1].set_ylabel("Initial Temperature [C]")

    plt.tight_layout()
    plt.savefig(f"isotherms_{mix_label.replace(' ', '_')}.png")
    plt.show()