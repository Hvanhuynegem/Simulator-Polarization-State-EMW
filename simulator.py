# Simulator of the polarization state of an EM wave using the polarization pattern method

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Parameters given
f = 9.3678e9  # Frequency of the EMW (Hz)
e_d = 2.4799 - 0.0142j  # Complex dielectric permittivity of the reflecting plate
d = 4.0497  # Distance between the antennas (m)
h = np.array([0.2903, 0.4674, 0.6445])  # Vector of height differences, simplified to the first element for clarity


# Polarization phasor function
def polarization_phasor(phi, tau):
    return (np.tan(phi) + 1j * np.tan(tau)) / (1 - np.tan(phi) * 1j * np.tan(tau))

# Fresnel coefficients function
def fresnel_coefficients(theta_i, e1, e2):
    horizontal = (np.cos(theta_i) - np.sqrt((e2/e1) - np.sin(theta_i)**2)) / (np.cos(theta_i) + np.sqrt((e2/e1) - np.sin(theta_i)**2))
    vertical = -((e2 * np.cos(theta_i) - e1 * np.sqrt((e2/e1) - np.sin(theta_i)**2)) / (e2 * np.cos(theta_i) + e1 * np.sqrt((e2/e1) - np.sin(theta_i)**2)))
    return horizontal, vertical

# Incident angle function
def calculate_incident_angle(d, h):
    return (np.pi / 2) - np.arctan(d / (2 * h))


def polarization_ratio_antenna_metallic(f, e_d, d, h, phi, tau):
    p = polarization_phasor(phi, tau)
    horizontal, vertical = 1, -1
    omega = 2 * np.pi * f
    c = 3e8
    R_L = 2 * np.sqrt(h**2 + (0.25 * d**2))
    k = omega / c
    dR = d - R_L
    p_r = p * (1 + vertical * np.exp(-1j * k * dR)) / (1 + horizontal * np.exp(-1j * k * dR))
    return p_r

def polarization_ratio_antenna_dielectric(f, e_d, d, h, phi, tau):
    p = polarization_phasor(phi, tau)
    theta_i = calculate_incident_angle(d, h)
    horizontal, vertical = fresnel_coefficients(theta_i, 1, e_d)
    omega = 2 * np.pi * f
    c = 3e8
    R_L = 2 * np.sqrt(h**2 + (0.25 * d**2))
    k = omega / c
    dR = d - R_L
    p_r = p * (1 + vertical * np.exp(-1j * k * dR)) / (1 + horizontal * np.exp(-1j * k * dR))
    return p_r

def length_between_antennas(h, angle):
    # angle = (np.pi / 2) - np.arctan(d / (2 * h))
    return np.tan(np.pi / 2 - angle) * 2 * h

# Function to calculate polarization state of the received wave
def polarization_ratio_antenna_brewster(f, e_d, d, h, phi, tau):
    p = polarization_phasor(phi, tau)
    brewster_angle = np.arctan(np.sqrt(e_d.real))
    horizontal, vertical = fresnel_coefficients(brewster_angle, 1, e_d)
    omega = 2 * np.pi * f
    c = 3e8
    d = length_between_antennas(h, brewster_angle)
    R_L = 2 * np.sqrt(h**2 + (0.25 * d**2))
    k = omega / c
    dR = d - R_L
    p_r = p * (1 + vertical * np.exp(1j * k * dR)) / (1 + horizontal * np.exp(1j * k * dR))

    return p_r

def polarization_power(p_t, p_r):
    p_t_conj = np.conjugate(p_t)
    p_r_conj = np.conjugate(p_r)
    return ((1+p_t*p_r)*(1+p_t_conj*p_r_conj))/((1+p_t*p_t_conj)*(1+p_r*p_r_conj))

# Update plot based on free space
def update_plot_freespace(event):
    phi_deg = phi_slider.get()
    tau_deg = tau_slider.get()
    phi_deg = round(phi_deg, 2)
    tau_deg = round(tau_deg, 2)
    phi = np.deg2rad(phi_deg)
    tau = np.deg2rad(tau_deg)
    
    theta = np.linspace(0, 2*np.pi, 360)
    power = []
    for theta_i in theta:
            p_t = polarization_phasor(phi, tau)
            p_r = polarization_phasor(theta_i, tau)
            power.append(polarization_power(p_t, p_r))
    # Update the plot
    ax.clear()
    ax.plot(theta, power, label="Polarization state")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(0)
    ax.set_title(f'Polarization Ellipse\nRotation angle phi = {phi_deg}°\nEllipticity angle tau = {tau_deg}°')
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    fig.tight_layout()
    canvas.draw()

# Update plot for metallic reflecting plate
def update_plot_metallic(event):
    phi_deg = phi_slider.get()
    tau_deg = tau_slider.get()
    phi_deg = round(phi_deg, 2)
    tau_deg = round(tau_deg, 2)
    phi = np.deg2rad(phi_deg)
    tau = np.deg2rad(tau_deg)
    height = h[0]
    
    theta = np.linspace(0, 2*np.pi, 360)
    power = []
    for theta_i in theta:
            p_t = polarization_phasor(phi, tau)
            p_r = polarization_ratio_antenna_metallic(f, e_d, d, height, theta_i, tau)
            power.append(polarization_power(p_t, p_r))
    # Update the plot
    ax.clear()
    ax.plot(theta, power, label="Polarization state")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(0)
    ax.set_title(f'Polarization Ellipse with height = {h[0]} \nRotation angle phi = {phi_deg}°\nEllipticity angle tau = {tau_deg}°')
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    fig.tight_layout()
    canvas.draw()

# Update plot for dielectric reflecting plate
def update_plot_dielectric(event):
    phi_deg = phi_slider.get()
    tau_deg = tau_slider.get()
    phi_deg = round(phi_deg, 2)
    tau_deg = round(tau_deg, 2)
    phi = np.deg2rad(phi_deg)
    tau = np.deg2rad(tau_deg)
    height = h[0]
    

    theta = np.linspace(0, 2*np.pi, 360)
    power = []
    for theta_i in theta:
            p_t = polarization_phasor(phi, tau)
            p_r = polarization_ratio_antenna_dielectric(f, e_d, d, height, theta_i, tau)
            power.append(polarization_power(p_t, p_r))
    # Update the plot
    ax.clear()
    ax.plot(theta, power, label="Polarization state")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(0)
    ax.set_title(f'Polarization Ellipse with height = {h[0]} \nRotation angle phi = {phi_deg}°\nEllipticity angle tau = {tau_deg}°')
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    fig.tight_layout()
    canvas.draw()

# update plot based on brewster angle and dielectric reflecting plate
def update_plot_brewster(event):
    phi_deg = phi_slider.get()
    tau_deg = tau_slider.get()
    phi_deg = round(phi_deg, 2)
    tau_deg = round(tau_deg, 2)
    phi = np.deg2rad(phi_deg)
    tau = np.deg2rad(tau_deg)
    height = h[0]
    
    theta = np.linspace(0, 2*np.pi, 360)
    power = []
    for theta_i in theta:
            p_t = polarization_phasor(phi, tau)
            p_r = polarization_ratio_antenna_brewster(f, e_d, d, height, theta_i, tau)
            power.append(polarization_power(p_t, p_r))
    # Update the plot
    ax.clear()
    ax.plot(theta, power, label="Polarization state")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(0)
    ax.set_title(f'Polarization Ellipse\nRotation angle phi = {phi_deg}°\nEllipticity angle tau = {tau_deg}°')
    ax.set_aspect('equal', 'box')
    ax.grid(True)
    fig.tight_layout()
    canvas.draw()

# initially set to one of the update functions to have a default plot.
update_plot = update_plot_freespace  # Default update function


# Create the main Tkinter window
root = tk.Tk()
root.title("Polarization State Simulator")

# Configure the grid to make the window scalable
root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)

# Create a frame for the plot to make it fixed size
plot_frame = ttk.Frame(root, width=800, height=800)
plot_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
plot_frame.grid_propagate(False)

# Create the matplotlib figure with a fixed size
fig, ax = plt.subplots(figsize=(4, 4), subplot_kw={'projection': 'polar'})

canvas = FigureCanvasTkAgg(fig, master=plot_frame)
widget = canvas.get_tk_widget()
widget.pack(fill=tk.BOTH, expand=True)

# Initialize a global variable for the current update function
current_update_function = None

# This function will be called by the sliders and will delegate to the currently selected plotting function
def update_plot_delegate(event):
    if current_update_function is not None:
        current_update_function(event)

# This function changes the current plotting function and updates the plot
def change_update_function(new_function):
    global current_update_function
    current_update_function = new_function
    # Trigger the update to reflect the change immediately
    update_plot_delegate(None)


# Sliders
slider_frame = ttk.Frame(root)
slider_frame.grid(row=1, column=0, sticky="ew")
root.grid_rowconfigure(1, weight=0)

# Connect the sliders to the delegate function instead of update_plot directly
phi_slider = ttk.Scale(slider_frame, from_=0, to_=180, orient='horizontal', command=update_plot_delegate)
phi_slider.grid(row=0, column=1, sticky="ew", padx=5)
ttk.Label(slider_frame, text="Phi (°)").grid(row=0, column=0, padx=5)


tau_slider = ttk.Scale(slider_frame, from_=0, to_=180, orient='horizontal', command=update_plot_delegate)
tau_slider.grid(row=1, column=1, sticky="ew", padx=5)
ttk.Label(slider_frame, text="Tau (°)").grid(row=1, column=0, padx=5)

slider_frame.grid_columnconfigure(1, weight=1)

# Add buttons for different polarization states
button_frame = ttk.Frame(root)
button_frame.grid(row=2, column=0, sticky="ew")
root.grid_rowconfigure(2, weight=0)

# Creating and placing buttons
freespace_button = ttk.Button(button_frame, text="Free Space", command=lambda: change_update_function(update_plot_freespace))
freespace_button.grid(row=0, column=0, padx=5, pady=5)

metallic_button = ttk.Button(button_frame, text="Metallic", command=lambda: change_update_function(update_plot_metallic))
metallic_button.grid(row=0, column=1, padx=5, pady=5)

dielectric_button = ttk.Button(button_frame, text="Dielectric", command=lambda: change_update_function(update_plot_dielectric))
dielectric_button.grid(row=0, column=2, padx=5, pady=5)

brewster_button = ttk.Button(button_frame, text="Brewster Angle", command=lambda: change_update_function(update_plot_brewster))
brewster_button.grid(row=0, column=3, padx=5, pady=5)

root.mainloop()
