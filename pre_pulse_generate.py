import numpy as np
import matplotlib.pyplot as plt

# Define the full time array: 2 ns total duration with 1 fs resolution
full_time = np.arange(0, 2e-9, 2e-15)  # 0 to 2 ns with 1 fs steps

# Convert FWHM to standard deviation
fwhm = 25e-15  # 25 fs
std_dev = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation

# Gaussian parameters for the profiles
amplitude1 = 1e16  # Amplitude in W/cm^2
mean1 = 0.2e-9     # Center at 1.6 ns
amplitude2 = 1e16  # Amplitude in W/cm^2
mean2 = 0.4e-9     # Center at 1.8 ns

# Define the Gaussian profiles with the calculated standard deviation
gaussian1 = amplitude1 * np.exp(-((full_time - mean1) ** 2) / (2 * std_dev ** 2))
gaussian2 = amplitude2 * np.exp(-((full_time - mean2) ** 2) / (2 * std_dev ** 2))

# Merge the two profiles into one array
merged_profile = gaussian1 + gaussian2

# Convert the profile from W/cm^2 to erg/cm^2
merged_profile_erg = merged_profile * 1e7  # Conversion factor

# Remove the first 1.5 ns and shift the time axis
relevant_time = full_time[full_time >= 0e-9]
relevant_profile = merged_profile_erg[full_time >= 0e-9]

print(relevant_profile.size)

# # Save the shifted time axis and the merged profile as txt files
# np.savetxt('shifted_time_axis.txt', relevant_time, header='Time (s)', fmt='%.5e')
# np.savetxt('shifted_merged_profile.txt', relevant_profile, header='Intensity (erg/cm^2)', fmt='%.5e')

# Save both time and profile into one file
output_data = np.column_stack((relevant_time, relevant_profile))
wkdir = '/Users/yao/Nextcloud/PROJECTS/Apollon/Shock_Oct2024/Apollon_shock2024_prepulse_1e16_25fs_CH/'
np.savetxt(wkdir+'time_and_profile.txt', output_data, header='Time (s)\tIntensity (erg/cm^2)', fmt='%.5e')


# Plotting the relevant portion of the merged Gaussian profile
plt.figure(figsize=(10, 6))
plt.plot(relevant_time * 1e12, relevant_profile, label='Merged Gaussian Profile', color='purple')  # Time in ps for plotting
plt.xlabel('Time (ps)')
plt.ylabel('Intensity (erg/cm²)')
plt.title('Merged Gaussian Profile as a Function of Shifted Time')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()