import matplotlib.pyplot as plt
import numpy as np
import h5py

file_name = "\\\\filer\\userfiles\\htb27\\windows_home\\ChipWhispherer\\first_capture_8bit.hdf5"
dataset_name = "0"

power_wave = None
trigger_wave = None

with h5py.File(file_name, "r") as file:
    power_wave = file["power_" + dataset_name]
    #power_wave_2 = file["power_2"]
    #trigger_wave = file["trigger_" + dataset_name]

    fig, axes = plt.subplots(figsize=(15,10), sharex=True)

#    axes[0].plot(trigger_wave[...])
 #   axes[0].set_title("Trigger")

    axes.plot(power_wave[...])
    axes.set_title("Power trace")
    #axes[1].plot(power_wave_2[...])
    
    print(f"Power gain: {power_wave.attrs['gain']}")
    print(f"Power offset: {power_wave.attrs['offset']}")
    print(f"Power key: {power_wave.attrs['key']}")
    print(f"Power nonce: {power_wave.attrs['nonce']}")
    print(f"Power counter: {power_wave.attrs['counter']}")
    print(f"Power plaintext: {power_wave.attrs['plaintext']}")

plt.show()
