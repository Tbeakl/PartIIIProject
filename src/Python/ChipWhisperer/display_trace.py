import matplotlib.pyplot as plt
import numpy as np

file_name = "test_captures/30_raw_capture.npy"

waves = np.load(file_name)

fig, axes = plt.subplots(figsize=(15,10), nrows=2, sharex=True)

axes[0].plot(waves[0])
axes[0].set_title("Trigger")

axes[1].plot(waves[1])
axes[1].set_title("Power trace")


plt.show()
