import CWController
import NIController
import niscope, nifgen
import numpy as np
import matplotlib.pyplot as plt
import time
from numpy.lib.format import open_memmap
import h5py

def convert_byte_array_to_integer_list(byte_array):
    return [int.from_bytes(byte_array[i:i + 4], 'little') for i in range(0, len(byte_array), 4)]


# ChipWhisperer Settings
firmware_path = None
F_CPU = int(5e6)  # Compile time option
target_frequency = int(5e6)

# NI Scope settings
sample_rate = int(2500e6)  # change sample_rate here
points_per_clock = sample_rate // target_frequency
number_of_cycles = 1500  # Scope length by cycle count
number_of_points = int(points_per_clock * number_of_cycles)
trace_channel = "0"
trigger_channel = "1"
channel_0_range = 0.03  # peak-to-peak voltage (by using high-pass filter) was 0.04 but looking at the data 0.003 looks like it should be able to work without clipping


def init(firmware_path=None, target_frequency=target_frequency):
    chip_whisperer_controller = setup_chip_whisperer(
        firmware_path
    )  # target communication, program firmware
    ni_controller = setup_NI_Fgen(chip_whisperer_controller, target_frequency)  # clock
    ni_controller = setup_NI_scope(ni_controller)  # oscilloscope
    return chip_whisperer_controller, ni_controller


def setup_chip_whisperer(firmware_path=None):
    print("Connect to ChipWhisperer")
    chip_whisperer_controller = CWController.CWController(f_cpu=F_CPU)
    time.sleep(0.5)

    if firmware_path is not None:
        print("Programming target")
        chip_whisperer_controller.program_target(firmware_path)

    chip_whisperer_controller.set_scope_clk(target_frequency)
    time.sleep(0.5)

    print("Reseting the target")
    chip_whisperer_controller.reset_target()

    print("Checking communication")
    com_check = chip_whisperer_controller.test_run(True)
    assert com_check
    return chip_whisperer_controller


def setup_NI_Fgen(chip_whisperer_controller, frequency=target_frequency):
    print("Switching to NI clock")
    chip_whisperer_controller.clk_off()
    #input("Turned off the scope clock and now connect up the waveform generator")
    ni_controller = NIController.NIController(fgen_freq=frequency)
    ni_controller.fgen.initiate()
    print("Check communication")
    chip_whisperer_controller.reset_target()
    com_check = chip_whisperer_controller.test_run()
    assert com_check, "CW target failed to recalibrate to new clock!"
    return ni_controller


def setup_NI_scope(ni_controller):
    print("Setting up NI scope")

    ni_controller.configure_trigger(trigger_channel, 1.6) # Looking at the data of the test captures I did this seems to offer the least amount of variability 

    # Configure NI Scope power measurement
    # Make sure the input is high pass filtered

    print("CH0 Impedance = 50 ohm, Make sure input is DC Blocked/High pass filtered!!")
    time.sleep(1.5)
    ni_controller.configure_ch0(channel_0_range, niscope.VerticalCoupling.DC, 50.0)

    ni_controller.configure_horizontal(
        min_sample_rate=sample_rate,
        min_num_pts=number_of_points,
        ref_position=1.0,
        num_records=1,
        enforce_realtime=True,
    )

    print("Finished setting up NI scope")
    time.sleep(0.5)
    return ni_controller


def capture_raw_trace(
    chip_whisperer_controller: CWController.CWController,
    ni_controller: NIController.NIController,
    plaintext,
    key,
    nonce,
    counter,
    group,
    dataset_name,
    record_trigger=False,
):
    # Open a memmap file, write directly to disk
    power_wave = np.ndarray(number_of_points, dtype=np.int16)
    if record_trigger:
        trigger_wave = np.ndarray(number_of_points, dtype=np.int16)
    chip_whisperer_controller.reset_cipher()
    chip_whisperer_controller.set_key(key)
    chip_whisperer_controller.set_nonce(nonce)
    chip_whisperer_controller.set_counter(counter)
    chip_whisperer_controller.set_plaintext(plaintext)

    ni_controller.arm()
    chip_whisperer_controller.perform_encryption()
    power_waveforms = ni_controller.scope.channels[0].fetch_into(power_wave, num_records = 1, timeout = 20.0)
    if record_trigger:
        trigger_waveforms = ni_controller.scope.channels[1].fetch_into(trigger_wave, num_records = 1, timeout = 20.0)
    
    power_trace = group.create_dataset("power_" + dataset_name, (number_of_points, ), dtype='i2', compression="gzip", compression_opts=9)
    power_trace[...] = power_waveforms[0].samples
    power_trace.attrs['gain'] = power_waveforms[0].gain
    power_trace.attrs['offset'] = power_waveforms[0].offset
    power_trace.attrs['key'] = key
    power_trace.attrs['nonce'] = nonce
    power_trace.attrs['counter'] = counter
    power_trace.attrs['plaintext'] = plaintext

    if record_trigger:
        trigger_trace = group.create_dataset("trigger_" + dataset_name, (number_of_points, ), dtype='i2', compression="gzip", compression_opts=9)
        trigger_trace[...] = trigger_waveforms[0].samples
        trigger_trace.attrs['gain'] = trigger_waveforms[0].gain
        trigger_trace.attrs['offset'] = trigger_waveforms[0].offset
        trigger_trace.attrs['key'] = key
        trigger_trace.attrs['nonce'] = nonce
        trigger_trace.attrs['counter'] = counter


def main(file_recording_number):
    chip_whisperer_controller, ni_controller = init()

    print("Capture several traces into the file")

    with h5py.File(f"\\\\filer\\userfiles\\htb27\\unix_home\\ChaChaRecordings\\recording_attack_counter_from_one_{file_recording_number}.hdf5", "a") as file:
        for i in range(100):
            print(f"Key number: {i}")
            key = convert_byte_array_to_integer_list(bytearray(np.random.bytes(32)))
            nonce = convert_byte_array_to_integer_list(bytearray(np.random.bytes(12)))

            for j in range(10):
                capture_raw_trace(
                    chip_whisperer_controller=chip_whisperer_controller,
                    ni_controller=ni_controller,
                    plaintext=convert_byte_array_to_integer_list(bytearray(np.random.bytes(64))),
                    key=key,
                    nonce=nonce,
                    counter=[j + 1],
                    group=file,
                    dataset_name=str(i) + "_" + str(j),
                    record_trigger=False
                )
    print("Captured raw trace")


if __name__ == "__main__":
    main(1)
