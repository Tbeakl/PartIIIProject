import time
import numpy as np
import chipwhisperer as cw

fw_path = "./firmware/simpleserial-base-CWLITEARM.hex"

def turn_integer_list_into_byte_array(integers):
    return bytearray(b''.join([(i).to_bytes(4, 'little') for i in integers]))

def convert_byte_array_to_integer_list(byte_array):
    return [int.from_bytes(byte_array[i:i + 4], 'little') for i in range(0, len(byte_array), 4)]


class CWController:
    def __init__(self, f_cpu=5E6):
        self.scope  = None
        self.target = None
        self.f_cpu = f_cpu  # the F_CPU in CW Makefile.platform, may influence target.baud
        self.connect()
        return

    def __del__(self):
        if self.target is not None:
            self.target.dis()
        if self.scope is not None:
            self.scope.dis()
        return

    def connect(self):
        self.scope = cw.scope()
        target_type = cw.targets.SimpleSerial

        try:
            self.target = cw.target(self.scope, target_type)
        except IOError:
            self.scope = cw.scope()
            self.target = cw.target(self.scope, target_type)
        time.sleep(0.05)
        self.scope.default_setup()
        self.reset_scope_clk(self.f_cpu)
        time.sleep(0.5)
        return

    def program_target(self, fw_path=fw_path):
        prog = cw.programmers.STM32FProgrammer
        cw.program_target(self.scope, prog, fw_path)
        time.sleep(0.5)
        self.reset_target()
        print('flush:', self.target.read())

    def reset_target(self):
        self.scope.io.nrst = 'low'
        time.sleep(0.25)
        self.scope.io.nrst = 'high_z'
        time.sleep(0.25)
        self.target.read()

    def reset_scope_clk(self, freq=5E6): # set freq to predesigned F_CPU
        self.scope.clock.adc_src = 'clkgen_x1'
        self.scope.clock.clkgen_freq = freq
        print(self.scope.clock.clkgen_mul, self.scope.clock.clkgen_div)
        self.scope.clock.freq_ctr_src = 'clkgen'
        time.sleep(0.5)
        print(self.scope.clock.freq_ctr)
        self.scope.clock.reset_dcms()
        self.scope.clock.reset_clkgen()
        self.scope.clock.reset_adc()
        time.sleep(0.01)
        self.reset_target()

    def set_scope_clk(self, freq, show_clk_info=False):
        old_freq = self.scope.clock.clkgen_freq
        self.scope.clock.clkgen_freq = int(freq)
        time.sleep(0.5)
        self.scope.clock.reset_dcms()
        self.scope.clock.reset_clkgen()
        self.scope.clock.reset_adc()
        time.sleep(0.5)
        self.target.baud = int(self.target.baud*freq/old_freq)
        #self.target.baud *= self.scope.clock.clkgen_freq/old_freq
        self.reset_target()
        print('flush:',self.target.read())
        if show_clk_info:
            print(self.scope.clock)

    def clk_off(self):
        self.scope.io.hs2 = "disabled"
    
    def close(self):
        self.target.dis()
        self.scope.dis()

    ### ChaCha ####################################

    def test_run(self, show_output=True):
        show_output and print("Running test script")
        flush = self.target.read()
        show_output and print(f"Flush: {flush}")
        show_output and print(f"Setting keys, nonce and counter")
        self.set_key([50462976, 117835012, 185207048, 252579084, 319951120, 387323156, 454695192, 522067228])
        self.set_nonce([150994944, 1241513984, 0])
        self.set_counter([1])
        self.perform_encryption()
        actual_output = self.read_final_state()
        expected_output = bytearray([0x10,0xf1,0xe7,0xe4,0xd1,0x3b,0x59,0x15,0x50,0x0f,0xdd,0x1f,0xa3,0x20,0x71,0xc4,0xc7,0xd1,0xf4,0xc7,0x33,0xc0,0x68,0x03,0x04,0x22,0xaa,0x9a,0xc3,0xd4,0x6c,0x4e,0xd2,0x82,0x64,0x46,0x07,0x9f,0xaa,0x09,0x14,0xc2,0xd7,0x05,0xd9,0x8b,0x02,0xa2,0xb5,0x12,0x9c,0xd1,0xde,0x16,0x4e,0xb9,0xcb,0xd0,0x83,0xe8,0xa2,0x50,0x3c,0x4e])
        show_output and print(f"Check final state == expected: {actual_output == expected_output}")
        self.reset_cipher()
        return actual_output == expected_output

    def reset_cipher(self):
        self.target.simpleserial_write('x', bytearray([]))
        self.target.simpleserial_wait_ack()
    
    def set_key(self, key):
        if type(key) is not bytearray:
            key = turn_integer_list_into_byte_array(key)
        self.target.simpleserial_write('k', key)
        self.target.simpleserial_wait_ack()

    def set_nonce(self, nonce):
        if type(nonce) is not bytearray:
            nonce = turn_integer_list_into_byte_array(nonce)
        self.target.simpleserial_write('n', nonce)
        self.target.simpleserial_wait_ack()

    def set_counter(self, counter):
        if type(counter) is not bytearray:
            counter = turn_integer_list_into_byte_array(counter)
        self.target.simpleserial_write('c', counter)
        self.target.simpleserial_wait_ack()

    def set_plaintext(self, plaintext):
        if type(plaintext) is not bytearray:
            plaintext = turn_integer_list_into_byte_array(plaintext)
        self.target.simpleserial_write('p', plaintext)
        self.target.simpleserial_wait_ack()

    def read_input_state(self):
        self.target.simpleserial_write('i', bytearray([]))
        return self.target.simpleserial_read('r', 64)
    
    def read_plaintext(self):
        self.target.simpleserial_write('m', bytearray([]))
        return self.target.simpleserial_read('r', 64)

    def perform_encryption(self):
        self.target.simpleserial_write('e', bytearray([]))
        self.target.simpleserial_wait_ack()

    def read_ciphertext(self):
        self.target.simpleserial_write('o', bytearray([]))
        return self.target.simpleserial_read('r', 64)

    def read_final_state(self):
        self.target.simpleserial_write('f', bytearray([]))
        return self.target.simpleserial_read('r', 64)
