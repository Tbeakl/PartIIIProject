# Taken from Eric's code
import niscope, nifgen

class NIController:
    
    def __init__(self, fgen_freq=5000000):
        self.scope = None
        self.fgen  = None
        self.scope_default_setup()
        self.fgen_default_setup(fgen_freq)

    def __del__(self):
        if self.fgen is not None:
            self.fgen.close()
        if self.scope is not None:
            self.scope.close()
    

    ### Fgen control ###
    def fgen_default_setup(self, freq=5000000):
        self.fgen = nifgen.Session('awg')
        self.fgen.reference_clock_source = nifgen.ReferenceClockSource.PXI_CLOCK
        self.fgen.output_mode = nifgen.OutputMode.FUNC
        self.fgen.configure_standard_waveform(waveform=nifgen.Waveform.SQUARE, 
                amplitude=1.6, frequency=freq, dc_offset=0.8) # Think that Shih-Chun has it with amplitude of 3 and offset of 1.5

    def fgen_init(self):
        self.fgen.initiate()

    def fgen_abort(self):
        self.fgen.abort()
    

    ### Scope control ###
    def scope_default_setup(self):
        self.scope = niscope.Session(resource_name='scope')
        self.scope.input_clock_source = 'VAL_PXI_CLOCK'
        self.configure_horizontal = self.scope.configure_horizontal_timing
        self.fetch_ch0 = self.scope.channels[0].fetch

        # configure trigger
        self.scope.channels[1].input_impedance = 1000000.0
        self.scope.channels[1].configure_vertical(range=5, 
                               coupling=niscope.VerticalCoupling.DC, 
                               offset=1.5, probe_attenuation=1.0) # Shih-Chun uses range=10 and no offset
        self.configure_trigger('1', 1.4)

        # configure power measurement
        self.scope.channels[0].input_impedance = 1000000.0
        self.scope.channels[0].configure_vertical(range=0.1,
                coupling=niscope.VerticalCoupling.DC)

        # configure sampling rate
        self.scope.configure_horizontal_timing(min_sample_rate=2.5E9, 
                   min_num_pts=100000000, ref_position=5.0, num_records=1, 
                   enforce_realtime=True) #May want to change the ref_position


    def configure_trigger(self, trigger_source, trigger_level, 
                          trigger_coupling=niscope.TriggerCoupling.DC):
        self.scope.configure_trigger_edge(trigger_source, level=trigger_level,
                   trigger_coupling=niscope.TriggerCoupling.DC)

    def configure_ch0(self, range, coupling, input_impedance=1E6):
        self.scope.channels[0].input_impedance = input_impedance # float
        self.scope.channels[0].configure_vertical(range=range,coupling=coupling)

    def configure_ch1(self, range, coupling, input_impedance=1E6):
        self.scope.channels[1].input_impedance = input_impedance # float
        self.scope.channels[1].configure_vertical(range=range,coupling=coupling)

    def arm(self):
        if self.scope.acquisition_status() == niscope.AcquisitionStatus.IN_PROGRESS:
            self.scope.abort()
        self.scope.initiate()
