import math
import random
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# ==========================================
# 1. Configuration Dataclasses
# ==========================================

# Original waveform struct[cite: 1]
@dataclass
class WaveformData:
    tot_n_samples: int = 0      #[cite: 1]
    init_ped_samples: int = 0   #[cite: 1]
    pedestal: int = 0           #[cite: 1]
    ped_range: int = 0          #[cite: 1]
    n_pulses: int = 0           #[cite: 1]
    pulse_amp: int = 0          #[cite: 1]
    pulse_width: int = 0        #[cite: 1]
    pulse_separation: int = 0   #[cite: 1]

# New Detector & Digitizer Structs (Constructor Injection)
@dataclass
class DetectorParams:
    tau_r: float = 2.0
    tau_f: float = 10.0
    tts_sigma: float = 1.5
    mean_gain: float = 5.0
    gain_sigma: float = 1.0
    ap_prob: float = 0.03
    ap_delay_mean: float = 80.0
    ap_delay_sigma: float = 5.0

@dataclass
class DigitizerParams:
    sample_rate: float = 4.0
    num_samples: int = 256
    baseline: float = 100.0
    noise_sigma: float = 2.5
    v_lsb: float = 1.0
    max_adc: int = 4095

# ==========================================
# 2. Waveform Generator Class
# ==========================================
class WaveformGenerator:
    def __init__(self, det_params: DetectorParams = None, daq_params: DigitizerParams = None):
        random.seed() #[cite: 1]
        
        # Initialize original data structure[cite: 1]
        self.waveform_data = WaveformData() 
        
        # Initialize physics configurations
        self.det_params = det_params if det_params else DetectorParams()
        self.daq_params = daq_params if daq_params else DigitizerParams()

    # Setter for the original waveform logic[cite: 1]
    def set_waveform_data(self, waveform_data: WaveformData):
        self.waveform_data = waveform_data #[cite: 1]

    # ---------------------------------------------------------
    # Algorithm 1: Original Multiple Pulse Waveform
    # ---------------------------------------------------------
    def multiple_pulse_waveform(self) -> list[int]:
        min_ped = self.waveform_data.pedestal - self.waveform_data.ped_range // 2 #[cite: 1]
        max_ped = self.waveform_data.pedestal + self.waveform_data.ped_range // 2 + 1 #[cite: 1]
        
        num_pulses = random.randint(2, self.waveform_data.n_pulses) #[cite: 1]
        pulse_info = [] #[cite: 1]
        
        # Calculate maximum spacing and starting points[cite: 1]
        max_total_space = (num_pulses * self.waveform_data.pulse_width) + \
                          ((num_pulses - 1) * self.waveform_data.pulse_separation) #[cite: 1]
                            
        max_first_start = self.waveform_data.tot_n_samples - max_total_space #[cite: 1]
        
        if max_first_start < self.waveform_data.init_ped_samples: #[cite: 1]
            max_first_start = self.waveform_data.init_ped_samples #[cite: 1]
            
        current_start = random.randint(self.waveform_data.init_ped_samples, max_first_start) #[cite: 1]
        
        # Generate random amplitudes and separations for each pulse[cite: 1]
        for _ in range(num_pulses):
            min_amp = self.waveform_data.pulse_amp // 2 #[cite: 1]
            current_amp = random.randint(min_amp, self.waveform_data.pulse_amp) #[cite: 1]
            pulse_info.append((current_start, current_amp)) #[cite: 1]
            
            # Random separation[cite: 1]
            current_separation = random.randint(self.waveform_data.pulse_separation//2, self.waveform_data.pulse_separation) #[cite: 1]
            current_start += self.waveform_data.pulse_width + current_separation #[cite: 1]
            
        wave = [] #[cite: 1]
        for i in range(self.waveform_data.tot_n_samples): #[cite: 1]
            in_pulse = False #[cite: 1]
            pulse_relative_idx = 0 #[cite: 1]
            active_amp = 0 #[cite: 1]
            
            for start_idx, amp in pulse_info: #[cite: 1]
                if start_idx <= i < start_idx + self.waveform_data.pulse_width: #[cite: 1]
                    in_pulse = True #[cite: 1]
                    pulse_relative_idx = i - start_idx #[cite: 1]
                    active_amp = amp #[cite: 1]
                    break #[cite: 1]
                    
            if in_pulse: #[cite: 1]
                wave.append(min_ped + active_amp // (2 ** pulse_relative_idx)) #[cite: 1]
            else: #[cite: 1]
                wave.append(random.randrange(min_ped, max_ped)) #[cite: 1]
                
        return wave #[cite: 1]

    # ---------------------------------------------------------
    # Algorithm 2: Physics-Based PE vs. Time Waveform
    # ---------------------------------------------------------
    def _spe_pulse(self, t: float, t0: float, amplitude: float) -> float:
        """Helper to calculate the voltage of a single photoelectron."""
        if t < t0:
            return 0.0
        return amplitude * (math.exp(-(t - t0) / self.det_params.tau_f) - 
                            math.exp(-(t - t0) / self.det_params.tau_r))

    def _generate_pe_times(self) -> list[float]:
        """Helper to simulate physics events: multi-hits and afterpulsing."""
        num_pulses = random.randint(1, 4)
        window_max = self.daq_params.num_samples * self.daq_params.sample_rate - 50.0
        
        all_pe_times = []
        for _ in range(num_pulses):
            hit_time = random.uniform(10.0, window_max)
            # Using numpy for Poisson distribution of PEs per event
            num_pes = np.random.poisson(15) 
            
            for _ in range(num_pes):
                all_pe_times.append(hit_time)
                
                # Afterpulse roll
                if random.random() < self.det_params.ap_prob:
                    delayed_time = hit_time + random.gauss(self.det_params.ap_delay_mean, self.det_params.ap_delay_sigma)
                    all_pe_times.append(delayed_time)
                    
        all_pe_times.sort()
        return all_pe_times

    def generate_pe_waveform(self) -> list[int]:
        """Main public method to generate the digitized FADC waveform."""
        pe_times = self._generate_pe_times()
        analog_waveform = [0.0] * self.daq_params.num_samples
        
        # Convolution & Stochastic Smearing
        for true_time in pe_times:
            smeared_time = true_time + random.gauss(0.0, self.det_params.tts_sigma)
            smeared_gain = random.gauss(self.det_params.mean_gain, self.det_params.gain_sigma)
            if smeared_gain < 0.0:
                smeared_gain = 0.0
                
            for i in range(self.daq_params.num_samples):
                current_time = i * self.daq_params.sample_rate
                # Optimization bounds
                if (smeared_time - 5.0) < current_time < (smeared_time + 50.0):
                    analog_waveform[i] += self._spe_pulse(current_time, smeared_time, smeared_gain)
                    
        # Noise Injection & Digitization
        adc_waveform = [0] * self.daq_params.num_samples
        for i in range(self.daq_params.num_samples):
            noisy_voltage = analog_waveform[i] + random.gauss(0.0, self.daq_params.noise_sigma)
            adc_count = int(math.floor(noisy_voltage / self.daq_params.v_lsb) + self.daq_params.baseline)
            # ADC Saturation Check
            adc_count = max(0, min(adc_count, self.daq_params.max_adc))
            adc_waveform[i] = adc_count
            
        return adc_waveform


# ==========================================
# 3. Main Execution and Plotting
# ==========================================
def main():
    # Setup for original logic[cite: 1]
    wave_data = WaveformData(
        tot_n_samples=256,
        init_ped_samples=4, #[cite: 1]
        pedestal=100,       #[cite: 1]
        ped_range=4,        #[cite: 1]
        n_pulses=4,         #[cite: 1]
        pulse_amp=300,      #[cite: 1]
        pulse_width=4,      #[cite: 1]
        pulse_separation=20 #[cite: 1]
    )

    # Setup for new PE logic
    det_params = DetectorParams(ap_prob=0.05, mean_gain=6.0)
    daq_params = DigitizerParams(num_samples=256)

    # Initialize generator and inject structs
    generator = WaveformGenerator(det_params=det_params, daq_params=daq_params)
    generator.set_waveform_data(wave_data) #[cite: 1]

    # Generate both types of waveforms
    orig_multiple_pulse = generator.multiple_pulse_waveform()
    physics_pe_waveform = generator.generate_pe_waveform()

    # Plotting Output
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Original Algorithm Plot
    ax1.plot(orig_multiple_pulse, marker='o', markersize=3, linestyle='-', color='g', label='Original Algorithm')
    ax1.set_title("Uploaded Algorithm: Simple Multiple Pulse Waveform")
    ax1.set_xlabel("Sample Number")
    ax1.set_ylabel("ADC Counts")
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()
    
    # Physics PE Algorithm Plot
    # Scaling X axis to time for the physics plot based on sample rate
    time_axis = [i * daq_params.sample_rate for i in range(daq_params.num_samples)]
    ax2.plot(time_axis, physics_pe_waveform, marker='o', markersize=3, linestyle='-', color='b', label='Physics PE Convolution')
    ax2.set_title("New Algorithm: Physics-Based PE vs. Time with Afterpulsing")
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("ADC Counts")
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()