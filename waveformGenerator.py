import random
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class WaveformData:
    tot_n_samples: int = 0
    init_ped_samples: int = 0
    pedestal: int = 0
    ped_range: int = 0
    n_pulses: int = 0
    pulse_amp: int = 0
    pulse_width: int = 0
    pulse_separation: int = 0

class WaveformGenerator:
    def __init__(self):
        random.seed()
        self.waveform_data = WaveformData()

    def set_waveform_data(self, waveform_data: WaveformData):
        self.waveform_data = waveform_data

    # Noise waveform
    def noise_waveform(self) -> list[int]:
        min_val = self.waveform_data.pedestal - self.waveform_data.ped_range // 2
        max_val = self.waveform_data.pedestal + self.waveform_data.ped_range // 2 + 1
        
        wave = []
        for _ in range(self.waveform_data.tot_n_samples):
            wave.append(random.randrange(min_val, max_val))
            
        return wave

    # Single pulse waveform
    def single_pulse_waveform(self) -> list[int]:
        min_ped = self.waveform_data.pedestal - self.waveform_data.ped_range // 2
        max_ped = self.waveform_data.pedestal + self.waveform_data.ped_range // 2 + 1
        
        # Generate a random start index for the pulse anywhere AFTER the initial pedestal buffer
        start_min = self.waveform_data.init_ped_samples
        start_max = self.waveform_data.tot_n_samples - 1
        
        # Ensure we don't error out if properties are configured weirdly
        if start_max < start_min:
            start_max = start_min
            
        random_start_idx = random.randint(start_min, start_max)
        
        wave = []
        for i in range(self.waveform_data.tot_n_samples):
            # Check if current index 'i' falls within the active pulse window
            if random_start_idx <= i < random_start_idx + self.waveform_data.pulse_width:
                j = i - random_start_idx
                
                # This is a deterministic curve without noise added on top. 
                wave.append(min_ped + self.waveform_data.pulse_amp // 2 ** (j))
            else:
                # Outside the pulse window, just generate background noise
                wave.append(random.randrange(min_ped, max_ped))
                
        return wave

    # Multiple pulse waveform
    def multiple_pulse_waveform(self) -> list[int]:
        min_ped = self.waveform_data.pedestal - self.waveform_data.ped_range // 2
        max_ped = self.waveform_data.pedestal + self.waveform_data.ped_range // 2 + 1
        
        num_pulses = random.randint(2, self.waveform_data.n_pulses)

        pulse_info = []
        
        max_total_space = (num_pulses * self.waveform_data.pulse_width) + \
                          ((num_pulses - 1) * self.waveform_data.pulse_separation)
                            
        max_first_start = self.waveform_data.tot_n_samples - max_total_space
        
        if max_first_start < self.waveform_data.init_ped_samples:
            max_first_start = self.waveform_data.init_ped_samples 
            
        current_start = random.randint(self.waveform_data.init_ped_samples, max_first_start)
        
        # Generate random amplitudes and separations for each pulse
        for _ in range(num_pulses):
            min_amp = self.waveform_data.pulse_amp // 2
            current_amp = random.randint(min_amp, self.waveform_data.pulse_amp)

            pulse_info.append((current_start, current_amp))
            
            # Random seperation
            current_separation = random.randint(self.waveform_data.pulse_separation//2, self.waveform_data.pulse_separation)
            current_start += self.waveform_data.pulse_width + current_separation
            
        wave = []
        for i in range(self.waveform_data.tot_n_samples):
            in_pulse = False
            pulse_relative_idx = 0
            active_amp = 0
            
            for start_idx, amp in pulse_info:
                if start_idx <= i < start_idx + self.waveform_data.pulse_width:
                    in_pulse = True
                    pulse_relative_idx = i - start_idx
                    active_amp = amp 
                    break
                    
            if in_pulse:
                wave.append(min_ped + active_amp // (2 ** pulse_relative_idx))
            else:
                wave.append(random.randrange(min_ped, max_ped))
                
        return wave

# Main function
def main():
    # Define WaveformData
    wave_data = WaveformData(
        tot_n_samples=120,
        init_ped_samples=4, # Minimum number of samples for initial pedestal
        pedestal=100,
        ped_range=4, # This is the sigma range for pedestal
        n_pulses=4,
        pulse_amp=300,
        pulse_width=4,
        pulse_separation=20
    )

    # Initialize generator and set data
    generator = WaveformGenerator()
    generator.set_waveform_data(wave_data)

    noise = generator.noise_waveform()
    single_pulse = generator.single_pulse_waveform()
    multiple_pulse = generator.multiple_pulse_waveform()

    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharey=True)
    
    # Pure Noise Plot
    ax1.plot(noise, marker='o', linestyle='-', color='b', label='Noise Only')
    ax1.set_title("Generated Noise Waveform")
    ax1.set_ylabel("Amplitude")
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()
    
    # Single Pulse Plot
    ax2.plot(single_pulse, marker='o', linestyle='-', color='r', label='Single Pulse + Noise')
    ax2.set_title("Generated Single Pulse Waveform")
    ax2.set_xlabel("Sample Index")
    ax2.set_ylabel("Amplitude")
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()
    
    # Multiple Pulse Plot
    ax3.plot(multiple_pulse, marker='o', linestyle='-', color='g', label='Multiple Pulse + Noise')
    ax3.set_title("Generated Multiple Pulse Waveform")
    ax3.set_xlabel("Sample Index")
    ax3.set_ylabel("Amplitude")
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.legend()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()