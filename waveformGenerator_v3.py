import math
import random
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class DetectorParams:
    tau_r: float = 2.0
    tau_f: float = 5.0
    tts_sigma: float = 0.5
    mean_gain: float = 10.0
    gain_sigma: float = 1.0
    ap_prob: float = 0.03
    ap_delay_mean: float = 80.0
    ap_delay_sigma: float = 5.0
    pulse_width: float = 5.0       # Square pulse width (samples) — from WaveformParams

@dataclass
class DigitizerParams:
    sample_rate: float = 4.0
    num_samples: int = 256
    baseline: float = 100.0
    noise_sigma: float = 2.5
    v_lsb: float = 1.0
    max_adc: int = 4095
    max_n_pulses: int = 4

class WaveformGenerator:
    def __init__(self, det_params: DetectorParams = None, daq_params: DigitizerParams = None):
        random.seed()

        self.det_params = det_params if det_params else DetectorParams()
        self.daq_params = daq_params if daq_params else DigitizerParams()

    # ---------------------------------------------------------
    # Pulse Shape Templates
    # ---------------------------------------------------------
    def _spe_pulse(self, t: float, t0: float, amplitude: float) -> float:
        """SPE pulse: exponential rise + decay (from spePulse in C++)."""
        if t < t0:
            return 0.0
        return amplitude * (math.exp(-(t - t0) / self.det_params.tau_f) -
                            math.exp(-(t - t0) / self.det_params.tau_r))

    def _square_pulse(self, t: float, t0: float, amplitude: float) -> float:
        """Square pulse: flat top between t0 and t0+width_ns (from squarePulse in C++)."""
        width_ns = self.det_params.pulse_width * (self.daq_params.sample_rate)  # samples -> ns
        if t < t0 or t > t0 + width_ns:
            return 0.0
        return amplitude

    def _triangle_pulse(self, t: float, t0: float, amplitude: float) -> float:
        """Triangle pulse: symmetric triangle centered on t0 (snapped to grid)."""
        dt = self.daq_params.sample_rate
        # Snap t0 to the nearest sample grid point so the peak is always sampled
        t0 = round(t0 / dt) * dt
        width_ns = self.det_params.pulse_width * dt  # samples -> ns
        half_w = width_ns / 2.0
        if t < t0 - half_w or t > t0 + half_w:
            return 0.0

        if t < t0:
            return amplitude * (t - (t0 - half_w)) / half_w
        else:
            return amplitude * ((t0 + half_w) - t) / half_w

    # ---------------------------------------------------------
    # PE Time Generator (shared by all waveform types)
    # ---------------------------------------------------------
    def _generate_pe_times(self) -> list[tuple[float, int]]:
        """Simulate physics events: multi-hits.

        Returns a list of (hit_time, n_pe) tuples — one entry per physics hit.
        TTS smearing is intentionally NOT applied here; it is applied once per
        hit in each waveform generator so that all PEs from the same hit share
        a single time offset, avoiding spurious secondary peaks on the rising edge.
        """
        num_pulses = random.randint(1, self.daq_params.max_n_pulses)
        window_max = self.daq_params.num_samples * self.daq_params.sample_rate - 50.0

        hits = []
        for _ in range(num_pulses):
            hit_time = random.uniform(10.0, window_max)
            num_pes = max(1, np.random.poisson(15))  # Avg 15 primary PEs per hit
            hits.append((hit_time, num_pes))

        # Sort by arrival time so the caller can process chronologically
        hits.sort(key=lambda h: h[0])
        return hits

    # --------------
    # SPE Waveform
    # --------------
    def generate_pe_waveform(self) -> list[int]:
        """Generate digitized FADC waveform using SPE pulse shape.

        TTS smearing is applied once per hit (not per individual PE) so that
        all PEs from the same physics event are time-shifted together. The
        per-PE gain fluctuations are summed into a single effective amplitude
        for the hit. This prevents spurious shoulders/secondary peaks on the
        rising edge that would otherwise arise from superimposing many
        individually time-shifted SPE pulses.
        """
        hits = self._generate_pe_times()   # list of (hit_time, n_pe)
        analog_waveform = [0.0] * self.daq_params.num_samples

        for hit_time, n_pe in hits:
            # One TTS smear for the whole hit — all PEs arrive together.
            smeared_time = hit_time + random.gauss(0.0, self.det_params.tts_sigma)

            # Sum individual PE gains into one effective amplitude.
            effective_amplitude = sum(
                random.gauss(self.det_params.mean_gain, self.det_params.gain_sigma)
                for _ in range(n_pe)
            )
            if effective_amplitude < 0.0:
                effective_amplitude = 0.0

            # Evaluate a single SPE-shaped pulse scaled by the total charge.
            for i in range(self.daq_params.num_samples):
                current_time = i * self.daq_params.sample_rate
                if (smeared_time - 5.0) < current_time < (smeared_time + 50.0):
                    analog_waveform[i] += self._spe_pulse(current_time, smeared_time, effective_amplitude)

        # Noise + Digitization
        adc_waveform = [0] * self.daq_params.num_samples
        for i in range(self.daq_params.num_samples):
            noisy_voltage = analog_waveform[i] + random.gauss(0.0, self.daq_params.noise_sigma)
            adc_count = int(math.floor(noisy_voltage / self.daq_params.v_lsb) + self.daq_params.baseline)
            adc_count = max(0, min(adc_count, self.daq_params.max_adc))
            adc_waveform[i] = adc_count

        return adc_waveform

    # -------------------
    # Square Waveform
    # -------------------
    def generate_square_waveform(self) -> list[int]:
        """Generate digitized FADC waveform using square pulse shape.
        
        This is a direct Python translation of 
        waveformGenerator_v2::squareWaveform() from the C++ code.
        """
        hits = self._generate_pe_times()   # list of (hit_time, n_pe)
        square_wf = [0.0] * self.daq_params.num_samples

        # Convolution with square pulse template
        for hit_time, n_pe in hits:
            smeared_time = hit_time          # no TTS for square (deterministic shape)
            smeared_gain = self.det_params.mean_gain * n_pe
            if smeared_gain < 0.0:
                smeared_gain = 0.0

            for i in range(self.daq_params.num_samples):
                current_time = i * self.daq_params.sample_rate
                # Optimization: only evaluate near the pulse
                if (smeared_time - 5.0) < current_time < (smeared_time + 50.0):
                    square_wf[i] += self._square_pulse(current_time, smeared_time, smeared_gain)

        # Digitize (no noise added — matches C++ squareWaveform)
        adc_waveform = [0] * self.daq_params.num_samples
        for i in range(self.daq_params.num_samples):
            adc_count = int(math.floor(square_wf[i] / self.daq_params.v_lsb) + self.daq_params.baseline)
            adc_count = max(0, min(adc_count, self.daq_params.max_adc))
            adc_waveform[i] = adc_count

        return adc_waveform

    # ----------------
    # Triangle Waveform
    # ----------------
    def generate_triangle_waveform(self) -> list[int]:
        """Generate digitized FADC waveform using triangle pulse shape."""
        hits = self._generate_pe_times()   # list of (hit_time, n_pe)
        square_wf = [0.0] * self.daq_params.num_samples

        # Convolution with triangle pulse template
        for hit_time, n_pe in hits:
            smeared_time = hit_time          # no TTS for triangle (deterministic shape)
            smeared_gain = self.det_params.mean_gain * n_pe
            if smeared_gain < 0.0:
                smeared_gain = 0.0

            for i in range(self.daq_params.num_samples):
                current_time = i * self.daq_params.sample_rate
                # Optimization: only evaluate near the pulse
                if (smeared_time - 5.0) < current_time < (smeared_time + 50.0):
                    square_wf[i] += self._triangle_pulse(current_time, smeared_time, smeared_gain)

        # Digitize (no noise added — matches C++ squareWaveform)
        adc_waveform = [0] * self.daq_params.num_samples
        for i in range(self.daq_params.num_samples):
            adc_count = int(math.floor(square_wf[i] / self.daq_params.v_lsb) + self.daq_params.baseline)
            adc_count = max(0, min(adc_count, self.daq_params.max_adc))
            adc_waveform[i] = adc_count

        return adc_waveform

def main():
    # Setup parameters
    det_params = DetectorParams(ap_prob=0.05, mean_gain=6.0, pulse_width=5.0)
    daq_params = DigitizerParams(num_samples=256)

    generator = WaveformGenerator(det_params=det_params, daq_params=daq_params)

    spe_waveform = generator.generate_pe_waveform()
    square_waveform = generator.generate_square_waveform()
    triangle_waveform = generator.generate_triangle_waveform()

    # Time axis in ns
    time_axis = [i * daq_params.sample_rate for i in range(daq_params.num_samples)]

    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

    # SPE Pulse Waveform
    ax1.plot(time_axis, spe_waveform, marker='o', markersize=3, linestyle='-',
             color='b', label='SPE Pulse (exp rise/decay)')
    ax1.set_title("SPE Waveform: Exponential Rise + Decay")
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("ADC Counts")
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()

    # Square Pulse Waveform
    ax2.plot(time_axis, square_waveform, marker='o', markersize=3, linestyle='-',
             color='r', label='Square Pulse')
    ax2.set_title("Square Waveform: Flat-Top Pulse Shape")
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("ADC Counts")
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()

    # Triangle Pulse Waveform
    ax3.plot(time_axis, triangle_waveform, marker='o', markersize=3, linestyle='-',
             color='g', label='Triangle Pulse')
    ax3.set_title("Triangle Waveform")
    ax3.set_xlabel("Time (ns)")
    ax3.set_ylabel("ADC Counts")
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
