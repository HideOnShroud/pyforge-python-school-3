import numpy as np
from scipy.io import wavfile
from scipy import signal

class SoundWaveFactory:
    SAMPLING_RATE = 44100  # 44.1 KHz
    DURATION_SECONDS = 5
    SOUND_ARRAY_LEN = SAMPLING_RATE * DURATION_SECONDS
    MAX_AMPLITUDE = 2 ** 13

    # From a list of https://en.wikipedia.org/wiki/Piano_key_frequencies
    NOTES = {
        '0': 0, 'e0': 20.60172, 'f0': 21.82676, 'f#0': 23.12465, 'g0': 24.49971, 'g#0': 25.95654, 'a0': 27.50000,
        'a#0': 29.13524, 'b0': 30.86771, 'c0': 32.70320, 'c#0': 34.64783, 'd0': 36.70810, 'd#0': 38.89087,
        'e1': 41.20344, 'f1': 43.65353, 'f#1': 46.24930, 'g1': 48.99943, 'g#1': 51.91309, 'a1': 55.00000,
        'a#1': 58.27047, 'b1': 61.73541, 'c1': 65.40639, 'c#1': 69.29566, 'd1': 73.41619, 'd#1': 77.78175,
        'e2': 82.40689, 'f2': 87.30706, 'f#2': 92.49861, 'g2': 97.99886, 'g#2': 103.8262, 'a2': 110.0000,
        'a#2': 116.5409, 'b2': 123.4708, 'c2': 130.8128, 'c#2': 138.5913, 'd2': 146.8324, 'd#2': 155.5635,
        'e3': 164.8138, 'f3': 174.6141, 'f#3': 184.9972, 'g3': 195.9977, 'g#3': 207.6523, 'a3': 220.0000,
        'a#3': 233.0819, 'b3': 246.9417, 'c3': 261.6256, 'c#3': 277.1826, 'd3': 293.6648, 'd#3': 311.1270,
        'e4': 329.6276, 'f4': 349.2282, 'f#4': 369.9944, 'g4': 391.9954, 'g#4': 415.3047, 'a4': 440.0000,
        'a#4': 466.1638, 'b4': 493.8833, 'c4': 523.2511, 'c#4': 554.3653, 'd4': 587.3295, 'd#4': 622.2540,
        'e5': 659.2551, 'f5': 698.4565, 'f#5': 739.9888, 'g5': 783.9909, 'g#5': 830.6094, 'a5': 880.0000,
        'a#5': 932.3275, 'b5': 987.7666, 'c5': 1046.502, 'c#5': 1108.731, 'd5': 1174.659, 'd#5': 1244.508,
        'e6': 1318.510, 'f6': 1396.913, 'f#6': 1479.978, 'g6': 1567.982, 'g#6': 1661.219, 'a6': 1760.000,
        'a#6': 1864.655, 'b6': 1975.533, 'c6': 2093.005, 'c#6': 2217.461, 'd6': 2349.318, 'd#6': 2489.016,
        'e7': 2637.020, 'f7': 2793.826, 'f#7': 2959.955, 'g7': 3135.963, 'g#7': 3322.438, 'a7': 3520.000,
        'a#7': 3729.310, 'b7': 3951.066, 'c7': 4186.009, 'c#7': 4434.922, 'd7': 4698.636, 'd#7': 4978.032,
    }


    def __init__(self, duration_seconds=DURATION_SECONDS, amplitude=MAX_AMPLITUDE):
        self.duration_seconds = duration_seconds
        self.amplitude = amplitude
        self.sound_array_len = self.SAMPLING_RATE * duration_seconds
        self.common_timeline = np.linspace(0, duration_seconds, num=self.sound_array_len)


    def get_sine_wave(self, frequency, timeline=None):
        """Generate a sine wave for a given frequency."""
        if timeline is None:
            timeline = self.common_timeline
        return self.amplitude * np.sin(2 * np.pi * frequency * timeline)


    def get_square_wave(self, frequency, timeline=None):
        """Generate a square wave for a given frequency."""
        if timeline is None:
            timeline = self.common_timeline
        return self.amplitude * signal.square(2 * np.pi * frequency * timeline)


    def get_triangle_wave(self, frequency, timeline=None):
        """Generate a triangular wave for a given frequency."""
        if timeline is None:
            timeline = self.common_timeline
        return self.amplitude * signal.sawtooth(2 * np.pi * frequency * timeline, 0.5)


    def create_wave(self, note="a4", wave_type="sine", name=None):
        """Create a sound wave based on the note and type (sine, triangle, square)."""
        frequency = self.NOTES.get(note, 440)  # Default to A4 if note not found
        if wave_type == "sine":
            sound_wave = self.get_sine_wave(frequency)
        elif wave_type == "square":
            sound_wave = self.get_square_wave(frequency)
        elif wave_type == "triangle":
            sound_wave = self.get_triangle_wave(frequency)
        else:
            raise ValueError("Unsupported wave type!")
        sound_wave = sound_wave.astype(np.int16)  # Convert to 16-bit PCM format
        if name is None:
            file_name = f"{note}_{wave_type}.wav".replace("#", "s")
        else:
            file_name = f"{name}.wav"
        wavfile.write(file_name, self.SAMPLING_RATE, sound_wave)
        return sound_wave


    def read_wave_from_txt(self, file_path):
        """Read a sound wave from a .txt file."""
        return np.loadtxt(file_path, dtype=np.int16)


    def save_wave(self, wave, file_name, file_type='txt'):
        """Save the wave as a .txt or .wav file."""
        if file_type == 'txt':
            np.savetxt(f"{file_name}.txt", wave)
        elif file_type == 'wav':
            wavfile.write(f"{file_name}.wav", self.SAMPLING_RATE, wave)
        else:
            raise ValueError("Unsupported file type!")


    def normalize_sound_waves(self, *waves):
        """Normalize multiple waves in both length (to the shortest) and amplitude."""
        min_length = min(len(wave) for wave in waves)
        normalized_waves = []
        for wave in waves:
            normalized_wave = np.interp(np.linspace(0, len(wave), min_length), np.arange(len(wave)), wave)
            normalized_wave = (self.amplitude / np.max(np.abs(normalized_wave))) * normalized_wave
            normalized_waves.append(normalized_wave.astype(np.int16))
        return normalized_waves


    def print_wave_details(self, wave):
        """Prints basic details about the wave."""
        print(f"Wave details:")
        print(f"  Length: {len(wave)} samples")
        print(f"  Max Amplitude: {np.max(wave)}")
        print(f"  Min Amplitude: {np.min(wave)}")
        print(f"  Mean Amplitude: {np.mean(wave)}")


    def apply_adsr(self, wave, attack=0.1, decay=0.1, sustain=0.7, release=0.1):
        """Apply ADSR to the wave."""
        length = len(wave)
        attack_len = int(attack * length)
        decay_len = int(decay * length)
        release_len = int(release * length)
        sustain_len = length - (attack_len + decay_len + release_len)
        envelope = np.concatenate([
            np.linspace(0, 1, attack_len),  # Attack
            np.linspace(1, sustain, decay_len),  # Decay
            np.full(sustain_len, sustain),  # Sustain
            np.linspace(sustain, 0, release_len)  # Release
        ])
        return (wave * envelope).astype(np.int16)


    def combine_waves(self, *waves):
        """Combine multiple waves into a sequence."""
        return np.concatenate(waves)


    def read_melody(self, melody_str):
        """Read a melody string and generate the waves."""
        melody = []
        parts = melody_str.split()
        i = 0
        while i < len(parts):
            note = parts[i]
            if note.startswith('('):  # Handle chords
                chord_notes = note[1:-1].split()  # Remove parentheses and split notes
                chord_waves = [self.create_wave(note=n, wave_type='sine') for n in chord_notes]
                combined_chord = sum(chord_waves) // len(chord_waves)  # Avg the chord waves
                melody.append(combined_chord)
                i += 1
            else:
                duration = float(parts[i + 1][:-1])  # Parse duration in seconds
                wave = self.create_wave(note=note, wave_type='sine')
                wave = wave[:int(self.SAMPLING_RATE * duration)]  # Trim to duration
                melody.append(wave)
                i += 2
        return self.combine_waves(*melody)


if __name__ == "__main__":
    factory = SoundWaveFactory()
    wave = factory.create_wave(note="a4", wave_type="sine")
    factory.print_wave_details(wave)
    adsr_wave = factory.apply_adsr(wave)
    factory.save_wave(adsr_wave, "a4_adsr", file_type='wav')
    melody = factory.read_melody("g4 0.2s b4 0.2s (g3 d5 g5) 0.5s")
    factory.save_wave(melody, "melody", file_type='wav')
