// simpleSignalProcessing.js

export function computeDFT(signal, sampleRate) {
	let N = signal.length;
	let real = new Array(N).fill(0);
	let imag = new Array(N).fill(0);
	
	// Handle odd N for frequency array indexing
	let halfN = Math.floor(N / 2);
	let frequencies = Array.from({ length: halfN + 1 }, (_, i) => (i * sampleRate) / N);

	for (let k = 0; k < N; k++) {
		for (let n = 0; n < N; n++) {
			let angle = (-2 * Math.PI * k * n) / N;
			real[k] += signal[n] * Math.cos(angle);
			imag[k] += signal[n] * Math.sin(angle);
		}
	}

	let psd = real.map((re, i) => 2/N*Math.sqrt(re * re + imag[i] * imag[i]));
	let phase = real.map((re, i) => Math.atan2(imag[i], re));

	return { frequencies, psd, phase, real, imag };
}

// Compute Inverse DFT (IDFT)
function analyticSignal(signal) {
	let N = signal.length;
	let { real, imag } = computeDFT(signal, 2*Math.PI); // outputs : { frequencies, psd, phase, real, imag };

	// Apply Hilbert transform filter in frequency domain
	let h = new Array(N).fill(0);
	h[0] = 1; // DC remains unchanged
	if (N % 2 === 0) {
		h[N / 2] = 1; // Nyquist frequency (only for even N)
	}
	for (let i = 1; i < Math.floor(N / 2); i++) {
		h[i] = 2; // Double positive frequencies
	}

	// Multiply FFT result by the Hilbert filter
	for (let i = 0; i < N; i++) {
		real[i] *= h[i];
		imag[i] *= h[i];
	}

	return inverseDFT(real, imag);
}

export function inverseDFT(real, imag) {
	let N = real.length;
	// Compute the inverse DFT (IDFT) to get the analytic signal
	let signalReal = new Array(N).fill(0);
	let signalImag = new Array(N).fill(0);

	for (let n = 0; n < N; n++) {
		for (let k = 0; k < N; k++) {
			let angle = (2 * Math.PI * k * n) / N;
			signalReal[n] += real[k] * Math.cos(angle) - imag[k] * Math.sin(angle);
			signalImag[n] += real[k] * Math.sin(angle) + imag[k] * Math.cos(angle);
		}
		signalReal[n] /= N; // Normalize
		signalImag[n] /= N; // Normalize
	}

	// Return the analytic signal as an array of [real, imag] components
	return signalReal.map((re, i) => (signalImag === 0) ? re : [re, signalImag[i]]);
}

export function generateSineWave(time, frequency=10, phase=0, noise=0, envFreq=0, envPhi=0) {
	if (envFreq === 0){
		return time.map((t, index) => Math.cos(2 * Math.PI * frequency * t + phase) + noise[index]);
	} else {
		const env = time.map((t, index) => .5*(.8*Math.cos(2 * Math.PI * envFreq * t + envPhi) + 1));
		return time.map((t, index) => env[index]*Math.cos(2 * Math.PI * frequency * t + phase) + noise[index]);
	}
}

export function pearsonCorrelation(x, y) {
	let n = x.length;
	let sumX = x.reduce((a, b) => a + b, 0);
	let sumY = y.reduce((a, b) => a + b, 0);
	let sumXY = x.map((xi, i) => xi * y[i]).reduce((a, b) => a + b, 0);
	let sumX2 = x.map(xi => xi * xi).reduce((a, b) => a + b, 0);
	let sumY2 = y.map(yi => yi * yi).reduce((a, b) => a + b, 0);

	let numerator = (n * sumXY) - (sumX * sumY);
	let denominator = Math.sqrt((n * sumX2 - sumX ** 2) * (n * sumY2 - sumY ** 2));
	let correlation = (denominator === 0) ? 0 : (numerator / denominator);
	
	// Compute Linear Regression inside the same function
    let slope = numerator / (n * sumX2 - sumX ** 2);
    let intercept = (sumY - slope * sumX) / n;

    return { correlation, slope, intercept };
}

export function computeEnvCorr(analytic1, analytic2) {
	// compute envelop
	const envp1 = analytic1.map(([real, imag]) => Math.sqrt(real ** 2 + imag ** 2));
			pearsonCorrelation(x, y)
	const envp2 = analytic2.map(([real, imag]) => Math.sqrt(real ** 2 + imag ** 2));
			pearsonCorrelation(x, y)

	return pearsonCorrelation(envp1, envp2)
}

export function computePLV(phase1, phase2) {
	let N = phase1.length;
	if (N !== phase2.length) {
		throw new Error("Phase arrays must have the same length.");
	}

	let sumReal = 0;
	let sumImag = 0;

	for (let i = 0; i < N; i++) {
		let deltaPhi = phase1[i] - phase2[i];
		sumReal += Math.cos(deltaPhi);
		sumImag += Math.sin(deltaPhi);
	}

	let magnitude = Math.sqrt(sumReal * sumReal + sumImag * sumImag) / N; // PLV magnitude
	let phaseLocking = Math.atan2(sumImag, sumReal); // PLV phase locking
	return { magnitude, phaseLocking};
}

export function computeiPLV(PLV) {

	let magnitude = Math.abs(PLV.magnitude * Math.sin(PLV.phaseLocking)); // iPLV magnitude
	let phaseLocking = Math.sign(Math.sin(PLV.phaseLocking)) * Math.PI/2; // PLV phase locking
	return { magnitude, phaseLocking};
}

export function computeciPLV(PLV) {
	let numerator = Math.abs(PLV.magnitude * Math.sin(PLV.phaseLocking));
	let denominator = Math.sqrt(1 -  Math.pow(PLV.magnitude * Math.cos(PLV.phaseLocking),2));
	let magnitude = (denominator === 0) ? 0 : (numerator/denominator); // iPLV magnitude
	let phaseLocking = Math.sign(Math.sin(PLV.phaseLocking)) * Math.PI/2; // PLV phase locking
	return { magnitude, phaseLocking};
}

export function computeSampleCoh(analytic1, analytic2) {
	// compute power
	const power1 = zWave1.reduce((sum, [re, im]) => sum + (re ** 2 + im ** 2), 0); // First sine wave
	const power2 = zWave2.reduce((sum, [re, im]) => sum + (re ** 2 + im ** 2), 0); // Second sine wave
	
	let sampleMagn = zWave1.map(([real, imag]) => N*(real ** 2 + imag ** 2)/Math.sqrt(power1 * power2));				
	let samplePhase = zWave1.map(([real, imag], i) => Math.atan2(imag * zWave2[i][0] - real * zWave2[i][1], real * zWave2[i][0] + imag * zWave2[i][1] ));
	
	return {sampleMagn, samplePhase};
}

export function computeCoherence(analytic1, analytic2) {
	if (analytic1.length !== analytic2.length) {
		throw new Error("Signals must have the same length.");
	}

	let N = analytic1.length;
	let coh = 0;
	let sumCrossRe = 0;
	let sumCrossIm = 0;
	let sumPower1 = 0;
	let sumPower2 = 0;
	
	for (let t = 0; t < N; t++) {
	
		// Analytic signal components
		let re1 = analytic1[t][0], im1 = analytic1[t][1];
		let re2 = analytic2[t][0], im2 = analytic2[t][1];

		// Compute cross-spectral real and imaginary parts
		sumCrossRe += re1 * re2 + im1 * im2;  // Real part of A1 * conj(A2)
		sumCrossIm += im1 * re2 - re1 * im2;  // Imaginary part of A1 * conj(A2)
		
		sumPower1 += re1 ** 2 + im1 ** 2;
		sumPower2 += re2 ** 2 + im2 ** 2;

	}
	
	// Compute coherence magnitude
	let crossMag = Math.sqrt(sumCrossRe ** 2 + sumCrossIm ** 2);
	let magnitude = crossMag / Math.sqrt(sumPower1 * sumPower2);
	// Compute preferred phase difference
	let preferredPhase = Math.atan2(sumCrossIm, sumCrossRe);

	return {magnitude, preferredPhase};
}

export function computeImCoh(Coh) {

	let magnitude = Math.abs(Coh.magnitude * Math.sin(Coh.preferredPhase)); // imCoh magnitude
	let preferredPhase = Math.sign(Math.sin(Coh.preferredPhase)) * Math.PI/2; // imCoh phase
	return { magnitude, preferredPhase};
}

export function computecImCoh(Coh) {
	let numerator = Math.abs(Coh.magnitude * Math.sin(Coh.preferredPhase));
	let denominator = Math.sqrt(1 -  Math.pow(Coh.magnitude * Math.cos(Coh.preferredPhase),2));
	let magnitude = (denominator === 0) ? 0 : (numerator/denominator); // cimCoh magnitude
	let preferredPhase = Math.sign(Math.sin(Coh.preferredPhase)) * Math.PI/2; // cimCoh phase locking
	return { magnitude, preferredPhase};
}

export function generateGaussianNoise(noiseLvl, length){
	const noise = [];
	for (let i = 0; i < length; i++) {
		const u1 = Math.random();
		const u2 = Math.random();
		const z0 = noiseLvl*Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2); // Generate standard normal variable
		noise.push(z0); // Scale and shift to desired mean and standard deviation
	}
	return noise
}
