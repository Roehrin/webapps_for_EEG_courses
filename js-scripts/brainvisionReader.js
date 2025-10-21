// brainVisionReader.js
export async function readBrainvisionEEGData(url, hdr, scaling = 1, doOffset = false) {
	try {
		const sampleSize = 4; // 4 bytes for a 32-bit float
		if (hdr["Binary Infos"]["BinaryFormat"] !== "IEEE_FLOAT_32"){
			throw new Error(`Only IEEE_FLOAT_32 is implemented.`);
		}
		const numberOfChannels = hdr["Common Infos"]["NumberOfChannels"];
		const Fs = 1e6 / hdr["Common Infos"]["SamplingInterval"];
		const chanInfos = hdr["Channel Infos"];
		const calib = 1; // Calibration factor

    
        // Fetch the file from the URL
        const response = await fetch(url);
        if (!response.ok) {
            throw new Error(`HTTP error! Status: ${response.status}`);
        }

        // Get the file as an ArrayBuffer
        const arrayBuffer = await response.arrayBuffer();
        const nSamples = arrayBuffer.byteLength / sampleSize / numberOfChannels;
		
		// create time vector
		let timeVec = Array.from({length: nSamples}, (_, i) => (i + 1)/Fs);
        // Create a DataView for reading binary data
        const dataView = new DataView(arrayBuffer);
		// write the data so they can be directly used by plot.ly
		// init data by pre-allocating ys as vector of zeros
		let data = [];
		// important because data are multiplex
		for (let i = 0; i < numberOfChannels; i++) {
			let dummyLine = {x:timeVec, y:Array(nSamples).fill(0), name:chanInfos["Ch"+(i+1)].label, mode: "lines", line: {color: 'rgb(0,0,0)'}}; 
			data.push(dummyLine);
		}
		for (let sampId = 0; sampId < nSamples; sampId++) {
			for (let elId = 0; elId < numberOfChannels; elId++) {
				const offset = (sampId * numberOfChannels + elId) * sampleSize;
                const value = dataView.getFloat32(offset, true); // true for little-endian
				data[elId].y[sampId] = scaling*value + doOffset*elId;
			}
		}
        console.log("EEG Data Loaded");
		return data
    } catch (error) {
        console.error("Error reading EEG file:", error);
		return null; // Return null in case of an error
    }
}
// Function to parse Brain Vision Header File when read as text
function parseBrainVisionHeader(headerText) {
	const header = {};

	// Split the file content into lines
	const lines = headerText.split(/\r?\n/);

	let currentSection = null;
	let channelInfos = ['label', 'refChannel', 'resolution', 'unit'];
	for (let line of lines) {
		// Trim the line to remove leading/trailing whitespace
		line = line.trim();

		// Skip empty lines or comments
		if (line === '' || line.startsWith(';')) {
			continue;
		}

		// Check for section headers (e.g., [Common Infos])
		const sectionMatch = line.match(/^\[(.+)\]$/);
		if (sectionMatch) {
			currentSection = sectionMatch[1];
			header[currentSection] = {}; // Initialize section as an object
			continue;
		}

		// Parse key-value pairs within a section
		if (currentSection && line.includes('=')) {
			const [key, value] = line.split('=').map(part => part.trim());
			if (currentSection === "Channel Infos"){
				const elements = value.split(',');
				let channelInf = {};
				for (let i = 0; i < elements.length; i++){
					channelInf[channelInfos[i]] = elements[i];
				}
				header[currentSection][key] = channelInf;
			} else{
				header[currentSection][key] = value;
			}
		}
	}

	return header;
}

export async function readVHDRfromURL(url) {
	const response = await fetch(url);
    if (!response.ok) {
        throw new Error(`Failed to fetch file: ${response.statusText}`);
    }
    const fileContent = await response.text();
	return  parseBrainVisionHeader(fileContent);
}

export async function readVHDRfromFile(file) {
	return new Promise((resolve, reject) => {
		const reader = new FileReader();

		reader.onload = function(e) {
			try {
				const fileContent = e.target.result;
				const parsedHeader = parseBrainVisionHeader(fileContent);
				resolve(parsedHeader); // Resolve the parsed data
			} catch (error) {
				reject(error); // Handle parsing errors
			}
		};
		reader.onerror = function(e) {
			reject(new Error("Error reading file")); // Handle file read errors
		};
		reader.readAsText(file);
	});
}