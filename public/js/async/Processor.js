"use strict";
var worker = this;

Math.logB = function(base, x) { // base, followed by number to find log of
	return Math.log(x) / Math.log(base);
};

function dB2Mag(dB) {
	return Math.pow(10, dB/20);
}

function mag2dB(mag) {
	return 20 * Math.logB(10, mag);
}

function zeroFill(arr) {
	var n = arr.length;
	while (n) {
		n--;
		arr[n] = 0;
	}
}

// takes a signal, the amount of volume gainAdj, and whether the amount is specified in magnitude or decibels
function gain(signal, result, gainAdj, isdB) {

	var length = signal.length;

	if ( length == 0 ) { 
		alert("Cannot change volume. No data.");
		return false;
	}
	
	if ( isdB ) gainAdj = dB2Mag(gainAdj);
		
	for ( var n = 0; n < length; n++) result[n] = signal[n] * gainAdj;

	console.log("Gain adjustment success.");

	return true;
}

// takes a signal and an optional decibel level to normalize to. if nothing is specified it defaults to full scale
function normalize( signal, dB ) {
	
	var mag = 1;
	var length = signal.length;
	var hVal = 0;
	var Xabs;
	
	if ( length == 0 ) { 
		alert("Cannot change volume. No data.");
		return false;
	}
	
	if (dB) mag = dB2Mag(dB); // if a value is specified for dB, change default magnitude
		
	for ( var n = 0; n < length; n++ ) {	
		Xabs = Math.abs(signal[n]);
		if (Xabs > hVal) hVal = Xabs;
	}
	
	hVal /= mag;
	
	for ( var n = 0; n < length; n++ ) signal[n] /= hVal;

	console.log("Volume adjusted by " + -mag2dB(hVal) + " dB.");

	return true;

}

// reverses audio
function reverse( signal ) {

	var length = signal.length;
	var mem;
	
	for ( var n = 0; n < length / 2; n++ ) {
		mem = signal[n];
		signal[n] = signal[length-1-n];
		signal[length-1-n] = mem;
	}
	
	console.log("Reverse success.");

	return true;

}

// rectifies audio, takes an AudioObj
function rectify( signal ) {

	console.log("Rectifying...");
	
	var length = signal.length;
	
	for ( var n = 0; n < length; n++ ) {
		signal[n] = Math.abs(signal[n]);
	}

	console.log("Done.");

	return true;
}

// creates linear fade, takes an AudioObj, sample to start fade, length of fade in samples, and a bool saying whether to fade in or fade out
function fadeLin( signal, start, length, fadeIn ) {

	var fadeRate = 1.0/length;

	for ( var n = 0; n < length; n++ ) {
		if ( fadeIn ) {
			signal[start+n] *= fadeRate * n;
		} else {
			signal[start+n] *= 1 - (fadeRate * n);
		}
	}

	return true;

}

// creates constant power fade, takes an AudioObj, sample to start fade, length of fade in samples, and a bool saying whether to fade in or fade out
function fadePwr( signal, start, length, fadeIn ) {

	var fadeRate = Math.PI/(2.0 * length);

	for ( var n = 0; n < length; n++ ) {
		if ( fadeIn ) {
			signal[start+n] *= Math.sin(fadeRate * n);
		} else {
			signal[start+n] *= Math.cos(fadeRate * n);
		}
	}

	return true;

}

// creates inverse power fade, takes an AudioObj, sample to start fade, length of fade in samples, and a bool saying whether to fade in or fade out
function fadePwrInv( signal, start, length, fadeIn ) {

	var fadeRate = Math.PI/(2.0 * length);

	for ( var n = 0; n < length; n++ ) {
		if ( fadeIn ) {
			signal[start+n] *= 1 - Math.cos(fadeRate * n);
		} else {
			signal[start+n] *= 1 - Math.sin(fadeRate * n);
		}
	}

	return true;

}


// finds the rms of a signal
function getRMS(signal) {
	
	var length = signal.length;
	
	if ( length == 0 ) { 
		alert("Cannot get RMS. No data.");
		return false;
	}

	var rms = 0;
	var rmsdB = 0;
	
	for ( var n = 0; n < length; n++) rms += Math.pow(signal[n], 2);	
	
	rms = Math.pow(rms / length, 0.5);	
	rmsdB = mag2dB(rms);

	console.log("RMS: " + rms + " (" + rmsdB + " dB).");

	return rmsdB;
	 
}

function toPolar( REX, IMX ) {

	var REXmem = 0;
	var IMXmem = 0;
	
	var length = REX.length;

	for ( var n = 0; n < length; n++ ) {

		REXmem = Math.sqrt(Math.pow(REX[n], 2) + Math.pow(IMX[n], 2));

		if ( REX[n] == 0 ) {
			IMXmem = ( IMX[n] < 0 ) ? -Math.PI/2 : Math.PI/2;
		} else {
			IMXmem = Math.atan(IMX[n]/REX[n]);
		}

		if ( REX[n] < 0 && IMX[n] < 0 ) {
			IMXmem -= Math.PI;
		} else if ( REX[n] < 0 && IMX[n] > 0 ) {
			IMXmem += Math.PI;
		}

		REX[n] = REXmem;
		IMX[n] = IMXmem;
	}

	return true;

}

function toRect( MAG, PHASE ) {

	var MAGmem = 0;
	
	var length = MAG.length;

	for ( var n = 0; n < length; n++ ) {

		MAGmem = MAG[n] * Math.cos(PHASE[n]);
		PHASE[n] =  MAG[n] * Math.sin(PHASE[n]);
		MAG[n] = MAGmem;

	}

	return true;

}

function FFTbase( REX, IMX, N ) { // DO NOT TRY TO FIND LENGTH OF REX IN HERE, Process.FFT sends N/2!

	// set constants
	var nd2 = N / 2;
	var NM1 = N - 1;
	var m = Math.floor( Math.log(N) / Math.log(2) );
	var j = nd2;
	var k = nd2;

	var tr, ti, ur, ui, sr, si;

	// bit reversal sorting
	for ( var i = 1; i < NM1; i++ ) {
		if ( i < j ) {
			tr = REX[j];
			ti = IMX[j];
			REX[j] = REX[i];
			IMX[j] = IMX[i];
			REX[i] = tr;
			IMX[i] = ti;
		}
		k = nd2;
		while ( k <= j ) {
			j -= k;
			k /= 2;
		}
		j += k;
	}

	for ( var l = 1; l <= m; l++ ) {

		var le = Math.floor(Math.pow(2, l));
		var le2 = le / 2;
		ur = 1;
		ui = 0;
		sr = Math.cos(Math.PI / le2);
		si = -Math.sin(Math.PI / le2);

		for ( var j = 1; j <= le2; j++ ) {

			var jm1 = j - 1;

			for ( var i = jm1; i < N; i+=le ) {
				var ip = i + le2;
				tr = REX[ip]*ur - IMX[ip]*ui;
				ti = REX[ip]*ui + IMX[ip]*ur;
				REX[ip] = REX[i] - tr;
				IMX[ip] = IMX[i] - ti;
				REX[i] += tr;
				IMX[i] += ti;
			}
	
			tr = ur;
			ur = tr*sr - ui*si;
			ui = tr*si + ui*sr;

		}
	}

	return;
}

function FFT( REX, IMX ) {

	var N = REX.length;
	
	var NH = N/2;
	var NM1 = N-1;
	var N4 = N/4;
	var l = Math.floor( Math.log(N) / Math.log(2) );
	var le = Math.floor(Math.pow(2, l));
	var le2 = le / 2;
	var jm1, im, ip2, ipm, ip;

	var tr, ti, ur = 1, ui = 0, sr = Math.cos(Math.PI / le2), si = -Math.sin(Math.PI / le2);

	for( var i = 0; i < NH; i++ ) {
		REX[i] = REX[2*i];
		IMX[i] = REX[2*i+1];
	}

	FFTbase(REX, IMX, NH);

	for ( var i = 1; i < N4; i++ ) {
		im = NH-i;
		ip2 = i+NH;
		ipm = im+NH;
		REX[ip2] = (IMX[i]+IMX[im])*0.5;
		REX[ipm] = REX[ip2];
		IMX[ip2] = -(REX[i]-REX[im])*0.5;
		IMX[ipm] = -IMX[ip2];
		REX[i] = (REX[i]+REX[im])*0.5;
		REX[im] = REX[i];
		IMX[i] = (IMX[i]-IMX[im])*0.5;
		IMX[im] = -IMX[i];
	}

	REX[N*3/4] = IMX[N4];
	REX[NH] = IMX[0];
	IMX[N*3/4] = 0;
	IMX[NH] = 0;
	IMX[N4] = 0;
	IMX[0] = 0;

	for ( var j = 1; j <= le2; j++ ) {

		jm1 = j - 1;

		for ( var i = jm1; i < NM1; i+=le ) {
			ip = i + le2;
			tr = REX[ip]*ur - IMX[ip]*ui;
			ti = REX[ip]*ui + IMX[ip]*ur;
			REX[ip] = REX[i] - tr;
			IMX[ip] = IMX[i] - ti;
			REX[i] += tr;
			IMX[i] += ti;
		}
	
		tr = ur;
		ur = tr*sr - ui*si;
		ui = tr*si + ui*sr;

	}

	return true;
}

function inverseFFT( REX,  IMX ) {

	var N = REX.length;

	for ( var n = N/2+1; n < N; n++ ) {
		REX[n] = REX[N-n];
		IMX[n] = -IMX[N-n];
	}

	for ( var n = 0; n < N; n++ ) REX[n] += IMX[n];

	FFT(REX, IMX);

	for ( var n = 0; n < N; n++ ) REX[n] = (REX[n]+IMX[n])/N;	

	return true;

}

function getFFT(signal) {

	console.log("Calculating FFT...");
	
	var length = signal.length;

	var N = Math.pow(2, Math.ceil(Math.logB(2, length)));

	/*
	for ( var n = 0; n < REX[0].length; n++ ) {
		for ( var chan = 0; chan < signal.channels; chan++ ) {

			REX[chan][n] *= Math.pow(N / REX[0].length, 0.5);

		}
	}
	*/
	
	var REX = new Float32Array(N);
	var IMX = new Float32Array(N);
	
	zeroFill(REX);
	zeroFill(IMX);

	REX.set(signal); // set() actually copies values, not just references them
	
	FFT(REX, IMX);
	
	console.log("FFT done.");

	return {REX: REX, IMX:IMX};
}

function getInverseFFT(REX, IMX) {

	if ( 0 == REX.length || 0 == IMX.length ) {
		alert("Inverse FFT failed. No data.");
		return false;
	}
	
	inverseFFT(REX, IMX);

	var signal = new Float32Array(REX.buffer.slice(0));

	return signal;
}

// takes an AudioObj and kernel length in samples
function createFilter( signal, M ) {

	if ( 0 == signal.length ) {
		alert("Cannot convert to filter. No data.");
		return false;
	}
	
	if ( M%2 != 0 ) M--;
	
	if ( M/2 > signal.length ) {
		alert("Cannot convert to filter. M value must be no more than 2X filter length. M = " + M + " , Filter length = " + signal.length + ".");
		return false;
	}

	var fft = getFFT(signal);

	zeroFill(fft[1]);  // fill IMX with zeros
	
	var signal = getInverseFFT(fft[0], fft[1]);

	var arr = signal.slice(-M/2);
	signal = arr.concat(signal);
	signal.length = M;

	var radians = 2 * Math.PI / M;
	var window;

	for ( var m = 0; m < M; m++ ) {

		window = (0.42 - (0.5 * Math.cos( radians * m )) + (0.08 * Math.cos( 2 * radians * m ))); //blackman window

		signal[m] *= window; 

	}

	normalize(signal);

	return signal;
	
}

// signal and result must be the same length on entering
function convolve( signal, filter, result ) {
	
	console.log("Convolving...");
	
	var sigLength = signal.length;
	var filterLength = filter.length;

	var N = Math.pow(2, Math.ceil(Math.logB(2, filterLength + 3))); // make sure N is a power of 2 and also longer than filter length + 3 (dont know why +3 though...)

	var filterREX = new Float32Array(N); // create arrays of N length to hold the fft
	var filterIMX = new Float32Array(N);
	
	zeroFill(filterREX); // fill them with zeros
	zeroFill(filterIMX);

	for ( var n = 0; n < filterLength; n++ ) {
		filterREX[n] = filter[n]; // copy filter into REX to prepare for fft (with the trailing zeros)
	}

	FFT(filterREX, filterIMX); // perform fft on filter

	var segmentSize = N - filterLength + 1; // each segment we convolve must be N - filter length + 1
	var numSegments = Math.floor(sigLength / segmentSize); // figure out how many segments there will be, not including the last partial segment

	zeroFill(result); // fill result with zeros

	//temporary buffers, they must all be the size of the filter's FFT
	var REX = new Float32Array(N);
	var IMX = new Float32Array(N);
	var TEMP = new Float32Array(N);

	var segStart;

	for ( var seg = 0; seg < numSegments + 1; seg++ ) { // +1 to allow for that 

		segStart = seg * segmentSize; // find the offset within signal that we will start at for each segment

		for ( var n = 0; n < N; n++ ) { // perhaps slightly more efficient than zeroFill() x3;
			REX[n] = 0;
			IMX[n] = 0;
			TEMP[n] = 0;
		}

		if ( seg == numSegments ) segmentSize = sigLength % segmentSize; // if this is the last segment, its probably going to be a fraction of N, so reset segmentSize

		for ( var n = 0; n < segmentSize; n++ ) REX[n] = signal[n+segStart]; // copy the section of signal into the temporary buffer for FFT

		FFT(REX, IMX); // perform FFT

		for ( var n = 0; n < N/2+1; n++ ) { // finally perform the actual convolution on the segment
			TEMP[n] = (REX[n] * filterREX[n]) - (IMX[n] * filterIMX[n]);
			IMX[n] = (IMX[n] * filterREX[n]) + (REX[n] * filterIMX[n]);
			REX[n] = TEMP[n];
		}

		inverseFFT(REX, IMX); // find inverse FFT to turn it back into a signal

		if ( seg == numSegments ) { // if this is the last segment, we don't want to copy to N because thats probably too far
			for(var n = 0; n < segmentSize; n++) {
				result[segStart + n] += REX[n];
			}
		} else {
			for(var n = 0; n < N; n++) { // copy the newly convolved segment from the temp buffer into the result buffer
				result[segStart + n] += REX[n];
			}
		}
		
		console.log((seg / numSegments) * 100);

	}

	//for ( var n = 0; n < sigLength; n++ ) signal[n] = result[n];
	

	console.log("Convolution complete. Convolved " + (numSegments+1) + " segments.");
	
	return result;

}