FreqEch = 44100;
signal = synSinus(1*FreqEch, 440,1,0,FreqEch);

snr = [];

for N = 0:8
    nbPaliers = 2^(N-1);
    bruit = (round(signal*nbPaliers)/nbPaliers)-signal;

    rms_bruit = sqrt(sum(bruit.^2));
    rms_signal = sqrt(sum(signal.^2));

    snr = [snr 20.0*log10(rms_signal/rms_bruit)];
end
