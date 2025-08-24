P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1;
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1;
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1;
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1;
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0;
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0;
    ];
blockSize = 27;
H = ldpcQuasiCyclicMatrix(blockSize,P);
[n, m] = size(H);
cfgLDPCEnc = ldpcEncoderConfig(H);
cfgLDPCDec_l = ldpcDecoderConfig(H, "layered-bp");
cfgLDPCDec = ldpcDecoderConfig(H);


M = 2;
maxnumiter = 5;
snr = linspace(0, 10, 10);
numframes = 100000;
data = randi([0 1], cfgLDPCEnc.NumInformationBits, numframes);
ndata = cfgLDPCDec.NumInformationBits * numframes;
encodedData = ldpcEncode(data, cfgLDPCEnc);

modSignal = pskmod(encodedData, M);

ber_flood = zeros(1, size(snr, 2));
ber_layered = zeros(1, size(snr, 2));

for i = 1 : size(snr, 2)
    [rxsig, noisevar] = awgn(modSignal, snr(i));
    
    llrout = pskdemod(rxsig, M, 'OutputType', 'llr', 'NoiseVariance', noisevar);
    
    rxbits_flood = ldpcDecode(llrout, cfgLDPCDec, maxnumiter);
    rxbits_layered = ldpcDecode(llrout, cfgLDPCDec_l, maxnumiter);
    
    ber_flood(i) = nnz(data ~= rxbits_flood) / ndata;
    ber_layered(i) = nnz(data ~= rxbits_layered) / ndata;
end
 
figure;
semilogy(snr, ber_flood, 'b-o');      
hold on;
semilogy(snr, ber_layered, 'r-s');
hold off;
xlabel('SNR (dB)');
ylabel('BER');
legend('flood', 'layered');
grid on;
