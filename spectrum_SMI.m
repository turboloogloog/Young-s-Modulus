function spectrum_SMI(result_g,n,fs)
    figure
    NFFT = 2^nextpow2(n); % Next power of 2 from length of y
    Y = fft(result_g,NFFT*20)/n;
    f = fs/2*linspace(0,1,NFFT*20/2);

    try_01 = ceil(NFFT/2*0.3); 

    x = f(1:try_01);
    y = 2*abs(Y(1:try_01));

    plot(x,y)

    title('Spectrum of g')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude of Fourier transformation')

end