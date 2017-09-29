function spectrum_emit_laser(emit_laser,fs,n)

    NFFT = 2^nextpow2(n); 
    NFFTy = 2^nextpow2(n); % Next power of 2 from length of y
    Yy = fft(emit_laser,NFFT*20)/n;
    fy = fs/2*linspace(0,1,NFFT*20/2);
    
    try_01 = ceil(NFFT/2*0.3); 
    x=fy(1:try_01);y=2*abs(Yy(1:try_01));
    
    figure
    plot(x,y)
    title('Spectrum of phi_0')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude of Fourier transformation')

end