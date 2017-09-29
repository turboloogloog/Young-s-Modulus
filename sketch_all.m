function sketch_all(emit_laser, result_y, result_g,n,fs,t)
    figure
    NFFT = 2^nextpow2(n); % Next power of 2 from length of y
    Y = fft(result_g,NFFT*20)/n;
    f = fs/2*linspace(0,1,NFFT*20/2);

    try_01 = ceil(NFFT/2*0.3); 
    subplot(2,2,4)

    x = f(1:try_01);y = 2*abs(Y(1:try_01));

    plot(x,y)

    title('Spectrum of g')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude of Fourier transformation')


    NFFTy = 2^nextpow2(n); % Next power of 2 from length of y
    Yy = fft(emit_laser,NFFT*20)/n;
    fy = fs/2*linspace(0,1,NFFT*20/2);
    subplot(2,2,2)

    x=fy(1:try_01);y=2*abs(Yy(1:try_01));

    plot(x,y)

    title('Spectrum of phi_0')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude of Fourier transformation')


    subplot(2,2,1)
    plot(t(1:n),emit_laser);
    title('Phase of laser without feedback')
    xlabel('Times (s)')
    ylabel('Amplitude of Phi_0 (rad)')
    set(gca,'xlim',[0 max(t)]);
    axis([0 max(t) -max(emit_laser)*1.5 max(emit_laser)*1.5]);
    grid on;


    subplot(2,2,3)
    plot(t,result_g);
    title('Normalized self-mixing signal')
    xlabel('Times (s)')
    ylabel('Amplitude of G(t)')
    set(gca,'xlim',[0 max(t)]);
    axis([0 max(t) -1.5 1.5]);
    grid on;
end