function sketch_SMI(result_g,t)
figure
    plot(t,result_g);
    title('Normalized self-mixing signal')
    xlabel('Times (s)')
    ylabel('Amplitude of G(t)')
    set(gca,'xlim',[0 max(t)]);
    axis([0 max(t) -1.5 1.5]);
    grid on;
end