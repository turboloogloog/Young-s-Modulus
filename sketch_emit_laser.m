function sketch_emit_laser(emit_laser,n,t)    
    figure

    plot(t(1:n),emit_laser);
    title('Phase of laser without feedback')
    xlabel('Times (s)')
    ylabel('Amplitude of Phi_0 (rad)')
    set(gca,'xlim',[0 max(t)]);
    axis([0 max(t) -max(emit_laser)*1.5 max(emit_laser)*1.5]);
    grid on;
end