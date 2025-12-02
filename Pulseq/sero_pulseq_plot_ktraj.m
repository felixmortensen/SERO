function sero_pulseq_plot_ktraj(seq)

[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();

% plot k-spaces
subplot(2,1,1)
plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; 
plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

subplot(2,1,2)
plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;
plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points