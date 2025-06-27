function sero_sampling_plot_sampling_stats_v3(tr, b, w)
%% SRP_PLOT_SAMPLING_STATS_V3   Visualize TR & b-value sampling statistics
%   srp_plot_sampling_stats_v3(tr, w, b) makes a 3×2 figure showing:
%     (1) shot×slice log-TR image
%     (2) sampling density vs. slice (line plot)
%     (3) TR distribution heatmap
%     (4) b-value distribution heatmap
%
%   Inputs:
%     tr – [M×Z] matrix of TRs (zeros = no sample)
%     w  – [M×Z] weight mask (nonzero = valid sample)
%     b  – [M×Z] corresponding b-values

%% 1) Data preparation
tr(tr==0) = nan;               % avoid –Inf in log()
[M, Z]    = size(tr);
n_plot    = min(150, M);       % how many shots to display

% build an integer "slice index" matrix
sliceIdx = repmat(1:Z, M, 1);

% **use the same mask everywhere**: wherever tr>0 is your real sample
valid    = (tr>0);

% sampling density (#samples per slice, in tr>0 sense)
dens     = sum(valid, 1);

%% 2) Plotting
figure

% 2.1) shot×slice log-TR image
subplot(1,2,1)
  imagesc(abs(log(tr(1:n_plot,:))))
  colormap(flip(fix_cmap_blackredwhite))
  cb = colorbar('eastoutside');
    cb.Title.String   = 'TR [s]';
    cb.Title.FontSize = 12;
  fix_axis
  ylabel('Shot index [1]')
  xlabel('Slice position [mm]')

% 2.2) sampling density line (dummy colorbar for alignment)
subplot(3,2,2)
  plot(dens(:), 'LineWidth',2, 'Color',[139 0 0]/255)
  axis tight
  fix_axis
  ylabel('Sample dens. (1/mm)')
  cb = colorbar('eastoutside');
    cb.Visible = 'off';

% 2.3) TR distribution heatmap
subplot(3,2,4)
  px = sliceIdx(valid);   % slice indices for each valid shot
  py = tr(valid);         % corresponding TR values
  histogram2(px, py, 1:(Z+1), linspace(0,5,16), ...
             'FaceColor','flat','FaceAlpha',1,'LineStyle','none')
  view(0,90); set(gca,'Layer','top')
  colormap(flip(fix_cmap_blackredwhite,1).^2)
  cb = colorbar('eastoutside');
    cb.Title.String   = 'Counts';
    cb.Title.FontSize = 12;
  fix_axis
  ylabel('TR dist. [s]')
  % xlabel('Slice position [mm]')

% 2.4) b-value distribution heatmap
subplot(3,2,6)
  pyb   = (tr>0).*b;
  pyb   = pyb(w>0);
  histogram2(px, pyb, 1:(Z+1), linspace(0,1.4,11), ...
             'FaceColor','flat','FaceAlpha',1,'LineStyle','none')
  view(0,90); set(gca,'Layer','top')
  colormap(flip(fix_cmap_blackredwhite,1).^2)
  cb = colorbar('eastoutside');
    cb.Title.String   = 'Counts';
    cb.Title.FontSize = 12;
  fix_axis
  set(gca,'YTick',[0.1 0.5 0.9 1.4])
  ylabel('b-val dist. [ms/μm^2]')
  xlabel('Slice position [mm]')

%% 3) Final formatting & tick scaling
fix_figure



% multiply every x-axis tick by 1.5
allAxes = findall(gcf,'Type','axes');
for ax = allAxes'
  % grab the current tick positions
  xt = ax.XTick;

  % pick up to 4 evenly spaced ticks
  if numel(xt) > 4
    idx = round(linspace(1, numel(xt), 4));
    xt = xt(idx);
  end

  % reassign those as the only ticks
  ax.XTick = xt;

  % now relabel them by 1.5×
  ax.XTickLabel = arrayfun(@(v) sprintf('%.1f', v*1.5), ...
                           xt, 'UniformOutput', false);
end

end
