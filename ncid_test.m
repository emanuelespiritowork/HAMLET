close all
clear variables
clc
warning off

load cmap_GnBrOr
load cmap_wvl
load cmap_rErr


%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sensor directory
SNS_dir = '\\10.0.1.243\projects\2023_HyperNet-POP\1_ACTIVITIES\MATCH-UPs\S2';
% sensor file name
SNS_fname = '20250313T101847_S2C_2x2_hypstar.txt';
% sensor file delimiter, starting row and column including wavelenghts
SNS_dlm = ',';
SNS_R1 = 1;
SNS_C1 = 1;
SNS_sf = 0.0001;
SNS_tag = 'S2C';

% URI of spectral response functions (SRFs)
srf_dir = '\\10.0.1.243\nr_working\shared\000_SENSORS_&_MANUALS\EO\S2';
srf_fname = 'S2C_SRF_10b.csv';
% SRF file delimiter, starting row and column including wavelenghts
srf_dlm = ',';
srf_R1 = 1;
srf_C1 = 1;

% HYPSTAR directory
HPS_dir = '\\10.0.1.243\projects\2023_HyperNet-POP\3-DATA\PROCESSED';
% HYPSTAR date and time
HPS_date = '20250313';
HPS_hh = '10';
HPS_mm = '00';
% Hypstar VAA (ranges = [0, 360])
vaa_min = 0;
vaa_max = 360;
% Hypstar VZA (ranges = [0 30])
vza_min = 9;
vza_max = 30;

% flag to remove bands affected by atmospheric influence according to
% bands set in atm_min and atm_max
% 0 = don't remove atm bands;
% 1 = remove atm bands.
atm_flg = 1;
atm_min = [580, 1340, 1800, 2350];
atm_max = [600, 1500, 2000, 2600];


%++ PLOT parameters ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% figure visibility flag ('on' or 'off')
fig_flg = 'on';
% plot info flag
% - 0 = do not show stats info
% - 1 = show stats info
plot_info_flg = 1;

% colormap:
% - Greens (LAI);
% - RdYlBu (LNC, CNC);
% - RdYlGn (LCC, CCC);
% - parula (CV).
% - GnOr   (NRLab)
cmap = cmap_GnBrOr;
% colormap inversion (0 = no inversion; 1 = inversion)
inv_flg = 0;

% Maximum value for y-axis (only for reflectance plot)
y_max = 0.5;

% minimum and maximum wavelengths
min_wvl = 400;
max_wvl = 1700;


%% SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('---------------------- START OF SCRIPT ----------------------\n')

% set output directory for plots
dir_plot = strcat(SNS_dir, '\PLOT');
mkdir(dir_plot)


%++ OPEN HYPSTAR DATA ++++++++++++++++++++++++++++++++++++++++++++++++++++

% % get information about the nc file
% ncinfo(path)
%
% % to display nc file
% ncdisp(path)

% list the content of the HPSstar folder
pattern = sprintf('**\\HYPERNETS_L_JSIT_L2A_REF_%sT%s%s*.nc', ...
    HPS_date, HPS_hh, HPS_mm);
list = dir(fullfile(HPS_dir, pattern));

% path of data source
HPS_path = fullfile(list.folder, list.name);

fprintf('- Reading file %s... ', list.name)

% open NetCDF data source
ncid = netcdf.open(HPS_path);
HPS_info = ncinfo(HPS_path);

% read variables exisiting in nc file
HPS_wvl = ncread(HPS_path,'wavelength');
HPS_ref = ncread(HPS_path,'reflectance');
HPS_vaa = ncread(HPS_path,'viewing_azimuth_angle');
HPS_vza = ncread(HPS_path,'viewing_zenith_angle');
HPS_urnd = 0.01 * ncread(HPS_path,'u_rel_random_reflectance');
HPS_usys = 0.01 * ncread(HPS_path,'u_rel_systematic_reflectance');
HPS_unc = HPS_urnd + HPS_usys;

% select reflectances and uncertainties accordingt to
% specified VAA and VZA range
ndx = logical((vaa_min <= HPS_vaa & HPS_vaa <= vaa_max) &...
    (vza_min <= HPS_vza & HPS_vza <= vza_max));
HPS_ref = HPS_ref(ndx,:);
HPS_urnd = HPS_urnd(ndx,:);
HPS_usys = HPS_usys(ndx,:);
HPS_unc = HPS_unc(ndx,:);


%++ NDVI POLAR PLOT +++++++++++++++++++++++++++++++++++++++++++++++++++++
% [red, red_ndx] = min(abs(HPS_wvl - 690));
% [nir, nir_ndx] = min(abs(HPS_wvl - 750));
% ndvi = (HPS_ref(:,nir_ndx) - HPS_ref(:,red_ndx)) ./ ...
%     (HPS_ref(:,nir_ndx) + HPS_ref(:,red_ndx));

% f = figure('Visible', fig_flg);
% hps = polaraxes;
% polarscatter(deg2rad(HPS_vaa), HPS_vza, ...
%     160, ndvi, 'filled')
% hps.ThetaZeroLocation = 'top';
% hps.ThetaDir = 'clockwise';
% hps.RLim = [0,45];

% % set axes
% ax = f.Children;
% ax.RAxisLocation = 0;
% ax.ThetaDir = 'clockwise';

fprintf('done!\n')


%++ OPEN SENSOR DATA +++++++++++++++++++++++++++++++++++++++++++++++++++++

fprintf('- Reading %s file %s... ', SNS_tag, SNS_fname)

% path of data source
SNS_path = fullfile(SNS_dir, SNS_fname);

% read sensor file
SNS_data = readmatrix(SNS_path, 'TreatAsMissing', 'NA');

% get sensor wavelengths and reflectance spectra
SNS_wvl = SNS_data(:,SNS_C1);
SNS_ref = SNS_sf * SNS_data(:,SNS_C1+1:end)';

fprintf('done!\n')


%++ PLOT REFLECTANCE +++++++++++++++++++++++++++++++++++++++++++++++++++++

fprintf('- Drawing HPSSTAR & %s reflectance plot...', SNS_tag)

% invert colormap (if inv_flg = 1)
if inv_flg == 1
    cmap = flipud(cmap);
end

% create figure handle
scrsz = get(0,'ScreenSize');
f = figure('Visible', fig_flg, 'Renderer', 'opengl', ...
    'Position',[100 100 1280 720]);
hold on;

% yyaxis left
% stdshade(amatrix, unc, flgshd, alpha, acolor, F, smth, lwidth)
p = HYPERNETS_stdshade(HPS_ref, mean(HPS_ref.*HPS_unc,1), 1, ...
    0.25, cmap(1, :), HPS_wvl, 1, 5);
p.LineStyle = '-';
if sum(strcmp(SNS_tag, {'S2A', 'S2B', 'S2C'})) >= 1
    errorbar(SNS_wvl, mean(SNS_ref, 1), std(SNS_ref, 1), ...
        'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', cmap(end, :))
else
    p = HYPERNETS_stdshade(SNS_ref, 0, 0, ...
        0.25, cmap(end, :), SNS_wvl, 1, 5);
    p.LineStyle = '-';
end
% plot(HPS_wvl, mean(HPS_ref, 1, 'omitnan'), ...
%     'color', cmap(1, :), 'LineWidth', 5)
% plot(SNS_wvl, mean(SNS_ref, 1, 'omitnan'), ...
%     'color', cmap(end, :), 'LineWidth', 5)

% yyaxis right
% % plot random & systematic uncertainties
% plot(HPS_wvl, mean(HPS_urnd,1), '-r')
% plot(HPS_wvl, mean(HPS_usys,1), '-b')
% % end
hold off

% set axes
ax = f.Children;
ax.Box = 'on';
ax.Color = [0.94 0.94 0.94];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = [1 1 1];
ax.GridAlpha = 1;
ax.FontSize = 30;
ax.XLim(1) = min_wvl;
ax.XLim(2) = max_wvl;
ax.YLim = [0 y_max];
% yyaxis left
% ylim([0 y_max]);
% yyaxis right
% ylim([0 20]);

% set labels
title(sprintf('%sT%s%s - %s vs HPS', ...
    HPS_date, HPS_hh, HPS_mm, SNS_tag));
xlabel('Wavelentgh [nm]', 'Fontsize', 30);
% yyaxis left
ylabel('R [-]', 'Fontsize', 30);
% yyaxis right
% ylabel('rel. uncertainty R [%]', 'Fontsize', 30);

% set legend
if sum(strcmp(SNS_tag, {'S2A', 'S2B', 'S2C'})) >= 1
    lgd = legend({'HPS unc.', 'HPS mean', SNS_tag});
else
    lgd = legend({'HPS unc.', 'HPS mean', ...
        sprintf('%s stdev', SNS_tag), sprintf('%s mean', SNS_tag)});
end
legend('Location', 'EastOutside')
legend('boxoff')
lgd.Title.String = 'Legend';

% save image to output file
exportgraphics(f, fullfile(dir_plot, ...
    sprintf('%sT%s%s_%s_vs_HPS_reflectance.png', ...
    HPS_date, HPS_hh, HPS_mm, SNS_tag)), ...
    'Resolution', 300)

% close figure
close(f)

fprintf('done!\n')


%++ SPECTRAL RESAMPLING ++++++++++++++++++++++++++++++++++++++++++++++++++

% set SRF path
srf_path = fullfile(srf_dir, srf_fname);

fprintf('- Reading SRF file %s... ', srf_fname)

% read SRF file (srf_R1-1, srf_C1-1 because dlmread starts @ 0,0)
srf_data = dlmread(srf_path, srf_dlm, srf_R1-1, srf_C1-1);

% get input SRF wavelengths
% (SRF wavelengths with input spectra)
srf_wvl = srf_data(srf_R1+1:end, srf_C1);

% get output SRF wavelengths
out_wvl = srf_data(srf_R1, srf_C1+1:end);

% get SRF values
srf = srf_data(srf_R1+1:end, srf_C1+1:end);

% get the size of SRF matrix
srf_size = size(srf);

fprintf('done!\n')

fprintf('- Resampling input file... ')

if length(HPS_wvl) ~= length(srf_wvl)

    % initialize the operator to resampling SRF to common wavelengths
    oper = zeros(length(HPS_wvl), length(srf_wvl));

    % evaluate operator
    for j = 1:length(HPS_wvl)
        for k = 1:length(srf_wvl)
            oper(j,k) = 1 / abs(HPS_wvl(j)-srf_wvl(k)+0.000001)^2;
        end
        oper(j,:) = oper(j,:) / sum(oper(j,:));
    end

    % new SRF resized on input wavelengths
    isrf = oper * srf;

end

% evaluate the sum of each colmumns in SRF matrix
srf_sum = sum(isrf);

% set normalization factor matrix
srf_norm = repmat(srf_sum, length(HPS_wvl), 1);

% evaluate spectral resampling
HPS_ref = HPS_ref * (isrf ./ srf_norm);

fprintf('done!\n')

% subset HPS bands to the same spectral interval
wvl_ndx = find((min_wvl < out_wvl) & (out_wvl < max_wvl));
HPS_wvl = out_wvl(wvl_ndx);
HPS_ref = HPS_ref(:,wvl_ndx);

% subset sensor bands to the same spectral interval
wvl_ndx = find((min_wvl < SNS_wvl) & (SNS_wvl < max_wvl));
SNS_wvl = SNS_wvl(wvl_ndx);
SNS_ref = SNS_ref(:,wvl_ndx);


%** ATMOSPHERIC BANDS ****************************************************

% remove atmospheric bands
if atm_flg == 1
    atm_ndx = logical((atm_min(1) < SNS_wvl & SNS_wvl < atm_max(1)) | ...
        (atm_min(2) < SNS_wvl & SNS_wvl < atm_max(2)) | ...
        (atm_min(3) < SNS_wvl & SNS_wvl < atm_max(3)));
    atm_ndx = ~atm_ndx;

    out_wvl = out_wvl(atm_ndx);
    HPS_ref = HPS_ref(:,atm_ndx);
    SNS_ref = SNS_ref(:,atm_ndx);


end

HPS_ref = mean(HPS_ref, 1, 'omitnan');
SNS_ref = mean(SNS_ref, 1, 'omitnan');


%++ STATISTICS +++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Error
ERR = SNS_ref - HPS_ref;
% relative error
rERR = 100 * (SNS_ref - HPS_ref) ./ HPS_ref;
% Mean Error ME
ME = mean(SNS_ref - HPS_ref);
% Mean Absolute Error MAE
MAE = mean(abs(SNS_ref - HPS_ref));
% Root Mean Squared Error RMSE
RMSE = sqrt(mean((SNS_ref - HPS_ref).^2));
% normalized Root Mean Squared Error nRMSE
nRMSE = sqrt(mean((SNS_ref - HPS_ref).^2)) /...
    (max(HPS_ref) - min(HPS_ref));
% relative Root Mean Squared Error rRMSE
rRMSE = sqrt(mean((SNS_ref - HPS_ref).^2)) /...
    mean(HPS_ref);
% Correlation Coefficient R
corr = corrcoef(HPS_ref, SNS_ref);
R = corr(1,2);
% Coefficient R2
R2 = R^2;

%-- STATISTICS -----------------------------------------------------------


%++ PLOTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%** SENSOR vs HYPSTAR ****************************************************

% create plot figure
f = figure('Renderer', 'opengl', 'Position', [100 100 860 720]);
hold on
plot([0, y_max], [0, y_max], ...
    'k--', 'LineWidth', 0.5)
scatter(HPS_ref, SNS_ref, ...
    100, out_wvl, 'filled');
hold off

% set Axes properties
ax = gca;
ax.Box = 'on';
ax.Color = [0.94 0.94 0.94];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = [1 1 1];
ax.GridAlpha = 1;
ax.FontSize = 30;
ax.XLim(1) = 0;
ax.XLim(2) = y_max;
ax.YLim = ax.XLim;
ax.YTick = ax.XTick;

% add colorbar and set colorbar limits
colormap(cmap_wvl)
c = colorbar;
c.Label.String = 'wavelengths [nm]';
clim([min_wvl max_wvl]);

% set labels
title(sprintf('%sT%s%s - %s vs HPS', ...
    HPS_date, HPS_hh, HPS_mm, SNS_tag));
xlabel('Hypstar R [-]', 'Fontsize', 30);
ylabel(sprintf('%s R [-]', SNS_tag), 'Fontsize', 30);

% set plot info (equation of linear regression, MAE, R2, nRMSE)
if plot_info_flg == 1
    text(0.05*y_max, 0.95*y_max, ...
        {sprintf('R^{2} = %3.2f', R2), ...
        sprintf('MAE = %3.2f', MAE), ...
        sprintf('nRMSE = %3.2f', nRMSE)}, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 24, 'Color', 'k');
end

% save image to output file
exportgraphics(f, fullfile(dir_plot, ...
    sprintf('%sT%s%s_%s_vs_HPS_scatterplot.png', ...
    HPS_date, HPS_hh, HPS_mm, SNS_tag)), ...
    'Resolution', 300)

% close figure
close(f)


%** RELATIVE ERROR *******************************************************

% create plot figure
f = figure('Renderer', 'opengl', 'Position', [100 100 1280 720]);
hold on
plot(min_wvl:1:max_wvl, zeros(1, max_wvl-min_wvl+1), ...
    'k--', 'LineWidth', 0.5)
plot(min_wvl:1:max_wvl, -20*ones(1, max_wvl-min_wvl+1), ...
    'k--', 'LineWidth', 0.5)
plot(min_wvl:1:max_wvl, 20*ones(1, max_wvl-min_wvl+1), ...
    'k--', 'LineWidth', 0.5)
scatter(out_wvl, rERR, ...
    100, rERR, 'filled');
hold off

% set Axes properties
ax = gca;
ax.Box = 'on';
ax.Color = [0.94 0.94 0.94];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = [1 1 1];
ax.GridAlpha = 1;
ax.FontSize = 30;
ax.XLim(1) = min_wvl;
ax.XLim(2) = max_wvl;
ax.YLim(1) = -100;
ax.YLim(2) = 100;
ytickformat('percentage');

% add colorbar and set colorbar limits
colormap(cmap_rERR)
c = colorbar;
c.Label.String = 'rERR [%]';
c.Ruler.TickLabelFormat = '%g%%';
clim([-20 20]);

% set labels
title(sprintf('%sT%s%s - %s vs HPS', ...
    HPS_date, HPS_hh, HPS_mm, SNS_tag));
xlabel('wavelengths [nm]', 'Fontsize', 30);
ylabel(sprintf('rERR [%%]'), 'Fontsize', 30);

% save image to output file
exportgraphics(f, fullfile(dir_plot, ...
    sprintf('%sT%s%s_%s_vs_HPS_rel_err.png', ...
    HPS_date, HPS_hh, HPS_mm, SNS_tag)), ...
    'Resolution', 300)

% close figure
close(f)

fprintf('----------------------- END OF SCRIPT -----------------------\n')


warning on
