function figure_end_to_end_variation()
% Fig 9: for each nucleus, a half-width example spectrum and half-width
% sampling histogram on top, and a full-page-width estimate-vs-noise trace below.
% 13C first (A-C), then 2H (D-F). 
% Bias % and sigma are annotated on the histograms.
here=fileparts(mfilename('fullpath'));
base=fileparts(fileparts(here));
addpath(genpath(fullfile(base,'OXSA','main')));
addpath(genpath(fullfile(base,'end_to_end_kinetic_simulation','hyperpolarized-mri-toolbox')));
S=load(fullfile(here,'noise_sweep_30x50.mat'));
mdir=here;
A=[0.80 0.10 0.10]; I=[0.00 0.35 0.75];
nz=S.noise(:).';
[~,hi]=min(abs(nz-0.25)); hnoise=nz(hi);   % histogram noise level (25% of max peak amplitude)
exnoise=hnoise; LWF=5;   % example spectra at the same noise as the histogram (22%)

cd(fullfile(base,'end_to_end_kinetic_simulation'));
f13=generate_fake_fids('13C','noise_level',exnoise,'k_true',0.05,'r1_pyr',1/30,'r1_bic',1/20,'pyr_0',1.0,'FlipAngle',15,'n_timepoints',60);
f2h=generate_fake_fids('2H','noise_level',exnoise,'a1',17.0,'kp',0.0416,'v',0.72,'kg',0.141,'Vmax',0.985,'kl',0.124);

f=figure('Position',[20 20 1160 1180],'Color','w');
% ---- 13C (A: spectrum, B: histogram, C: full-width estimate) ----
example_spectrum(subplot(4,2,1), f13, {'1-13C Pyruvate','1-13C Bicarbonate','1-13C Lactate','1-13C PyrHydrate','1-13C UreaRef'}, ...
    LWF, '(A)  ^{13}C example spectrum', exnoise, [-14 14], 2);
histpanel(subplot(4,2,2), S.kpb_a(hi,:), S.kpb_i(hi,:), S.kpbT, 'k_{pb} (s^{-1})','^{13}C', A, I, 'B', hnoise);
estpanel(subplot(4,2,[3 4]), S.kpb_a, S.kpb_i, S.kpbT, nz, hnoise, 'k_{pb} (s^{-1})','^{13}C HP', A, I, 'C', [1 1000], [0.02 0.08]);
% ---- 2H (D: spectrum, E: histogram, F: full-width estimate) ----
example_spectrum(subplot(4,2,5), f2h, {'D2O','Glucose','Glx','Lactate'}, ...
    LWF, '(D)  ^{2}H example spectrum', exnoise, [0 6], 4);
histpanel(subplot(4,2,6), S.vmx_a(hi,:), S.vmx_i(hi,:), S.dmiT.Vmax, 'V_{max} (mM/min)','^{2}H', A, I, 'E', hnoise);
estpanel(subplot(4,2,[7 8]), S.vmx_a, S.vmx_i, S.dmiT.Vmax, nz, hnoise, 'V_{max} (mM/min)','^{2}H DMI', A, I, 'F', [1 1000], [0 3]);

try, exportgraphics(f, fullfile(mdir,'end-to-end-variation.pdf'),'ContentType','vector');
     fprintf('wrote Fig 9 (redrawn); histogram noise = %.4g\n', hnoise); catch e, fprintf('%s\n',e.message); end
end

% ---------------------------------------------------------------- sampling histogram
function histpanel(ax, va, vi, tru, xlab, tag, A, I, L, hnoise)
axes(ax); hold on; box on;
va=va(isfinite(va)); vi=vi(isfinite(vi));
edges=linspace(min([va vi tru])*0.9, max([va vi tru])*1.1, 24);
histogram(va,edges,'FaceColor',A,'FaceAlpha',0.55,'EdgeColor','none','DisplayName','AMARES'); hold on;
histogram(vi,edges,'FaceColor',I,'FaceAlpha',0.55,'EdgeColor','none','DisplayName','Integration');
xline(tru,'k--','LineWidth',1.5,'DisplayName','True Value');
ba=100*(mean(va)-tru)/tru; bi=100*(mean(vi)-tru)/tru;
text(0.035,0.96,sprintf('AMARES: %+.1f%%, \\sigma=%.2g',ba,std(va)),'Units','normalized','Color',A*0.85,'FontSize',8,'VerticalAlignment','top');
text(0.035,0.85,sprintf('Integration: %+.1f%%, \\sigma=%.2g',bi,std(vi)),'Units','normalized','Color',I*0.85,'FontSize',8,'VerticalAlignment','top');
xlabel(xlab,'Interpreter','tex'); ylabel('Count');
title(sprintf('(%c)  %s: Distribution at %.1f%% Noise',L,tag,100*hnoise),'Interpreter','tex');
legend('Location','east'); grid on; xlim([edges(1) edges(end)]);
end

% ---------------------------------------------------------------- full-width estimate
function estpanel(ax, Xa, Xi, tru, nz, hnoise, ylab, tag, A, I, L, xr, yr)
axes(ax); hold on; box on;
band=@(x,lo,hi_,c)fill([x fliplr(x)],[lo fliplr(hi_)],c,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
nzp=nz*100; hnp=hnoise*100;
ma=nanmean(Xa,2).'; mi=nanmean(Xi,2).';
cav=[prctile(Xa,2.5,2).'; prctile(Xa,97.5,2).']; civ=[prctile(Xi,2.5,2).'; prctile(Xi,97.5,2).'];
band(nzp,cav(1,:),cav(2,:),A); band(nzp,civ(1,:),civ(2,:),I);
semilogx(nzp,ma,'-','Color',A,'LineWidth',2,'DisplayName','AMARES');
semilogx(nzp,mi,'-','Color',I,'LineWidth',2,'DisplayName','Integration');
yline(tru,'k--','DisplayName','True Value'); set(gca,'XScale','log');
xline(hnp,':','Color',[.4 .4 .4],'HandleVisibility','off');
xlabel('Noise Level (% of max peak amplitude)'); ylabel(ylab,'Interpreter','tex');
title(sprintf('(%c)  %s: Estimate vs Noise',L,tag),'Interpreter','tex');
legend('Location','best'); grid on; xlim(xr); ylim(yr);
end

% ---------------------------------------------------------------- example spectrum
function example_spectrum(ax, fs, shade_names, LWF, ttl, noiselevel, xl, prod_idx)
axes(ax); hold on; box on;
[~,tp] = max(fs.true_amplitudes(prod_idx,:));  % timepoint where the product resonance peaks
fid = fs.noisy_fids(:,tp);
spec= real(specFft(fid));
n   = numel(fid); imf = fs.base_freq*1e-6;
freqAxis = ((-n/2):(n/2-1)).'/(fs.dwell_time*n);
ppm = freqAxis/imf;
spec = spec/max(spec);
iv = fs.prior_knowledge.values.IVCellArray;
pal = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.30 0.75 0.93; 0.64 0.08 0.18];
ylo = min(spec)-0.08; yhi = 1.22;
allc=cell2mat(iv(:,2)); allhw=LWF*(cell2mat(iv(:,3))/2)/imf;   % non-overlap cap (mirrors pipeline)
for k=1:numel(shade_names)
    idx = find(strcmp(iv(:,1), shade_names{k}),1);
    c = iv{idx,2}; lw = iv{idx,3}; hwhm=(lw/2)/imf;
    lo=c-LWF*hwhm; hi=c+LWF*hwhm;
    bl=find(allc<c-1e-6); if ~isempty(bl), [cb,ib]=max(allc(bl)); br=cb+allhw(bl(ib)); if br>lo, lo=(br+lo)/2; end; end
    ab=find(allc>c+1e-6); if ~isempty(ab), [ca,ia]=min(allc(ab)); al=ca-allhw(ab(ia)); if al<hi, hi=(al+hi)/2; end; end
    w = [lo hi];
    col = pal(mod(k-1,size(pal,1))+1,:);
    fill([w(1) w(2) w(2) w(1)],[ylo ylo yhi yhi], col, 'FaceAlpha',0.20,'EdgeColor',col,'LineStyle',':','LineWidth',0.5,'HandleVisibility','off');
    plot([c c],[ylo yhi-0.07],'-','Color',col,'LineWidth',0.7,'HandleVisibility','off');
    yl_txt = 1.16 - 0.075*mod(k,2);
    text(c, yl_txt, shortname(shade_names{k}),'HorizontalAlignment','center','FontSize',7.5,'Color',col*0.75,'FontWeight','bold');
end
plot(ppm, spec, 'k','LineWidth',0.5);
set(gca,'XDir','reverse'); xlim(xl); ylim([ylo yhi]);
xlabel('Chemical Shift (ppm)'); ylabel('Signal (a.u.)');
title(sprintf('%s (%.1f%% noise)', ttl, 100*noiselevel),'Interpreter','tex');
end

function s=shortname(n)
switch n
    case '1-13C Pyruvate',   s='Pyr';
    case '1-13C Bicarbonate',s='Bic';
    case '1-13C Lactate',    s='Lac';
    case '1-13C PyrHydrate', s='Pyr-H';
    case '1-13C Alanine',    s='Ala';
    case '1-13C UreaRef',    s='Urea';
    case 'D2O',              s='HDO';
    case 'Glucose',          s='Glc';
    case 'Glx',              s='Glx';
    case 'Lactate',          s='Lac';
    case 'Lipid',            s='Lip';
    case 'CDCl3',            s='CDCl_3';
    otherwise,               s=n;
end
end
