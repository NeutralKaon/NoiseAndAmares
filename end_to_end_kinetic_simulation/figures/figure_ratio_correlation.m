function figure_ratio_correlation(matfile)
% Fig 10: Columns 1-2 show example 13C and 2H spectra at each noise level (ppm
% axis; peaks labelled in the top row), titled with the noise level -- this
% replaces a numeric row label. Columns 3-6: parameter recovery and ratio
% surrogate for 13C and 2H. AMARES red, integration blue, mean +/- 95% CI.
% Rows are noise levels.
if nargin<1, matfile='parameter_sweep_8x5x15.mat'; end
here=fileparts(mfilename('fullpath'));
base=fileparts(fileparts(here));
addpath(genpath(fullfile(base,'OXSA','main')));
S=load(fullfile(here,matfile));
mdir=here;
cA=[0.80 0.10 0.10]; cI=[0.00 0.35 0.75];
nz=S.noise; nN=numel(nz);
% data columns 2-5: {x, Fa, Fi, isrec, xlabel, ylabel, ylim}
cols = { S.kpbR,S.hp.fit_a,S.hp.fit_i,1,'True k_{pb} (s^{-1})','Fitted k_{pb} (s^{-1})',[0.02 0.08];
         S.vmxR,S.dm.fit_a,S.dm.fit_i,1,'True V_{max} (mM/min)','Fitted V_{max} (mM/min)',[0 2];
         S.kpbR,S.hp.rat_a,S.hp.rat_i,0,'True k_{pb} (s^{-1})','Bic/Pyr',[0 1];
         S.vmxR,S.dm.rat_a,S.dm.rat_i,0,'True V_{max} (mM/min)','Lac/Total',[0 0.08] };
ctitle = {'^{13}C: k_{pb}','^{2}H: V_{max}','^{13}C: ratio','^{2}H: ratio'};

cd(fullfile(base,'end_to_end_kinetic_simulation'));
nC=6; % 2 spectrum columns + 4 data columns
rng(1); % reproducible example spectra
f=figure('Position',[10 5 1440 1320],'Color','w');
for r=1:nN
  % cols 1-2: example 13C and 2H spectra at this noise level
  fx=generate_fake_fids('13C','noise_level',nz(r),'k_true',0.05,'r1_pyr',1/30,'r1_bic',1/20,'pyr_0',1.0,'FlipAngle',15);
  examplespectrum(subplot(nN,nC,(r-1)*nC+1), fx, 2, [-14 14], ...
      sprintf('^{13}C, %.3g%% noise',100*nz(r)), r==1, r==nN, '\delta - \delta_{Pyr} (ppm)');
  fd=generate_fake_fids('2H','noise_level',nz(r),'a1',S.dmiT.a1,'kp',S.dmiT.kp, ...
      'v',S.dmiT.v,'kg',S.dmiT.kg,'Vmax',S.dmiT.Vmax,'kl',S.dmiT.kl);
  examplespectrum(subplot(nN,nC,(r-1)*nC+2), fd, 4, [-0.5 8], ...
      sprintf('^{2}H, %.3g%% noise',100*nz(r)), r==1, r==nN, '\delta (ppm)');
  % cols 3-6: data
  for c=1:4
    subplot(nN,nC,(r-1)*nC+2+c); hold on; box on; grid on;
    x=cols{c,1}(:).'; Fa=cols{c,2}(:,r,:); Fi=cols{c,3}(:,r,:); isrec=cols{c,4};
    drawcell(x,Fi,cI); drawcell(x,Fa,cA);
    if isrec, plot([min(x) max(x)],[min(x) max(x)],'k--','LineWidth',0.75); end
    xlim([min(x) max(x)]); ylim(cols{c,7});
    if r==1, title(ctitle{c},'Interpreter','tex','FontSize',9); end
    if r==nN, xlabel(cols{c,5},'Interpreter','tex','FontSize',8); else, set(gca,'XTickLabel',[]); end
    ylabel(cols{c,6},'Interpreter','tex','FontSize',7.5);
    set(gca,'FontSize',7);
  end
end
axL=axes('Position',[0.34 0.965 0.32 0.03],'Visible','off'); hold(axL,'on');
hA=plot(axL,nan,nan,'-','Color',cA,'LineWidth',2); hI=plot(axL,nan,nan,'-','Color',cI,'LineWidth',2);
hD=plot(axL,nan,nan,'k--'); legend(axL,[hA hI hD],{'AMARES','Integration','identity'},'Orientation','horizontal','Box','off');
try, exportgraphics(f, fullfile(mdir,'ratios-vs-parameter-correlation-stacked.pdf'),'ContentType','vector');
     fprintf('wrote stacked Fig 10\n'); catch e, fprintf('%s\n',e.message); end
end

function examplespectrum(ax, fs, prodRow, xl, ttl, doLabels, doXlabel, xlab)
% Example spectrum at the timepoint where the product metabolite (row prodRow of
% true_amplitudes: 13C bicarbonate, 2H lactate) is maximal. Real part, scaled to
% the maximum of the noisy spectrum. Peaks are labelled (top row) at their
% simulated chemical shifts.
axes(ax);
[~,tp]=max(fs.true_amplitudes(prodRow,:));
spec=real(specFft(fs.noisy_fids(:,tp))); clean=real(specFft(fs.clean_fids(:,tp)));
n=numel(spec); imf=fs.base_freq*1e-6;
ppm=((-n/2):(n/2-1))'/(fs.dwell_time*n)/imf;
sc=max(spec); spec=spec/sc; clean=clean/sc;
plot(ppm,spec,'k','LineWidth',0.4); set(gca,'XDir','reverse'); xlim(xl);
set(gca,'YColor','none','Box','off','TickDir','out','FontSize',7);
if doXlabel, xlabel(xlab,'Interpreter','tex','FontSize',8); end
title(ttl,'FontSize',8,'Interpreter','tex');
if doLabels
  short=containers.Map( ...
    {'1-13C Pyruvate','1-13C Bicarbonate','1-13C Lactate','1-13C PyrHydrate','1-13C Alanine','1-13C UreaRef', ...
     'D2O','Glucose','Glx','Lactate','Lipid','CDCl3'}, ...
    {'Pyr','Bic','Lac','Hyd','Ala','Urea','HDO','Glc','Glx','Lac','Lip','Ref'});
  iv=fs.prior_knowledge.values.IVCellArray;
  for k=1:numel(fs.metabolite_names)
    nm=fs.metabolite_names{k}; j=find(strcmp(iv(:,1),nm),1);
    if isempty(j) || fs.true_amplitudes(k,tp)<=0, continue; end
    cs=iv{j,2}; if cs<min(xl) || cs>max(xl), continue; end
    w=abs(ppm-cs)<0.3; yk=max([clean(w); spec(w)]);
    text(cs,min(yk,0.92)+0.06,short(nm),'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',6.5,'Clipping','off');
  end
  ylim([min(spec) 1.12]);
end
end

function drawcell(x, M3, col)
m=nanmean(M3,3); m=m(:).'; lo=prctile(M3,2.5,3); lo=lo(:).'; hi=prctile(M3,97.5,3); hi=hi(:).';
ok=isfinite(m)&isfinite(lo)&isfinite(hi);
if any(ok), fill([x(ok) fliplr(x(ok))],[lo(ok) fliplr(hi(ok))],col,'FaceAlpha',0.15,'EdgeColor','none'); end
plot(x,m,'-','Color',col,'LineWidth',1.6);
end
