function simulate_noise_sweep(n_noise, n_mc, n_param, n_mc_param)
% function simulate_noise_sweep(n_noise, n_mc, n_param, n_mc_param)
%   simulate_noise_sweep                 % full run (30 noise, 50 MC, 20 param, 10 MC)
%   simulate_noise_sweep(3,3,5,2)        % quick validation run
%
% Fig 9  = AMARES vs integration: kinetic parameter and % bias vs noise level,
%          for 13C HP (kPB) and 2H DMI (Vmax).
% Fig 10 = metabolite-ratio vs true kinetic parameter (Pearson r), at low noise.
%
% Uses the same generate_fake_fids / quantify_spectra_pipline / fit_* functions
% and the same true parameters and fit settings as end_to_end_simulation_example.m.

if nargin<1||isempty(n_noise),    n_noise=30;  end
if nargin<2||isempty(n_mc),       n_mc=50;     end
if nargin<3||isempty(n_param),    n_param=20;  end
if nargin<4||isempty(n_mc_param), n_mc_param=10; end

here = fileparts(mfilename('fullpath'));
tag  = sprintf('%dx%d', n_noise, n_mc);
diary(fullfile(here, sprintf('noise_sweep_%s.out', tag))); diary on;
LWF = 5;
try
    base = fileparts(fileparts(here));
    addpath(genpath(fullfile(base,'OXSA','main')));
    addpath(genpath(fullfile(base,'end_to_end_kinetic_simulation','hyperpolarized-mri-toolbox')));
    cd(fullfile(base,'end_to_end_kinetic_simulation'));

    % ---- true parameters (match the driver) ----
    dmiT = struct('a1',17.0,'kp',0.0416,'v',0.72,'kg',0.141,'Vmax',0.985,'kl',0.124);
    kpbT = 0.05; r1pyr=1/30; r1bic=1/20; pyr0=1.0; flip=15;
    noise = logspace(-3,1,n_noise);

    %% ===================== FIG 9: noise sweep =====================
    fprintf('\n===== FIG 9 noise sweep: %d levels x %d MC =====\n', n_noise, n_mc);
    kpb_a=nan(n_noise,n_mc); kpb_i=nan(n_noise,n_mc);
    vmx_a=nan(n_noise,n_mc); vmx_i=nan(n_noise,n_mc);
    for ni=1:n_noise
        fprintf('  noise %d/%d (%.4g)\n', ni, n_noise, noise(ni));
        for mc=1:n_mc
            % 13C HP
            try
                f = generate_fake_fids('13C','noise_level',noise(ni),'k_true',kpbT, ...
                    'r1_pyr',r1pyr,'r1_bic',r1bic,'pyr_0',pyr0,'FlipAngle',flip,'n_timepoints',60);
                [ka,~,ki,~] = fit_hp_both(f, LWF);
                kpb_a(ni,mc)=ka; kpb_i(ni,mc)=ki;
            catch e, fprintf('   13C n%d m%d: %s\n',ni,mc,e.message); end
            % 2H DMI
            try
                f = generate_fake_fids('2H','noise_level',noise(ni),'a1',dmiT.a1,'kp',dmiT.kp, ...
                    'v',dmiT.v,'kg',dmiT.kg,'Vmax',dmiT.Vmax,'kl',dmiT.kl);
                [va,~,vi,~] = fit_dmi_both(f, LWF);
                vmx_a(ni,mc)=va; vmx_i(ni,mc)=vi;
            catch e, fprintf('   2H n%d m%d: %s\n',ni,mc,e.message); end
        end
    end
    f9 = plot_fig9(noise, kpb_a, kpb_i, kpbT, vmx_a, vmx_i, dmiT.Vmax);
    savepdf(f9, fullfile(here, sprintf('noise_sweep_preview_%s.pdf', tag)));

    %% ===================== FIG 10: ratio vs true parameter =====================
    fprintf('\n===== FIG 10 correlation: %d params x %d MC (noise=1e-3) =====\n', n_param, n_mc_param);
    kpbR = linspace(kpbT*0.3, kpbT*1.7, n_param);
    vmxR = linspace(dmiT.Vmax*0.3, dmiT.Vmax*1.7, n_param);
    hp_true=[]; hp_ra=[]; hp_ri=[]; dm_true=[]; dm_ra=[]; dm_ri=[];
    for j=1:n_param
        for mc=1:n_mc_param
            try
                f = generate_fake_fids('13C','noise_level',1e-3,'k_true',kpbR(j), ...
                    'r1_pyr',r1pyr,'r1_bic',r1bic,'pyr_0',pyr0,'FlipAngle',flip,'n_timepoints',60);
                evalc('qa=quantify_spectra_pipline(f,''amares'');');
                evalc('qi=quantify_spectra_pipline(f,''integration'',''window_mode'',''linewidth'',''linewidth_factor'',LWF,''baseline_poly_order'',0);');
                hp_true(end+1)=kpbR(j); hp_ra(end+1)=qa.ratios_of_interest(2); hp_ri(end+1)=qi.ratios_of_interest(2);
            catch e, fprintf('   F10 13C p%d m%d: %s\n',j,mc,e.message); end
            try
                f = generate_fake_fids('2H','noise_level',1e-3,'a1',dmiT.a1,'kp',dmiT.kp, ...
                    'v',dmiT.v,'kg',dmiT.kg,'Vmax',vmxR(j),'kl',dmiT.kl);
                evalc('qa=quantify_spectra_pipline(f,''amares'');');
                evalc('qi=quantify_spectra_pipline(f,''integration'',''window_mode'',''linewidth'',''linewidth_factor'',LWF,''baseline_poly_order'',0);');
                dm_true(end+1)=vmxR(j); dm_ra(end+1)=qa.ratios_of_interest(2); dm_ri(end+1)=qi.ratios_of_interest(2);
            catch e, fprintf('   F10 2H p%d m%d: %s\n',j,mc,e.message); end
        end
    end
    f10 = plot_fig10(hp_true,hp_ra,hp_ri,'^{13}C HP: BIC/PYR vs k_{PB}', ...
                     dm_true,dm_ra,dm_ri,'^{2}H DMI: LAC/total vs V_{max}');
    savepdf(f10, fullfile(here, sprintf('param_sweep_preview_%s.pdf', tag)));

    save(fullfile(here, sprintf('noise_sweep_%s.mat', tag)), ...
        'noise','kpb_a','kpb_i','kpbT','vmx_a','vmx_i','dmiT', ...
        'hp_true','hp_ra','hp_ri','dm_true','dm_ra','dm_ri','LWF');
    fprintf('\nDONE (%s).\n', tag);
catch ME
    fprintf('FATAL: %s\n', ME.message);
    for k=1:numel(ME.stack), fprintf('  at %s line %d\n', ME.stack(k).name, ME.stack(k).line); end
end
diary off;
end

% ------------------------------------------------------------------ HP fit
function [kpb_a,r1b_a,kpb_i,r1b_i] = fit_hp_both(f, LWF)
kpb_a=NaN; r1b_a=NaN; kpb_i=NaN; r1b_i=NaN;    % independent fits: one failing must not lose the other
evalc('qa = quantify_spectra_pipline(f,''amares'');');   % evalc: suppress PK-struct flood -- dodgy hack to suppress verbose output with amares and speed things up...
evalc('qi = quantify_spectra_pipline(f,''integration'',''window_mode'',''linewidth'',''linewidth_factor'',LWF,''baseline_poly_order'',0);');
names=f.metabolite_names;
pI=find(strcmp(names,'1-13C Pyruvate')); bI=find(strcmp(names,'1-13C Bicarbonate')); rI=find(strcmp(names,'1-13C UreaRef'));
nf = robust_nf(qa.quantified_amplitudes(rI,:), qi.quantified_amplitudes(rI,:));  
nt=length(f.time_points); TR=f.time_points(2)-f.time_points(1);
fl=[deg2rad(15)*ones(1,nt); deg2rad(15)*ones(1,nt)];
pf=struct('R1P',1/30); pe=struct('kPL',0.01,'R1L',1/25);  % match Fig 7: start kPL=0.01, default (unbounded) kPL, no imposed lower bounds
try
    sa=zeros(1,2,nt); sa(1,1,:)=qa.quantified_amplitudes(pI,:); sa(1,2,:)=qa.quantified_amplitudes(bI,:);
    evalc('fa=fit_pyr_kinetics(sa,TR,fl,pf,pe);'); kpb_a=fa.kPL; r1b_a=fa.R1L;
catch, end
try
    si=zeros(1,2,nt); si(1,1,:)=qi.quantified_amplitudes(pI,:)*nf; si(1,2,:)=qi.quantified_amplitudes(bI,:)*nf;
    evalc('fi=fit_pyr_kinetics(si,TR,fl,pf,pe);'); kpb_i=fi.kPL; r1b_i=fi.R1L;
catch, end
end

function nf = robust_nf(a_ref, i_ref)
% scale integration onto the AMARES scale via a constant reference peak,
% using only timepoints where the reference integral is positive (the small
% reference can be baseline-zeroed under overlap); NaN if none are usable.
v = isfinite(a_ref) & isfinite(i_ref) & (i_ref > 0);
if any(v), nf = mean(a_ref(v)./i_ref(v)); else, nf = NaN; end
end

% ------------------------------------------------------------------ DMI fit
function [vmax_a,kg_a,vmax_i,kg_i] = fit_dmi_both(f, LWF)
vmax_a=NaN; kg_a=NaN; vmax_i=NaN; kg_i=NaN;
evalc('qa = quantify_spectra_pipline(f,''amares'');');
evalc('qi = quantify_spectra_pipline(f,''integration'',''window_mode'',''linewidth'',''linewidth_factor'',LWF,''baseline_poly_order'',0);');
names=f.metabolite_names;
gI=find(strcmp(names,'Glucose')); lI=find(strcmp(names,'Lactate'));
% D2O (natural-abundance water) is the large, reliable 2H concentration
% reference; CDCl3 is too small to integrate robustly under overlap. nf DOES
% matter here because the DMI Michaelis-Menten Km is fixed in absolute units.
rI=find(strcmp(names,'D2O'));
nf = robust_nf(qa.quantified_amplitudes(rI,:), qi.quantified_amplitudes(rI,:));
tp=f.time_points/60;
pf=struct('f',0.6,'Km',2.8); pe=struct('a1',15,'kp',0.05,'v',0.6,'kg',0.1,'Vmax',0.8,'kl',0.1);
% generate_fake_fids now scales glucose and lactate by a COMMON factor (both = 2
% MR-visible deuterons: [6,6'-2H2]glucose -> [3,3'-2H2]lactate, cf. De Feyter 2018 
% 10.1126/sciadv.aat7314), so signals are proportional to concentration with equal weight
% and the fixed-Km fit recovers Vmax without further calibration. 
try
    ga=qa.quantified_amplitudes(gI,:); la=qa.quantified_amplitudes(lI,:);
    evalc('fa=fit_dmi_kinetics(ga'',la'',tp,pf,pe);'); vmax_a=fa.Vmax; kg_a=fa.kg;
catch, end
try
    gi_=qi.quantified_amplitudes(gI,:)*nf; li_=qi.quantified_amplitudes(lI,:)*nf;
    evalc('fi=fit_dmi_kinetics(gi_'',li_'',tp,pf,pe);'); vmax_i=fi.Vmax; kg_i=fi.kg;
catch, end
end

% ------------------------------------------------------------------ Fig 9
function fh = plot_fig9(noise, kpb_a, kpb_i, kpbT, vmx_a, vmx_i, vmxT)
fh=figure('Position',[80 80 1180 760],'Color','w');
ci=@(X)[prctile(X,2.5,2) prctile(X,97.5,2)];
band=@(x,lo,hi,c)fill([x fliplr(x)],[lo' fliplr(hi')],c,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
A=[0.85 0.33 0.10]; I=[0.00 0.45 0.74];
% (A) kPB
subplot(2,2,1); ma=nanmean(kpb_a,2); mi=nanmean(kpb_i,2); ca=ci(kpb_a); cci=ci(kpb_i);
band(noise,ca(:,1),ca(:,2),A); hold on; band(noise,cci(:,1),cci(:,2),I);
semilogx(noise,ma,'-','Color',A,'LineWidth',2,'DisplayName','AMARES');
semilogx(noise,mi,'-','Color',I,'LineWidth',2,'DisplayName','integration');
yline(kpbT,'k--','DisplayName','true'); set(gca,'XScale','log');
xlabel('noise level'); ylabel('k_{PB} (s^{-1})'); title('(A) ^{13}C HP: k_{PB}'); legend('Location','best'); grid on;
% (B) kPB bias
subplot(2,2,2);
semilogx(noise,100*(ma-kpbT)/kpbT,'-','Color',A,'LineWidth',2,'DisplayName','AMARES'); hold on;
semilogx(noise,100*(mi-kpbT)/kpbT,'-','Color',I,'LineWidth',2,'DisplayName','integration');
yline(0,'k--'); xlabel('noise level'); ylabel('k_{PB} bias (%)'); title('(B) ^{13}C HP: bias'); legend('Location','best'); grid on;
% (C) Vmax
subplot(2,2,3); ma=nanmean(vmx_a,2); mi=nanmean(vmx_i,2); ca=ci(vmx_a); cci=ci(vmx_i);
band(noise,ca(:,1),ca(:,2),A); hold on; band(noise,cci(:,1),cci(:,2),I);
semilogx(noise,ma,'-','Color',A,'LineWidth',2,'DisplayName','AMARES');
semilogx(noise,mi,'-','Color',I,'LineWidth',2,'DisplayName','integration');
yline(vmxT,'k--','DisplayName','true'); set(gca,'XScale','log');
xlabel('noise level'); ylabel('V_{max} (mM/min)'); title('(C) ^{2}H DMI: V_{max}'); legend('Location','best'); grid on;
% (D) Vmax bias
subplot(2,2,4);
semilogx(noise,100*(ma-vmxT)/vmxT,'-','Color',A,'LineWidth',2,'DisplayName','AMARES'); hold on;
semilogx(noise,100*(mi-vmxT)/vmxT,'-','Color',I,'LineWidth',2,'DisplayName','integration');
yline(0,'k--'); xlabel('noise level'); ylabel('V_{max} bias (%)'); title('(D) ^{2}H DMI: bias'); legend('Location','best'); grid on;
sgtitle('Figure 9 (regenerated): AMARES vs integration [linewidth windows + non-negativity]','Interpreter','tex');
end

% ------------------------------------------------------------------ Fig 10
function fh = plot_fig10(xa,ra,ri,ta, xb,rb,rib,tb)
fh=figure('Position',[80 80 1150 470],'Color','w');
subplot(1,2,1); corr_panel(xa,ra,ri,ta,'true k_{PB} (s^{-1})');
subplot(1,2,2); corr_panel(xb,rb,rib,tb,'true V_{max} (mM/min)');
sgtitle('Figure 10 (regenerated): metabolite ratio vs true parameter','Interpreter','tex');
end
function corr_panel(x,ra,ri,ttl,xl)
x=x(:); ra=ra(:); ri=ri(:); ok=isfinite(x)&isfinite(ra)&isfinite(ri); x=x(ok);ra=ra(ok);ri=ri(ok);
cA=corrcoef(x,ra); rA=cA(1,2); cI=corrcoef(x,ri); rI=cI(1,2);
plot(x,ra,'o','MarkerFaceColor',[0.85 0.33 0.10],'MarkerEdgeColor','k','DisplayName',sprintf('AMARES (r=%.3f)',rA)); hold on;
plot(x,ri,'s','MarkerFaceColor',[0.00 0.45 0.74],'MarkerEdgeColor','k','DisplayName',sprintf('integration (r=%.3f)',rI));
xlabel(xl,'Interpreter','tex'); ylabel('metabolite ratio'); title(ttl,'Interpreter','tex'); legend('Location','northwest'); grid on;
end

% ------------------------------------------------------------------ export
function savepdf(fh,fname)
try
    if exist('export_fig','file')==2, export_fig(fh,fname,'-pdf','-painters','-nocrop'); fprintf('Saved (export_fig) %s\n',fname);
    elseif exist('exportgraphics','file'), exportgraphics(fh,fname,'ContentType','vector'); fprintf('Saved (exportgraphics) %s\n',fname);
    else, set(fh,'PaperPositionMode','auto'); print(fh,fname,'-dpdf','-bestfit'); fprintf('Saved (print) %s\n',fname); end
catch e, fprintf('savepdf failed %s: %s\n',fname,e.message); end
end
