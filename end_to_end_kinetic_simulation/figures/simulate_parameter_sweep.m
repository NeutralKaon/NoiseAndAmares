function simulate_parameter_sweep(nP, nN, nMC)
% Grid over (true parameter x noise level) capturing, for BOTH AMARES and
% integration: the fitted kinetic parameter (kPB / Vmax) and the metabolite
% ratio (BIC/PYR ; LAC/total). Feeds the revised Fig 10:
%   (b) fitted vs true parameter, as a function of noise level
%   (c) ratio  vs true parameter, at different noise levels
% Run: simulate_parameter_sweep            (full: 10 params x 5 noise x 8 MC)
%      simulate_parameter_sweep(3,2,2)      (quick validation)
if nargin<1||isempty(nP),  nP=10; end
if nargin<2||isempty(nN),  nN=5;  end
if nargin<3||isempty(nMC), nMC=8; end
here=fileparts(mfilename('fullpath'));
tag=sprintf('%dx%dx%d',nP,nN,nMC);
diary(fullfile(here,sprintf('parameter_sweep_%s.out',tag))); diary on;
LWF=5;
try
    base=fileparts(fileparts(here));
    addpath(genpath(fullfile(base,'OXSA','main')));
    addpath(genpath(fullfile(base,'end_to_end_kinetic_simulation','hyperpolarized-mri-toolbox')));
    cd(fullfile(base,'end_to_end_kinetic_simulation'));

    dmiT=struct('a1',17.0,'kp',0.0416,'v',0.72,'kg',0.141,'Vmax',0.985,'kl',0.124);
    kpbT=0.05; r1pyr=1/30; r1bic=1/20; pyr0=1.0; flip=15;

    noise = [0.005, 0.05, 0.2, 0.5, 1.0];  noise = noise(1:nN);  % span high-SNR -> failure (SNR~1)
    kpbR  = linspace(0.015, 0.085, nP);          % true kPB sweep (13C)
    vmxR  = linspace(0.30, dmiT.Vmax*1.7, nP);   % true Vmax sweep (2H)

    hp = init(nP,nN,nMC); hp.true=kpbR; hp.noise=noise;
    dm = init(nP,nN,nMC); dm.true=vmxR; dm.noise=noise;

    for in=1:nN
        for ip=1:nP
            for m=1:nMC
                % ---- 13C ----
                try
                    f=generate_fake_fids('13C','noise_level',noise(in),'k_true',kpbR(ip), ...
                        'r1_pyr',r1pyr,'r1_bic',r1bic,'pyr_0',pyr0,'FlipAngle',flip,'n_timepoints',60);
                    [ka,ki,ra,ri]=hp_fit_ratio(f,LWF);
                    hp.fit_a(ip,in,m)=ka; hp.fit_i(ip,in,m)=ki; hp.rat_a(ip,in,m)=ra; hp.rat_i(ip,in,m)=ri;
                catch e, fprintf('13C p%d n%d m%d: %s\n',ip,in,m,e.message); end
                % ---- 2H ----
                try
                    f=generate_fake_fids('2H','noise_level',noise(in),'a1',dmiT.a1,'kp',dmiT.kp, ...
                        'v',dmiT.v,'kg',dmiT.kg,'Vmax',vmxR(ip),'kl',dmiT.kl);
                    [va,vi,ra,ri]=dmi_fit_ratio(f,LWF);
                    dm.fit_a(ip,in,m)=va; dm.fit_i(ip,in,m)=vi; dm.rat_a(ip,in,m)=ra; dm.rat_i(ip,in,m)=ri;
                catch e, fprintf('2H p%d n%d m%d: %s\n',ip,in,m,e.message); end
            end
        end
        fprintf('noise %d/%d (%.4g) done\n', in, nN, noise(in));
    end
    save(fullfile(here,sprintf('parameter_sweep_%s.mat',tag)),'hp','dm','noise','kpbR','vmxR','kpbT','dmiT');
    fprintf('DONE %s\n',tag);
catch ME
    fprintf('FATAL: %s\n',ME.message);
    for k=1:numel(ME.stack), fprintf('  at %s line %d\n',ME.stack(k).name,ME.stack(k).line); end
end
diary off;
end

function s=init(nP,nN,nMC)
z=nan(nP,nN,nMC); s=struct('fit_a',z,'fit_i',z,'rat_a',z,'rat_i',z);
end

function [ka,ki,ra,ri]=hp_fit_ratio(f,LWF)
ka=NaN;ki=NaN;ra=NaN;ri=NaN;
evalc('qa=quantify_spectra_pipline(f,''amares'');');
evalc('qi=quantify_spectra_pipline(f,''integration'',''window_mode'',''linewidth'',''linewidth_factor'',LWF,''baseline_poly_order'',0);');
nm=f.metabolite_names; pI=find(strcmp(nm,'1-13C Pyruvate')); bI=find(strcmp(nm,'1-13C Bicarbonate')); rI=find(strcmp(nm,'1-13C UreaRef'));
nt=length(f.time_points); TR=f.time_points(2)-f.time_points(1);
fl=[deg2rad(15)*ones(1,nt); deg2rad(15)*ones(1,nt)];
pf=struct('R1P',1/30); pe=struct('kPL',0.01,'R1L',1/25);  % match Fig 7 / new Fig 9: start kPL=0.01, default (unbounded) kPL
% CRITICAL: fit_pyr_kinetics is NOT scale-invariant in kPL, so integration must be
% referenced onto the AMARES (concentration) scale via the constant urea peak before
% fitting -- otherwise the raw integration scale gives a spurious kPB (matches Fig 9).
ar=qa.quantified_amplitudes(rI,:); ir=qi.quantified_amplitudes(rI,:);
vok=isfinite(ar)&isfinite(ir)&(ir>0); nf=NaN; if any(vok), nf=mean(ar(vok)./ir(vok)); end
try, sa=zeros(1,2,nt); sa(1,1,:)=qa.quantified_amplitudes(pI,:); sa(1,2,:)=qa.quantified_amplitudes(bI,:);
     evalc('g=fit_pyr_kinetics(sa,TR,fl,pf,pe);'); ka=g.kPL; catch, end
try, si=zeros(1,2,nt); si(1,1,:)=qi.quantified_amplitudes(pI,:)*nf; si(1,2,:)=qi.quantified_amplitudes(bI,:)*nf;
     evalc('g=fit_pyr_kinetics(si,TR,fl,pf,pe);'); ki=g.kPL; catch, end
try, ra=qa.ratios_of_interest(2); ri=qi.ratios_of_interest(2); catch, end
end

function [va,vi,ra,ri]=dmi_fit_ratio(f,LWF)
va=NaN;vi=NaN;ra=NaN;ri=NaN;
evalc('qa=quantify_spectra_pipline(f,''amares'');');
evalc('qi=quantify_spectra_pipline(f,''integration'',''window_mode'',''linewidth'',''linewidth_factor'',LWF,''baseline_poly_order'',0);');
nm=f.metabolite_names; gI=find(strcmp(nm,'Glucose')); lI=find(strcmp(nm,'Lactate')); rI=find(strcmp(nm,'D2O'));
tp=f.time_points/60; pf=struct('f',0.6,'Km',2.8); pe=struct('a1',15,'kp',0.05,'v',0.6,'kg',0.1,'Vmax',0.8,'kl',0.1);
% D2O reference: rescale integration onto the AMARES (concentration) scale, since
% the fixed-Km Michaelis-Menten fit requires absolute-scale inputs (nf cancels for
% the linear HP model and for ratios, but matters here). Matches simulate_noise_sweep.
ar=qa.quantified_amplitudes(rI,:); ir=qi.quantified_amplitudes(rI,:);
vok=isfinite(ar)&isfinite(ir)&(ir>0); nf=NaN; if any(vok), nf=mean(ar(vok)./ir(vok)); end
try, ga=qa.quantified_amplitudes(gI,:); la=qa.quantified_amplitudes(lI,:);
     evalc('g=fit_dmi_kinetics(ga'',la'',tp,pf,pe);'); va=g.Vmax; catch, end
try, gi=qi.quantified_amplitudes(gI,:)*nf; li=qi.quantified_amplitudes(lI,:)*nf;
     evalc('g=fit_dmi_kinetics(gi'',li'',tp,pf,pe);'); vi=g.Vmax; catch, end
try, ra=qa.ratios_of_interest(2); ri=qi.ratios_of_interest(2); catch, end
end
