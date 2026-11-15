# Non-Gaussian noise in AMARES amplitudes may bias kinetic estimates in X-nuclear MRS

This repo provides (I hope!) code and data for:

> C. Laustsen and J. J. Miller, *Non-Gaussian distributed AMARES quantified amplitudes of spectral noise may bias kinetic estimations in X-nuclear MRS.* Magnetic Resonance in Medicine (2026, in review). DOI: *pending*

## The point

Kinetic models in hyperpolarized 13C MRS and 2H deuterium metabolic imaging (DMI) are usually fitted by least squares to metabolite amplitudes from a spectral quantification algorithm such as AMARES. Least squares assumes the errors on those amplitudes are Gaussian and zero-mean. Because amplitudes fitted by AMARES are positive semidefinite, this assumption is violated. We have explored (in detail) the distribution of experimental noise-only spectra found by both jMRUI and OXSA and it does not appear to be consistent with any particular distribution, but rather a distribution of distributions (!) when selected for by AIC. 

This matters because this noise **biases kinetic rate constants, and the bias depends on SNR.** We have shown here that if we add noise drawn from these empirical distributions to known kinetic curves and fit them by least squares, the recovered k<sub>pb</sub> (13C) and V<sub>max</sub> (2H) parameters of interest are biased, and the bias is worse at low SNR. In practice, most X-nuclear spectroscopy data is low SNR. 

Also, in end-to-end simulations we added complex Gaussian noise to the FID and then quantified it. Both AMARES and frequency-domain integration were biased at low SNR, in different ways. Where the resonances are well separated (13C) the two methods agree at high SNR, but where they overlap heavily (e.g. in 2H DMI) integration underestimates V<sub>max</sub> at any SNR.

Rather frustratingly to those of us who like kinetic modelling, **temporally summed metabolite ratios ("AUC ratios") are robust surrogates** of the underlying rate constants. They are approximately linear in the rate constant, but their slope depends on the quantification method when resonances overlap. 

**Some practical suggestions:**
- Characterise the noise distribution of *your* quantification pipeline from noise-only acquisitions.
- Report spectral SNR.
- Don't assume that least-squares kinetic fits on AMARES amplitudes are unbiased at low SNR.

## Contents  

```
StatisticalAnalysis/                 R + MATLAB: distribution fitting of quantified noise
  ExampleAmps-jMRUI.csv              experimental noise-only 13C amplitudes, quantified by jMRUI
  ExampleAmps_for_upload_oxsa.csv    the same acquisitions, quantified by OXSA
  DistributionFitting.R              fits candidate distributions (by `fitdistrplus`), and produces the Cullen-Frey graph, AIC/BICs 
  noise_statistics_test.m            simulated noise-only spectra quantified by AMARES / integration
  noise_statistics_test_results.mat  its output (so the R step runs immediately)
  noise_statistics_analysis.R        distribution analysis of the simulated noise-only amplitudes

end_to_end_kinetic_simulation/       MATLAB: synthetic FIDs -> quantification -> kinetic fitting
  generate_fake_fids.m               13C HP pyruvate/bicarbonate and 2H DMI glucose/lactate FIDs
  quantify_spectra_pipline.m         OXSA-AMARES or frequency-domain integration
  noise_added_to_metabolite_amplitudes_{carbon,deuterium}.m
                                     empirical non-Gaussian noise added to amplitudes; kinetic fits
  end_to_end_simulation_example.m    end-to-end example (Fig. 3)
  fit_dmi_kinetics.m, solve_dmi_model.m  DMI kinetic model (Kreis et al. 2019)
  PK_3T_*.m                          copies of the prior-knowledge files (see below)
  figures/                           noise / parameter sweeps and the plotting scripts for Figs 9-10
  hyperpolarized-mri-toolbox/        UCSF Hyperpolarized MRI Toolbox v1.7 (kinetic fitting; own LICENSE)

OXSA/                                OXSA AMARES toolbox (Purvis et al. 2017) add-ons 
```

### AMARES prior knowledge

- **End-to-end simulations (Figs 3, 9, 10):** `quantify_spectra_pipline.m` loads the 6-peak files
  `OXSA/main/+AMARES/+priorKnowledge/PK_3T_1C_13C_6Peaks.m` and `PK_3T_2H_6Peaks.m`
  (copies are in `end_to_end_kinetic_simulation/`). Every resonance is a Lorentzian singlet;
  amplitude, chemical shift and linewidth are fitted within bounds; one zero-order phase is shared
  by all peaks (the ²H lipid peak has its own); the begin time is fixed at 0. The values are listed
  in the Methods of the paper.
- **Noise-only simulations (Fig. 8):** `noise_statistics_test.m` uses the 5-peak files
  `end_to_end_kinetic_simulation/PK_3T_1C_13C_5Peaks.m` and `PK_3T_2H_5Peaks.m`, and keeps
  signed (unclipped) integrals.
- **Experimental noise-only data (Figs 1, 2):** quantified in jMRUI/OXSA with the prior knowledge
  described in the Methods of the paper.

## Requirements

- OXSA: run `./get_oxsa.sh` once. It clones OXSA (Purvis et al. 2017) into `OXSA/`, alongside the two prior-knowledge files already there.
- MATLAB (tested on R2022a) with the Optimization Toolbox and the Statistics and Machine Learning Toolbox. The scripts add OXSA and the UCSF toolbox to the path themselves.
- R with:

  ```r
  install.packages(c("tidyverse", "fitdistrplus", "EnvStats", "VGAM", "R.matlab", "gridExtra"))
  ```

## Reproducing the figures

Several figures were assembled by hand from the component plots provided by these scripts 

| Figure | Script(s) |
| --- | --- |
| 1 | `StatisticalAnalysis/DistributionFitting.R`, reading `ExampleAmps-jMRUI.csv` |
| 2 | `StatisticalAnalysis/DistributionFitting.R`, reading `ExampleAmps_for_upload_oxsa.csv` |
| 3 (B, C) | `end_to_end_kinetic_simulation/end_to_end_simulation_example.m` (panel A hand-drawn) |
| 4 | `StatisticalAnalysis/DistributionFitting.R` |
| 5 | `end_to_end_kinetic_simulation/noise_added_to_metabolite_amplitudes_carbon.m` |
| 6 | `end_to_end_kinetic_simulation/noise_added_to_metabolite_amplitudes_deuterium.m` |
| 7 | both of the above (¹³C and ²H panels) |
| 8 | `StatisticalAnalysis/noise_statistics_analysis.R` |
| 9 | `figures/simulate_noise_sweep.m` (slow), then `figures/figure_end_to_end_variation.m` |
| 10 | `figures/simulate_parameter_sweep.m` (very slow), then `figures/figure_ratio_correlation.m` |


The summarised Monte Carlo results for Figs 9 and 10 (`figures/noise_sweep_30x50.mat`, `figures/parameter_sweep_8x5x15.mat`) are included, so the plotting scripts run straight away -- they take a few hours otherwise. To regenerate them from scratch, delete them and re-run the `simulate_*` scripts. 

`figure_ratio_correlation.m` regenerates the example ¹³C and ²H spectra in the first two columns of Fig. 10 with a fixed random seed.

## Licence and third-party code

OXSA is not redistributed here (fetch it with `get_oxsa.sh`); `end_to_end_kinetic_simulation/hyperpolarized-mri-toolbox/` is a third-party toolbox redistributed under its own terms. Please cite them if you use them. This is provided under GPLv3. 

## Contact

Jack J. Miller — jack.miller@physics.org. Please email me if I've been an idiot.
