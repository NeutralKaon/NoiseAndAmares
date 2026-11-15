function outStruct = PK_3T_1C_13C_5Peaks(ppm_off,sgn)
%% .M file to assemble the bounds, priorKnowledge and initialValues structs for the matlab implementation of AMARES

%Each of B, PK and IV is a 1xN struct, where N is the number of peaks. Note
%multiplets are counted as one peak.
%The fields are as follows:
%bounds           initialValues          priorKnowledge

%peakName         peakName               peakName
%chemShift        chemShift              multiplet
%linewidth        linewidth              chemShiftDelta
%amplitude        amplitude              amplitudeRatio
%phase            phase                  G_linewidth
%chemShiftDelta                          G_amplitude
%amplitudeRatio                          G_phase
%                                        G_chemShiftDelta
%                                        refPeak

%% 3T Deuterium PK file
% version 1. Created by Liam Young. University of Oxford.
% version 2: Rolf Schulte, GEHC
% version 3: Michael Væggemose, GEHC/AU - Change pk file from 2C-13C-10Peaks to 1C-13C-5Peaks 

%% Bounds
fields.Bounds = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio'};
values.boundsCellArray = {...
    '1-13C Pyruvate',   [-10,10],           [0,120],     [0,inf],     [0,360], [],               [];
    '1-13C Bicarbonate',[-11.1,-9.1],        [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C Lactate',    [11.1,13.1],         [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C PyrHydrate', [7.2,9.2],           [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C Alanine',    [4.6,6.6],           [0,20],     [0,inf],     [0,360], [],               [];    
    };

%% initialValues
fields.IV = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase'};
values.IVCellArray = {...
    '1-13C Pyruvate',     0,         50,          10,          180;
    '1-13C Bicarbonate',  -10.1,     10,          0.1,           0;    
    '1-13C Lactate',      12.1,      10,          0.1,           0;
    '1-13C PyrHydrate',   8.2,       10,          0.1,           0;
    '1-13C Alanine',      5.6,       10,          0.1,           0;
    };

%%
fields.PK = {
    'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak'};
values.PKCellArray = {...
    '1-13C Pyruvate',   [],          [],               [],               [],            [],            1,         [],         [],                 1;
    '1-13C Bicarbonate',[],          [],               [],               [],            [],            1,         [],         [],                 0;    
    '1-13C Lactate',    [],          [],               [],               [],            [],            1,         [],         [],                 0;
    '1-13C PyrHydrate', [],          [],               [],               [],            [],            1,         [],         [],                 0;
    '1-13C Alanine',    [],          [],               [],               [],            [],            1,         [],         [],                 0;    
    };


%% centre to f0 + adjust sign
if ((abs(ppm_off)>1d-10) || (sgn<0))
    for l=1:size(values.boundsCellArray,1)
        if sgn<0
            values.boundsCellArray{l,2} = sgn*(values.boundsCellArray{l,2}([2,1]) - ppm_off);
        else
            values.boundsCellArray{l,2} = values.boundsCellArray{l,2} - ppm_off;
        end
        values.IVCellArray{l,2} = sgn*(values.IVCellArray{l,2} - ppm_off);
    end
end


%% Pass to the function which assembles the constraints into structs and saves them
outStruct = AMARES.priorKnowledge.preparePriorKnowledge(fields,values);
outStruct.svnVersion = '$Rev: 8034 $';
outStruct.svnHeader = '$Header: RFS $';
