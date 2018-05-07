% sample script to run MakeModRaster

%% set the variables
% import sampleData.mat, this contains Chan1_ds and Chan3_ds data.
dam.signal.NAcc = Chan1_ds;
dam.signal.PFC = Chan3_ds;
startLastRangeSamples = [1 numel(Chan1_ds)];
dataAcq = 'TDT';


%% compute modulation index (Phase amplitude coupling)
chanNameCellPhaseFreq = {'PFC','NAcc'};
chanNameCellAmpFreq = {'NAcc','PFC'};
[modStruct,rasterWindowTimesSamplesStruct,flow,fhigh] = MakeModRaster(dam,startLastRangeSamples,chanNameCellPhaseFreq,chanNameCellAmpFreq,dataAcq,'specifyBehavs',0);

%% compute coherence
chanNameCell1 = {'NAcc'};
chanNameCell2 = {'PFC'};

[coherenceStruct,S1Struct,S2Struct,rasterWindowTimesSamplesStruct,fSpectOut] = MakeCoherenceRaster(dam,startLastRangeSamples,chanNameCell1,chanNameCell2,dataAcq,'specifyBehavs',0)