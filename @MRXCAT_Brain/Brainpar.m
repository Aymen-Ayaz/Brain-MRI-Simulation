function Brainpar( MRX, filename,filepath )
% This function is the parameter file for MRXCAT_Brain.
% Change parameters in section "MRXCAT settings" only
%
% Note: Not all combinations of any parameter values are possible.
%       Some parameter changes require changes in the XCAT files.
%       E.g. if you want to increase the number of segments, you need
%       more XCAT heart phases for the additional segments,
%       i.e. additional	*.bin files.


% --------------------------------------------------------------------
%   MRXCAT settings
% --------------------------------------------------------------------

% brain Rho
RhoEye = 0 ;
RhoLens = 0 ;
RhoSinus = 0 ;
RhoCorticalBone = 0.08 ;
RhoBoneMarrow = 1.0 ;
RhoSkin = 0 ;
RhoSpine = 0.91 ;
RhoSpinalCord = 0.58 ;
RhoSalivary = 0 ;
RhoBody = 0 ;
RhoCartilage = 0.08 ;
RhoMuscle = 1.0 ;
RhoGreyMatter = 0.86 ;
RhoWhiteMatter = 0.77 ;
RhoInternalCapsule = 0.77 ;
RhoAnteriorCommissure = 0.77 ;
Rhocerebellum = 0.77 ;
RhoCerebralPeduncles = 0.77 ;
RhoSuperiorCerebellarPeduncle = 0.77 ;
RhoMiddleCerebellarPeduncle = 0.77 ;
RhoCorpusCallosum = 0.77 ;
RhoHippocampus = 0.86 ;
RhoMidbrain = 0.77 ;
RhoTegmentumOfMidbrain = 0.77 ;
RhoMedulla = 0.77 ;
RhoMedullaryPyramids = 0.77 ;
RhoInferiorOlive = 0.77 ;
RhoSuperiorColliculus = 0.77 ;
RhoInferiorColliculus = 0.77 ;
RhoPeriacquaductalGreyOuter = 0.77 ;
RhoPeriacquaductalGrey = 0.77 ;
RhoSubstantiaNigra = 0.77 ;
RhoPons = 0.77 ;
RhoAmygdala = 0.86 ;
RhoFornix = 0.77 ;
RhoPutamen = 0.86 ;
RhoThalamus = 0.86 ;
RhoGlobusPallidus = 0.86 ;
RhoCaudate = 0.86 ;
RhoThirdVentricle = 1.0 ;
RhoLateralVentricle = 1.0 ;
RhoFourthVentricle = 1.0 ;
RhoBrain = 1.0 ;
RhoMamillaryBodies = 1.0 ;
RhoPinealGland = 1.0 ;
RhoCerebralAqueduct = 1.0 ;
RhoSeptumPellicudum = 0.77 ;

% brain T1
T1Eye = 0 ;
T1Lens = 0 ;
T1Sinus = 0 ;
T1CorticalBone = 269 ;
T1BoneMarrow = 350 ;
T1Skin = 0 ;
T1Spine = 515 ;
T1SpinalCord = 780 ;
T1Salivary = 0 ;
T1Body = 0 ;
T1Cartilage = 1240 ;
T1Muscle = 900 ;
T1GreyMatter = 833 ;
T1WhiteMatter = 500 ;
T1InternalCapsule = 500 ;
T1AnteriorCommissure = 500 ;
T1cerebellum = 500 ;
T1CerebralPeduncles = 500 ;
T1SuperiorCerebellarPeduncle = 500 ;
T1MiddleCerebellarPeduncle = 500 ;
T1CorpusCallosum = 725 ;
T1Hippocampus = 833 ;
T1Midbrain = 500 ;
T1TegmentumOfMidbrain = 500 ;
T1Medulla = 500 ;
T1MedullaryPyramids = 500 ;
T1InferiorOlive = 500 ;
T1SuperiorColliculus = 500 ;
T1InferiorColliculus = 500 ;
T1PeriacquaductalGreyOuter = 500 ;
T1PeriacquaductalGrey = 500 ;
T1SubstantiaNigra = 500 ;
T1Pons = 500 ;
T1Amygdala = 833 ;
T1Fornix = 500 ;
T1Putamen = 981 ;
T1Thalamus = 972 ;
T1GlobusPallidus = 746 ;
T1Caudate = 1083 ;
T1ThirdVentricle = 2569 ;
T1LateralVentricle = 2569 ;
T1FourthVentricle = 2569 ;
T1Brain = 2569 ;
T1MamillaryBodies = 2569 ;
T1PinealGland = 2569 ;
T1CerebralAqueduct = 2569 ;
T1SeptumPellicudum = 500 ;

% brain T2
T2Eye = 0 ;
T2Lens = 0 ;
T2Sinus = 0 ;
T2CorticalBone = 55 ;
T2BoneMarrow = 70 ;
T2Skin = 0 ;
T2Spine = 101 ;
T2SpinalCord = 101 ;
T2Salivary = 0 ;
T2Body = 0 ;
T2Cartilage = 34 ;
T2Muscle = 47 ;
T2GreyMatter = 83 ;
T2WhiteMatter = 70 ;
T2InternalCapsule = 70 ;
T2AnteriorCommissure = 70 ;
T2cerebellum = 70 ;
T2CerebralPeduncles = 70 ;
T2SuperiorCerebellarPeduncle = 70 ;
T2MiddleCerebellarPeduncle = 70 ;
T2CorpusCallosum = 70 ;
T2Hippocampus = 83 ;
T2Midbrain = 70 ;
T2TegmentumOfMidbrain = 70 ;
T2Medulla = 70 ;
T2MedullaryPyramids = 70 ;
T2InferiorOlive = 70 ;
T2SuperiorColliculus = 70 ;
T2InferiorColliculus = 70 ;
T2PeriacquaductalGreyOuter = 70 ;
T2PeriacquaductalGrey = 70 ;
T2SubstantiaNigra = 70 ;
T2Pons = 70 ;
T2Amygdala = 83 ;
T2Fornix = 70 ;
T2Putamen = 80 ;
T2Thalamus = 83 ;
T2GlobusPallidus = 65 ;
T2Caudate = 95 ;
T2ThirdVentricle = 329 ;
T2LateralVentricle = 329 ;
T2FourthVentricle = 329 ;
T2Brain = 329 ;
T2MamillaryBodies = 329 ;
T2PinealGland = 329 ;
T2CerebralAqueduct = 329 ;
T2SeptumPellicudum = 70 ;
%%
TR          = 400;                         % Repetition time       [ms]
TE          = 10;                          % Echo time             [ms]
Flip        = 90;                          % Flip                  [deg] modified
Frames      = 1;                             % Number of heart phases (default: 24; 0=use # XCAT frames (.bin files))
Segments    = 1;                             % Number of segments

%%
% file name to save as Nifiti
Sequence = 'Ernst';
% implemented equations for Sequence:
%               'Ernst' : Ernst Equation
%               'FID'   : The FID or S+ signal in unspoiled steady state seqs.
%               'ECHO'  : The ECHO or S- signal in unspoiled steady state seqs.
%               'SSFP'  : fully ballanced steady state free precession

Name = [date '_' Sequence ];


BoundingBox = [0,1;0,1;0,1]; % BoundingBox in rel. units
RotationXYZ = [115.0;35.0;240.0];           % Rotations about x,y,z [deg]
% x=(RL) y=(AP) z=(FH)

LowPassFilt = 1;  % low-pass filter images
FilterStr   = [0.7, 0.7, 0.7];    %[FOV,res] [resx , resy, resz]
ResizeImage = 320;

SNR         = 200 ;                        % signal-to-noise ratio 20
Coils       = 2;                            % number coils (Biot-Savart)
CoilDist    = 450;                          % body radius  450         [mm]    = distance of coil centres from origin
CoilsPerRow = 8;                            % number of coils on 1 "ring" or row of coil array (default: 8)

Trajectory  = 'Cartesian';                  % k-space trajectory (Cartesian, Radial, GoldenAngle)
Undersample = 1;                              % undersampling factor (only for Radial/GoldenAngle right now)

% --------------------------------------------------------------------
%   Read log file
% --------------------------------------------------------------------
MRX.Filename = [pathname filename];
MRX.readLogFile;

% --------------------------------------------------------------------
%   Store tissue, contrast and sequence parameters
% --------------------------------------------------------------------
MRX.Par.tissue.rhoeye = RhoEye ;
MRX.Par.tissue.rholens = RhoLens ;
MRX.Par.tissue.rhosinus = RhoSinus ;
MRX.Par.tissue.rhocorticalbone = RhoCorticalBone ;
MRX.Par.tissue.rhobonemarrow = RhoBoneMarrow ;
MRX.Par.tissue.rhoskin = RhoSkin ;
MRX.Par.tissue.rhospine = RhoSpine ;
MRX.Par.tissue.rhospinalcord = RhoSpinalCord ;
MRX.Par.tissue.rhosalivary = RhoSalivary ;
MRX.Par.tissue.rhobody = RhoBody ;
MRX.Par.tissue.rhocartilage = RhoCartilage ;
MRX.Par.tissue.rhomuscle = RhoMuscle ;
MRX.Par.tissue.rhogreymatter = RhoGreyMatter ;
MRX.Par.tissue.rhowhitematter = RhoWhiteMatter ;
MRX.Par.tissue.rhointernalcapsule = RhoInternalCapsule ;
MRX.Par.tissue.rhoanteriorcommissure = RhoAnteriorCommissure ;
MRX.Par.tissue.rhocerebellum = Rhocerebellum ;
MRX.Par.tissue.rhocerebralpeduncles = RhoCerebralPeduncles ;
MRX.Par.tissue.rhosuperiorcerebellarpeduncle = RhoSuperiorCerebellarPeduncle ;
MRX.Par.tissue.rhomiddlecerebellarpeduncle = RhoMiddleCerebellarPeduncle ;
MRX.Par.tissue.rhocorpuscallosum = RhoCorpusCallosum ;
MRX.Par.tissue.rhohippocampus = RhoHippocampus ;
MRX.Par.tissue.rhomidbrain = RhoMidbrain ;
MRX.Par.tissue.rhotegmentumofmidbrain = RhoTegmentumOfMidbrain ;
MRX.Par.tissue.rhomedulla = RhoMedulla ;
MRX.Par.tissue.rhomedullarypyramids = RhoMedullaryPyramids ;
MRX.Par.tissue.rhoinferiorolive = RhoInferiorOlive ;
MRX.Par.tissue.rhosuperiorcolliculus = RhoSuperiorColliculus ;
MRX.Par.tissue.rhoinferiorcolliculus = RhoInferiorColliculus ;
MRX.Par.tissue.rhoperiacquaductalgreyouter = RhoPeriacquaductalGreyOuter ;
MRX.Par.tissue.rhoperiacquaductalgrey = RhoPeriacquaductalGrey ;
MRX.Par.tissue.rhosubstantianigra = RhoSubstantiaNigra ;
MRX.Par.tissue.rhopons = RhoPons ;
MRX.Par.tissue.rhoamygdala = RhoAmygdala ;
MRX.Par.tissue.rhofornix = RhoFornix ;
MRX.Par.tissue.rhoputamen = RhoPutamen ;
MRX.Par.tissue.rhothalamus = RhoThalamus ;
MRX.Par.tissue.rhoglobuspallidus = RhoGlobusPallidus ;
MRX.Par.tissue.rhocaudate = RhoCaudate ;
MRX.Par.tissue.rhothirdventricle = RhoThirdVentricle ;
MRX.Par.tissue.rholateralventricle = RhoLateralVentricle ;
MRX.Par.tissue.rhofourthventricle = RhoFourthVentricle ;
MRX.Par.tissue.rhobrain = RhoBrain ;
MRX.Par.tissue.rhomamillarybodies = RhoMamillaryBodies ;
MRX.Par.tissue.rhopinealgland = RhoPinealGland ;
MRX.Par.tissue.rhocerebralaqueduct = RhoCerebralAqueduct ;
MRX.Par.tissue.rhoseptumpellicudum = RhoSeptumPellicudum ;


MRX.Par.tissue.t1eye = T1Eye ;
MRX.Par.tissue.t1lens = T1Lens ;
MRX.Par.tissue.t1sinus = T1Sinus ;
MRX.Par.tissue.t1corticalbone = T1CorticalBone ;
MRX.Par.tissue.t1bonemarrow = T1BoneMarrow ;
MRX.Par.tissue.t1skin = T1Skin ;
MRX.Par.tissue.t1spine = T1Spine ;
MRX.Par.tissue.t1spinalcord = T1SpinalCord ;
MRX.Par.tissue.t1salivary = T1Salivary ;
MRX.Par.tissue.t1body = T1Body ;
MRX.Par.tissue.t1cartilage = T1Cartilage ;
MRX.Par.tissue.t1muscle = T1Muscle ;
MRX.Par.tissue.t1greymatter = T1GreyMatter ;
MRX.Par.tissue.t1whitematter = T1WhiteMatter ;
MRX.Par.tissue.t1internalcapsule = T1InternalCapsule ;
MRX.Par.tissue.t1anteriorcommissure = T1AnteriorCommissure ;
MRX.Par.tissue.t1cerebellum = T1cerebellum ;
MRX.Par.tissue.t1cerebralpeduncles = T1CerebralPeduncles ;
MRX.Par.tissue.t1superiorcerebellarpeduncle = T1SuperiorCerebellarPeduncle ;
MRX.Par.tissue.t1middlecerebellarpeduncle = T1MiddleCerebellarPeduncle ;
MRX.Par.tissue.t1corpuscallosum = T1CorpusCallosum ;
MRX.Par.tissue.t1hippocampus = T1Hippocampus ;
MRX.Par.tissue.t1midbrain = T1Midbrain ;
MRX.Par.tissue.t1tegmentumofmidbrain = T1TegmentumOfMidbrain ;
MRX.Par.tissue.t1medulla = T1Medulla ;
MRX.Par.tissue.t1medullarypyramids = T1MedullaryPyramids ;
MRX.Par.tissue.t1inferiorolive = T1InferiorOlive ;
MRX.Par.tissue.t1superiorcolliculus = T1SuperiorColliculus ;
MRX.Par.tissue.t1inferiorcolliculus = T1InferiorColliculus ;
MRX.Par.tissue.t1periacquaductalgreyouter = T1PeriacquaductalGreyOuter ;
MRX.Par.tissue.t1periacquaductalgrey = T1PeriacquaductalGrey ;
MRX.Par.tissue.t1substantianigra = T1SubstantiaNigra ;
MRX.Par.tissue.t1pons = T1Pons ;
MRX.Par.tissue.t1amygdala = T1Amygdala ;
MRX.Par.tissue.t1fornix = T1Fornix ;
MRX.Par.tissue.t1putamen = T1Putamen ;
MRX.Par.tissue.t1thalamus = T1Thalamus ;
MRX.Par.tissue.t1globuspallidus = T1GlobusPallidus ;
MRX.Par.tissue.t1caudate = T1Caudate ;
MRX.Par.tissue.t1thirdventricle = T1ThirdVentricle ;
MRX.Par.tissue.t1lateralventricle = T1LateralVentricle ;
MRX.Par.tissue.t1fourthventricle = T1FourthVentricle ;
MRX.Par.tissue.t1brain = T1Brain ;
MRX.Par.tissue.t1mamillarybodies = T1MamillaryBodies ;
MRX.Par.tissue.t1pinealgland = T1PinealGland ;
MRX.Par.tissue.t1cerebralaqueduct = T1CerebralAqueduct ;
MRX.Par.tissue.t1septumpellicudum = T1SeptumPellicudum ;

MRX.Par.tissue.t2eye = T2Eye ;
MRX.Par.tissue.t2lens = T2Lens ;
MRX.Par.tissue.t2sinus = T2Sinus ;
MRX.Par.tissue.t2corticalbone = T2CorticalBone ;
MRX.Par.tissue.t2bonemarrow = T2BoneMarrow ;
MRX.Par.tissue.t2skin = T2Skin ;
MRX.Par.tissue.t2spine = T2Spine ;
MRX.Par.tissue.t2spinalcord = T2SpinalCord ;
MRX.Par.tissue.t2salivary = T2Salivary ;
MRX.Par.tissue.t2body = T2Body ;
MRX.Par.tissue.t2cartilage = T2Cartilage ;
MRX.Par.tissue.t2muscle = T2Muscle ;
MRX.Par.tissue.t2greymatter = T2GreyMatter ;
MRX.Par.tissue.t2whitematter = T2WhiteMatter ;
MRX.Par.tissue.t2internalcapsule = T2InternalCapsule ;
MRX.Par.tissue.t2anteriorcommissure = T2AnteriorCommissure ;
MRX.Par.tissue.t2cerebellum = T2cerebellum ;
MRX.Par.tissue.t2cerebralpeduncles = T2CerebralPeduncles ;
MRX.Par.tissue.t2superiorcerebellarpeduncle = T2SuperiorCerebellarPeduncle ;
MRX.Par.tissue.t2middlecerebellarpeduncle = T2MiddleCerebellarPeduncle ;
MRX.Par.tissue.t2corpuscallosum = T2CorpusCallosum ;
MRX.Par.tissue.t2hippocampus = T2Hippocampus ;
MRX.Par.tissue.t2midbrain = T2Midbrain ;
MRX.Par.tissue.t2tegmentumofmidbrain = T2TegmentumOfMidbrain ;
MRX.Par.tissue.t2medulla = T2Medulla ;
MRX.Par.tissue.t2medullarypyramids = T2MedullaryPyramids ;
MRX.Par.tissue.t2inferiorolive = T2InferiorOlive ;
MRX.Par.tissue.t2superiorcolliculus = T2SuperiorColliculus ;
MRX.Par.tissue.t2inferiorcolliculus = T2InferiorColliculus ;
MRX.Par.tissue.t2periacquaductalgreyouter = T2PeriacquaductalGreyOuter ;
MRX.Par.tissue.t2periacquaductalgrey = T2PeriacquaductalGrey ;
MRX.Par.tissue.t2substantianigra = T2SubstantiaNigra ;
MRX.Par.tissue.t2pons = T2Pons ;
MRX.Par.tissue.t2amygdala = T2Amygdala ;
MRX.Par.tissue.t2fornix = T2Fornix ;
MRX.Par.tissue.t2putamen = T2Putamen ;
MRX.Par.tissue.t2thalamus = T2Thalamus ;
MRX.Par.tissue.t2globuspallidus = T2GlobusPallidus ;
MRX.Par.tissue.t2caudate = T2Caudate ;
MRX.Par.tissue.t2thirdventricle = T2ThirdVentricle ;
MRX.Par.tissue.t2lateralventricle = T2LateralVentricle ;
MRX.Par.tissue.t2fourthventricle = T2FourthVentricle ;
MRX.Par.tissue.t2brain = T2Brain ;
MRX.Par.tissue.t2mamillarybodies = T2MamillaryBodies ;
MRX.Par.tissue.t2pinealgland = T2PinealGland ;
MRX.Par.tissue.t2cerebralaqueduct = T2CerebralAqueduct ;
MRX.Par.tissue.t2septumpellicudum = T2SeptumPellicudum ;

MRX.Par.scan.seq           = Sequence;
MRX.Par.scan.name          = Name;
MRX.Par.scan.bbox           = BoundingBox;
MRX.Par.scan.lowpass        = LowPassFilt;
MRX.Par.scan.snr            = SNR;
MRX.Par.scan.coils          = Coils;
MRX.Par.scan.coildist       = CoilDist;
MRX.Par.scan.coilsperrow    = CoilsPerRow;
MRX.Par.scan.rotation       = pi*RotationXYZ/180;
MRX.Par.scan.resizeimage    = ResizeImage;
if Frames>0 % only overwrite Par.scan.frames, if Frames ~= 0
    MRX.Par.scan.frames     = Frames;
    xcat_segments =1;
    MRX.Par.scan.frames_xcat=1;
    frames_max              = MRX.Par.scan.frames_xcat/xcat_segments;
    MRX.Par.scan.phases     = round( linspace(1,frames_max,Frames) );
    
else
    MRX.Par.scan.frames     = MRX.Par.scan.frames_xcat/MRX.Par.scan.segments;
    MRX.Par.scan.phases     = 1:MRX.Par.scan.frames;
end
MRX.Par.scan.trajectory     = Trajectory;
MRX.Par.scan.undersample    = Undersample;

% --------------------------------------------------------------------
%   Error checks
% --------------------------------------------------------------------
if mod(max(MRX.Par.scan.frames),1)>0 % check if #frames is whole number
    error('Number of frames must be an integer value. Check number of segments in Brainpar.m and number of XCAT .bin files!')
end
if exist('frames_max','var') && frames_max < MRX.Par.scan.frames
    error('Number of XCAT phases < desired number of phases. Set Frames <= %d in Brainpar.m',frames_max)
end

end