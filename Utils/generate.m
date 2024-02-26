function pars = generate(filepath,filename,outname,c,sequence)


%----------------------------
% path and filename for XCAT subject
%----------------------------
% filepath = filepath;
% filename = filename;   % 'VSM_bodysc_1_act_1.bin';
subject = outname;             %'_VSM_bodysc_1_ED';
Sequence=sequence;
% c =  number of contrast you would like to generate on one virutal subject
%----------------------------
% tissue parameters reference is the deliverable report
%----------------------------
ntis = c ; % number  of samples for tissue property value
T1.eye =  0 + 0*randn(ntis,1);
T2.eye =  0 + 0*randn(ntis,1);
PD.eye =  0 + 0*randn(ntis,1);

T1.lens =  T1.eye;
T2.lens =  T2.eye;
PD.lens =  PD.eye;

T1.sinus =  T1.eye;
T2.sinus =  T2.eye;
PD.sinus =  PD.eye;

T1.corticalbone =  269 + 2*randn(ntis,1);
T2.corticalbone =  55 + 2*randn(ntis,1);
PD.corticalbone =  0.08 + 0*randn(ntis,1);

% T1.bonemarrow =  350 + 73*randn(ntis,1);
% T2.bonemarrow =  237 + 50*randn(ntis,1);
% PD.bonemarrow =  1.0 + 0*randn(ntis,1);

T1.bonemarrow =   T1.eye;
T2.bonemarrow =   T2.eye;
PD.bonemarrow =   PD.eye;

T1.skin =  T1.eye;
T2.skin =  T2.eye;
PD.skin =  PD.eye;

T1.spine =  515 + 2*randn(ntis,1);
T2.spine =  101 + 2*randn(ntis,1);
PD.spine =  0.4 + 0*randn(ntis,1);

T1.spinalcord =  780 + 8*randn(ntis,1);
T2.spinalcord =  101 + 2*randn(ntis,1);
PD.spinalcord =  0.5 + 0*randn(ntis,1);

T1.salivary =  T1.eye;
T2.salivary =  T2.eye;
PD.salivary =  PD.eye;

T1.body =  T1.eye;
T2.body =  T2.eye;
PD.body =  PD.eye;

T1.cartilage =  1168 + 18*randn(ntis,1);
T2.cartilage =  27 + 3*randn(ntis,1);
PD.cartilage =  0.08 + 0*randn(ntis,1);

T1.muscle =  1008  + 20*randn(ntis,1);
T2.muscle =  44 + 6*randn(ntis,1);
PD.muscle =  1 + 0*randn(ntis,1);

T1.greymatter =  1433 + 18*randn(ntis,1);
T2.greymatter =  92 + 2*randn(ntis,1);
PD.greymatter =  0.86 + 0*randn(ntis,1);

T1.whitematter =  999 + 27*randn(ntis,1);
% T1.whitematter =  894 + 27*randn(ntis,1); % elongated WM T1 for MS patient
T2.whitematter =  75 + 1*randn(ntis,1);
PD.whitematter =  0.77 + 0*randn(ntis,1);


T1.cerebellumWM = T1.whitematter;
T2.cerebellumWM = T2.whitematter;
PD.cerebellumWM = PD.whitematter;

T1.cerebellumGM = T1.greymatter;
T2.cerebellumGM = T2.greymatter;
PD.cerebellumGM = PD.greymatter;

T1.internalcapsule =  T1.whitematter;
T2.internalcapsule =  T2.whitematter;
PD.internalcapsule =  PD.whitematter;

T1.anteriorcommisure =  T1.whitematter;
T2.anteriorcommisure =  T2.whitematter;
PD.anteriorcommisure =  PD.whitematter;

T1.cerebellum =  T1.whitematter;
T2.cerebellum =  T2.whitematter;
PD.cerebellum =  PD.whitematter;

T1.cerebralpeduncles =  T1.whitematter;
T2.cerebralpeduncles =  T2.whitematter;
PD.cerebralpeduncles =  PD.whitematter;

T1.superiorcerebellarpeduncle =  T1.whitematter;
T2.superiorcerebellarpeduncle =  T2.whitematter;
PD.superiorcerebellarpeduncle =  PD.whitematter;

T1.middlecerebellarpeduncle =  T1.whitematter;
T2.middlecerebellarpeduncle =  T2.whitematter;
PD.middlecerebellarpeduncle =  PD.whitematter;

T1.corpuscallosum =  T1.whitematter;
T2.corpuscallosum =  T2.whitematter;
PD.corpuscallosum =  PD.whitematter;

T1.hippocampus =  T1.greymatter;
T2.hippocampus =  T2.greymatter;
PD.hippocampus =  PD.greymatter;

T1.midbrain =  T1.whitematter;
T2.midbrain =  T2.whitematter;
PD.midbrain =  PD.whitematter;

T1.tegmentumofmidbrain =  T1.whitematter;
T2.tegmentumofmidbrain =  T2.whitematter;
PD.tegmentumofmidbrain =  PD.whitematter;

T1.medulla =  T1.whitematter;
T2.medulla =  T2.whitematter;
PD.medulla =  PD.whitematter;

T1.medullarypyramids =  T1.whitematter;
T2.medullarypyramids =  T2.whitematter;
PD.medullarypyramids =  PD.whitematter;

T1.inferiorolive =  T1.whitematter;
T2.inferiorolive =  T2.whitematter;
PD.inferiorolive =  PD.whitematter;

T1.superiorcolliculus =  T1.whitematter;
T2.superiorcolliculus =  T2.whitematter;
PD.superiorcolliculus =  PD.whitematter;

T1.inferiorcolliculus =  T1.whitematter;
T2.inferiorcolliculus =  T2.whitematter;
PD.inferiorcolliculus =  PD.whitematter;

T1.periacquaductalgreyouter =  T1.whitematter;
T2.periacquaductalgreyouter =  T2.whitematter;
PD.periacquaductalgreyouter =  PD.whitematter;

T1.periacquaductalgrey =  T1.whitematter;
T2.periacquaductalgrey =  T2.whitematter;
PD.periacquaductalgrey =  PD.whitematter;

T1.substantianigra =  T1.whitematter;
T2.substantianigra =  T2.whitematter;
PD.substantianigra =  PD.whitematter;

T1.pons =  T1.whitematter;
T2.pons =  T2.whitematter;
PD.pons =  PD.whitematter;

T1.amygdala =  T1.greymatter;
T2.amygdala =  T2.greymatter;
PD.amygdala =  PD.greymatter;

T1.fornix =  T1.whitematter;
T2.fornix =  T2.whitematter;
PD.fornix =  PD.whitematter;


T1.putamen = 1433 + 18*randn(ntis,1);
T2.putamen =  81 + 3*randn(ntis,1);
PD.putamen =  PD.greymatter;


T1.thalamus =  1433 + 18*randn(ntis,1);
T2.thalamus =  95 + 1.4*randn(ntis,1);
PD.thalamus =  PD.greymatter;


T1.globuspallidus = 1433 + 18*randn(ntis,1);
T2.globuspallidus =  66 + 3*randn(ntis,1);
PD.globuspallidus =  PD.greymatter;


T1.caudate = 1433 + 18*randn(ntis,1);
T2.caudate =  111 + 7*randn(ntis,1);
PD.caudate =  PD.greymatter;


T1.brain =  4200 + 263*randn(ntis,1);
T2.brain =  1500 + 7*randn(ntis,1);
PD.brain =  1 + 0*randn(ntis,1);

T1.thirdventricle =  T1.brain;
T2.thirdventricle =  T2.brain;
PD.thirdventricle =  PD.brain;

T1.lateralventricle =  T1.brain;
T2.lateralventricle =  T2.brain;
PD.lateralventricle =  PD.brain;

T1.fourthventricle =  T1.brain;
T2.fourthventricle =  T2.brain;
PD.fourthventricle =  PD.brain;


T1.mamillarybodies =  T1.brain;
T2.mamillarybodies =  T2.brain;
PD.mamillarybodies =  PD.brain;

T1.pinealgland =  T1.brain;
T2.pinealgland =  T2.brain;
PD.pinealgland =  PD.brain;

T1.cerebralaqueduct =  T1.brain;
T2.cerebralaqueduct =  T2.brain;
PD.cerebralaqueduct =  PD.brain;

T1.septumpellicudum =  T1.whitematter;
T2.septumpellicudum =  T2.whitematter;
PD.septumpellicudum =  PD.whitematter;


T1.lesn =  1353 + 25*randn(ntis,1);
T2.lesn =  237 + 1*randn(ntis,1);
PD.lesn =  0.76 + 0*randn(ntis,1);

%%
%----------------------------
% sequence parameters
%----------------------------
%----------------------------
% SNR level
%----------------------------
nsnr = c ; % number  of samples for SNR value
nseq = c ; % number  of samples for sequence paramters value

if (Sequence == "SPINECHO")
    % spin echo t2w
    TR = 3000 + 40*randn(nseq,1);
    TE = 100 + 1*randn(nseq,1);
    FA = 90 + 0*randn(nseq,1);
    TI = 0 + 0*randn(nseq,1);
    SNR = 400 + 50 * randn(nsnr,1);
end
if (Sequence == "Ernst")
    % fixed for SR
    TR = 18;
    TE = 10;
    FA = 30;
    
    TI = 0;
    SNR = 8000;
    
end
if (Sequence == "IR")
    % inverion recovery t1w flair
    TR = 2100 + 0*randn(nseq,1);
    TE = 10 + 0*randn(nseq,1);
    FA = 90 + 0*randn(nseq,1);
    TI = 900 + 0*randn(nseq,1);
    SNR = 1300 + 50 * randn(nsnr,1);
    
end


%  Res2 = [0.7 0.7 0.7];
Res2 = [1 1 1];
%   Res2 = [2 2 2];

pars.SNR = SNR;
pars.TR = TR;
pars.TE = TE;
pars.T1 = TI;
pars.FA = FA;
pars.Res = Res2;
pars.T1 = T1;
pars.T2 = T2;
pars.PD = PD;




for i=1:c
    
    MRXCAT_CMR_CINE(filename,filepath,'snr',SNR(i),'flip',FA(i),'tr',TR(i),'te',TE(i),'ti',TI(i) ...
        ,'resolution',Res2 ...
        ,'t1eye',T1.eye(i), 't2eye',T2.eye(i),'pdeye', PD.eye(i) ...
        ,'t1lens',T1.lens(i), 't2lens',T2.lens(i),'pdlens', PD.lens(i) ...
        ,'t1sinus',T1.sinus(i), 't2sinus',T2.sinus(i),'pdsinus', PD.sinus(i) ...
        ,'t1corticalbone',T1.corticalbone(i), 't2corticalbone',T2.corticalbone(i),'pdcorticalbone', PD.corticalbone(i) ...
        ,'t1bonemarrow',T1.bonemarrow(i), 't2bonemarrow',T2.bonemarrow(i),'pdbonemarrow', PD.bonemarrow(i) ...
        ,'t1skin',T1.skin(i), 't2skin',T2.skin(i),'pdskin', PD.skin(i) ...
        ,'t1spine',T1.spine(i), 't2spine',T2.spine(i),'pdspine', PD.spine(i) ...
        ,'t1spinalcord',T1.spinalcord(i), 't2spinalcord',T2.spinalcord(i),'pdspinalcord', PD.spinalcord(i) ...
        ,'t1salivary',T1.salivary(i), 't2salivary',T2.salivary(i),'pdsalivary', PD.salivary(i) ...
        ,'t1body',T1.body(i), 't2body',T2.body(i),'pdbody', PD.body(i) ...
        ,'t1cartilage',T1.cartilage(i), 't2cartilage',T2.cartilage(i),'pdcartilage', PD.cartilage(i) ...
        ,'t1muscle',T1.muscle(i), 't2muscle',T2.muscle(i),'pdmuscle', PD.muscle(i) ...
        ,'t1greymatter',T1.greymatter(i), 't2greymatter',T2.greymatter(i),'pdgreymatter', PD.greymatter(i) ...
        ,'t1whitematter',T1.whitematter(i), 't2whitematter',T2.whitematter(i),'pdwhitematter', PD.whitematter(i) ...
        ,'t1internalcapsule',T1.internalcapsule(i), 't2internalcapsule',T2.internalcapsule(i),'pdinternalcapsule', PD.internalcapsule(i) ...
        ,'t1anteriorcommisure',T1.anteriorcommisure(i), 't2anteriorcommisure',T2.anteriorcommisure(i),'pdanteriorcommisure', PD.anteriorcommisure(i) ...
        ,'t1cerebellum',T1.cerebellum(i), 't2cerebellum',T2.cerebellum(i),'pdcerebellum', PD.cerebellum(i) ...
        ,'t1cerebralpeduncles',T1.cerebralpeduncles(i), 't2cerebralpeduncles',T2.cerebralpeduncles(i),'pdcerebralpeduncles', PD.cerebralpeduncles(i) ...
        ,'t1superiorcerebellarpeduncle',T1.superiorcerebellarpeduncle(i), 't2superiorcerebellarpeduncle',T2.superiorcerebellarpeduncle(i),'pdsuperiorcerebellarpeduncle', PD.superiorcerebellarpeduncle(i) ...
        ,'t1middlecerebellarpeduncle',T1.middlecerebellarpeduncle(i), 't2middlecerebellarpeduncle',T2.middlecerebellarpeduncle(i),'pdmiddlecerebellarpeduncle', PD.middlecerebellarpeduncle(i) ...
        ,'t1corpuscallosum',T1.corpuscallosum(i), 't2corpuscallosum',T2.corpuscallosum(i),'pdcorpuscallosum', PD.corpuscallosum(i) ...
        ,'t1hippocampus',T1.hippocampus(i), 't2hippocampus',T2.hippocampus(i),'pdhippocampus', PD.hippocampus(i) ...
        ,'t1midbrain',T1.midbrain(i), 't2midbrain',T2.midbrain(i),'pdmidbrain', PD.midbrain(i) ...
        ,'t1tegmentumofmidbrain',T1.tegmentumofmidbrain(i), 't2tegmentumofmidbrain',T2.tegmentumofmidbrain(i),'pdtegmentumofmidbrain', PD.tegmentumofmidbrain(i) ...
        ,'t1medulla',T1.medulla(i), 't2medulla',T2.medulla(i),'pdmedulla', PD.medulla(i) ...
        ,'t1medullarypyramids',T1.medullarypyramids(i), 't2medullarypyramids',T2.medullarypyramids(i),'pdmedullarypyramids', PD.medullarypyramids(i) ...
        ,'t1inferiorolive',T1.inferiorolive(i), 't2inferiorolive',T2.inferiorolive(i),'pdinferiorolive', PD.inferiorolive(i) ...
        ,'t1superiorcolliculus',T1.superiorcolliculus(i), 't2superiorcolliculus',T2.superiorcolliculus(i),'pdsuperiorcolliculus', PD.superiorcolliculus(i) ...
        ,'t1inferiorcolliculus',T1.inferiorcolliculus(i), 't2inferiorcolliculus',T2.inferiorcolliculus(i),'pdinferiorcolliculus', PD.inferiorcolliculus(i) ...
        ,'t1periacquaductalgreyouter',T1.periacquaductalgreyouter(i), 't2periacquaductalgreyouter',T2.periacquaductalgreyouter(i),'pdperiacquaductalgreyouter', PD.periacquaductalgreyouter(i) ...
        ,'t1periacquaductalgrey',T1.periacquaductalgrey(i), 't2periacquaductalgrey',T2.periacquaductalgrey(i),'pdperiacquaductalgrey', PD.periacquaductalgrey(i) ...
        ,'t1substantianigra',T1.substantianigra(i), 't2substantianigra',T2.substantianigra(i),'pdsubstantianigra', PD.substantianigra(i) ...
        ,'t1pons',T1.pons(i), 't2pons',T2.pons(i),'pdpons', PD.pons(i) ...
        ,'t1amygdala',T1.amygdala(i), 't2amygdala',T2.amygdala(i),'pdamygdala', PD.amygdala(i) ...
        ,'t1fornix',T1.fornix(i), 't2fornix',T2.fornix(i),'pdfornix', PD.fornix(i) ...
        ,'t1putamen',T1.putamen(i), 't2putamen',T2.putamen(i),'pdputamen', PD.putamen(i) ...
        ,'t1thalamus',T1.thalamus(i), 't2thalamus',T2.thalamus(i),'pdthalamus', PD.thalamus(i) ...
        ,'t1globuspallidus',T1.globuspallidus(i), 't2globuspallidus',T2.globuspallidus(i),'pdglobuspallidus', PD.globuspallidus(i) ...
        ,'t1caudate',T1.caudate(i), 't2caudate',T2.caudate(i),'pdcaudate', PD.caudate(i) ...
        ,'t1thirdventricle',T1.thirdventricle(i), 't2thirdventricle',T2.thirdventricle(i),'pdthirdventricle', PD.thirdventricle(i) ...
        ,'t1lateralventricle',T1.lateralventricle(i), 't2lateralventricle',T2.lateralventricle(i),'pdlateralventricle', PD.lateralventricle(i) ...
        ,'t1fourthventricle',T1.fourthventricle(i), 't2fourthventricle',T2.fourthventricle(i),'pdfourthventricle', PD.fourthventricle(i) ...
        ,'t1brain',T1.brain(i), 't2brain',T2.brain(i),'pdbrain', PD.brain(i) ...
        ,'t1mamillarybodies',T1.mamillarybodies(i), 't2mamillarybodies',T2.mamillarybodies(i),'pdmamillarybodies', PD.mamillarybodies(i) ...
        ,'t1pinealgland',T1.pinealgland(i), 't2pinealgland',T2.pinealgland(i),'pdpinealgland', PD.pinealgland(i) ...
        ,'t1cerebralaqueduct',T1.cerebralaqueduct(i), 't2cerebralaqueduct',T2.cerebralaqueduct(i),'pdcerebralaqueduct', PD.cerebralaqueduct(i) ...
        ,'t1septumpellicudum',T1.septumpellicudum(i), 't2septumpellicudum',T2.septumpellicudum(i),'pdseptumpellicudum', PD.septumpellicudum(i) ...
        ,'t1lesn',T1.lesn(i), 't2lesn',T2.lesn(i),'pdlesn', PD.lesn(i) ...
        ,'t1cerebellumWM',T1.cerebellumWM(i), 't2cerebellumWM',T2.cerebellumWM(i),'pdcerebellumWM', PD.cerebellumWM(i) ...
        ,'t1cerebellumGM',T1.cerebellumGM(i), 't2cerebellumGM',T2.cerebellumGM(i),'pdcerebellumGM', PD.cerebellumGM(i) ...
        ,'simname',['_' subject '_contrast_' num2str(i)])
    
    
end
end
