%=========================================================================
%
% MRXCAT_Brain MRXCAT class for structural brain imaging
%
% The MRXCAT_Brain class contains specific MRXCAT methods.
% Methods common with other classes are in the MRXCAT superclass.
%
% PROPERTIES:
%
%           All properties are listed in MRXCAT.m
%
% METHODS:  (MRX refers to the MRXCAT instance)
%
%           Brainpar( MRX, filename )    External parameter function (Brainpar.m)
%                                       Use Brainpar.m to modify MRXCAT parameters
%           computeNoiseStdDev( MRX, sen )
%                                       Compute standard deviation of object noise based on desired SNR
%                                       in MRX.Par.scan.snr. Signals are averaged over the brain white matter only
%                                       to get "brain SNR".
%           mapTissueProps( MRX, data ) Assign tissue properties (T1, T2, rho) and apply signal model.
%                                       For brain, the signal model is a earnest angle solution.
%                                       Create tissue masks with different XCAT mask values.
%           extractSegment( MRX, segm, img )
%                                       Extract k-space segment segm from image img (in image space).
%                                       Total no. of segments: MRX.Par.scan.segments
%           radialResample( MRX, img )  Calculate and apply radial trajectory. The Cartesian input
%                                       is resampled to a radial trajectory (MRX.Par.scan.trajectory)
%                                       including optional undersampling (factor MRX.Par.scan.undersample).
%
% STATIC METHODS:
%
%           getRadWeights2D( no_samples, no_profiles, alt_prof, gafl, normfl )
%                                       Calculate sampling weights for radial trajectory. For more details,
%                                       cf. function help.
%           buildRadTraj2D( no_samples, no_profiles, alt_prof, gafl, normfl, dim_t, t_offset, dim_z, z_offset)
%                                       Create radial trajectory. For more details, cf. function help.
%
%
% WEBSITE: 	http://www.biomed.ee.ethz.ch/mrxcat
%
%=========================================================================

%=========================================================================
% VERSION HISTORY:
%		130129SK INITIAL VERSION
%		130326LW DXCAT2 W/ ANGULATION - v0.1
%       130625LW ADAPTATION TO MRXCAT CMR PERF - v0.8
%       130828LW RADIAL TRAJECTORIES - v0.9
%       140130LW OO IMPLEMENTATION - v1.0
%       140202LW SEGMENTED ACQ, LOW-PASS FILTER - v1.1
%
% AUTHORS:  Lukas Wissmann, Sebastian Kozerke, Claudio Santelli
%           Institute for Biomedical Engineering, University and ETH Zurich
%
%=========================================================================

classdef MRXCAT_Brain < MRXCAT
    
    properties
        
    end
    
    methods
        function MRX = MRXCAT_Brain( filename, filepath, varargin )
            
            MRX = MRX@MRXCAT();
            
            % ------------------------------------------------------------
            % Load and assign parameters (file
            % ------------------------------------------------------------
            if ~exist('filename','var'), filename = ''; end
            MRX.Brainpar( filename, filepath );
            
            % ------------------------------------------------------------
            % Check for additional inputs
            % ------------------------------------------------------------
            for i = 1:length( varargin )
                if any( strcmpi( varargin{i}, {'snr'} ))
                    MRX.Par.scan.snr = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'tr'} ))
                    MRX.Par.scan.tr = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'te'} ))
                    MRX.Par.scan.te = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'ti'} ))
                    MRX.Par.scan.ti = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'frames'} ))
                    MRX.Par.scan.frames = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'coils'} ))
                    MRX.Par.scan.coils = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'resolution'} ))
                    MRX.Par.scan.lowpass_str = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'flip'} ))
                    MRX.Par.scan.flip = pi*varargin{i+1}/180;
                end
                if any( strcmpi( varargin{i}, {'t1eye'} ))
                    MRX.Par.tissue.t1eye  = varargin{i+1};
                    MRX.Par.tissue.t2eye  = varargin{i+3};
                    MRX.Par.tissue.rhoeye  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1lens'} ))
                    MRX.Par.tissue.t1lens  = varargin{i+1};
                    MRX.Par.tissue.t2lens  = varargin{i+3};
                    MRX.Par.tissue.rholens  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1sinus'} ))
                    MRX.Par.tissue.t1sinus  = varargin{i+1};
                    MRX.Par.tissue.t2sinus  = varargin{i+3};
                    MRX.Par.tissue.rhosinus  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1corticalbone'} ))
                    MRX.Par.tissue.t1corticalbone  = varargin{i+1};
                    MRX.Par.tissue.t2corticalbone  = varargin{i+3};
                    MRX.Par.tissue.rhocorticalbone  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1bonemarrow'} ))
                    MRX.Par.tissue.t1bonemarrow  = varargin{i+1};
                    MRX.Par.tissue.t2bonemarrow  = varargin{i+3};
                    MRX.Par.tissue.rhobonemarrow  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1skin'} ))
                    MRX.Par.tissue.t1skin  = varargin{i+1};
                    MRX.Par.tissue.t2skin  = varargin{i+3};
                    MRX.Par.tissue.rhoskin  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1spine'} ))
                    MRX.Par.tissue.t1spine  = varargin{i+1};
                    MRX.Par.tissue.t2spine  = varargin{i+3};
                    MRX.Par.tissue.rhospine  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1spinalcord'} ))
                    MRX.Par.tissue.t1spinalcord  = varargin{i+1};
                    MRX.Par.tissue.t2spinalcord  = varargin{i+3};
                    MRX.Par.tissue.rhospinalcord  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1salivary'} ))
                    MRX.Par.tissue.t1salivary  = varargin{i+1};
                    MRX.Par.tissue.t2salivary  = varargin{i+3};
                    MRX.Par.tissue.rhosalivary  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1body'} ))
                    MRX.Par.tissue.t1body  = varargin{i+1};
                    MRX.Par.tissue.t2body  = varargin{i+3};
                    MRX.Par.tissue.rhobody  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1cartilage'} ))
                    MRX.Par.tissue.t1cartilage  = varargin{i+1};
                    MRX.Par.tissue.t2cartilage  = varargin{i+3};
                    MRX.Par.tissue.rhocartilage  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1muscle'} ))
                    MRX.Par.tissue.t1muscle  = varargin{i+1};
                    MRX.Par.tissue.t2muscle  = varargin{i+3};
                    MRX.Par.tissue.rhomuscle  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1greymatter'} ))
                    MRX.Par.tissue.t1greymatter  = varargin{i+1};
                    MRX.Par.tissue.t2greymatter  = varargin{i+3};
                    MRX.Par.tissue.rhogreymatter  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1whitematter'} ))
                    MRX.Par.tissue.t1whitematter  = varargin{i+1};
                    MRX.Par.tissue.t2whitematter  = varargin{i+3};
                    MRX.Par.tissue.rhowhitematter  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1internalcapsule'} ))
                    MRX.Par.tissue.t1internalcapsule  = varargin{i+1};
                    MRX.Par.tissue.t2internalcapsule  = varargin{i+3};
                    MRX.Par.tissue.rhointernalcapsule  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1anteriorcommisure'} ))
                    MRX.Par.tissue.t1anteriorcommisure  = varargin{i+1};
                    MRX.Par.tissue.t2anteriorcommisure  = varargin{i+3};
                    MRX.Par.tissue.rhoanteriorcommisure  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1cerebellum'} ))
                    MRX.Par.tissue.t1cerebellum  = varargin{i+1};
                    MRX.Par.tissue.t2cerebellum  = varargin{i+3};
                    MRX.Par.tissue.rhocerebellum  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1cerebralpeduncles'} ))
                    MRX.Par.tissue.t1cerebralpeduncles  = varargin{i+1};
                    MRX.Par.tissue.t2cerebralpeduncles  = varargin{i+3};
                    MRX.Par.tissue.rhocerebralpeduncles  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1superiorcerebellarpeduncle'} ))
                    MRX.Par.tissue.t1superiorcerebellarpeduncle  = varargin{i+1};
                    MRX.Par.tissue.t2superiorcerebellarpeduncle  = varargin{i+3};
                    MRX.Par.tissue.rhosuperiorcerebellarpeduncle  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1middlecerebellarpeduncle'} ))
                    MRX.Par.tissue.t1middlecerebellarpeduncle  = varargin{i+1};
                    MRX.Par.tissue.t2middlecerebellarpeduncle  = varargin{i+3};
                    MRX.Par.tissue.rhomiddlecerebellarpeduncle  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1corpuscallosum'} ))
                    MRX.Par.tissue.t1corpuscallosum  = varargin{i+1};
                    MRX.Par.tissue.t2corpuscallosum  = varargin{i+3};
                    MRX.Par.tissue.rhocorpuscallosum  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1hippocampus'} ))
                    MRX.Par.tissue.t1hippocampus  = varargin{i+1};
                    MRX.Par.tissue.t2hippocampus  = varargin{i+3};
                    MRX.Par.tissue.rhohippocampus  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1midbrain'} ))
                    MRX.Par.tissue.t1midbrain  = varargin{i+1};
                    MRX.Par.tissue.t2midbrain  = varargin{i+3};
                    MRX.Par.tissue.rhomidbrain  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1tegmentumofmidbrain'} ))
                    MRX.Par.tissue.t1tegmentumofmidbrain  = varargin{i+1};
                    MRX.Par.tissue.t2tegmentumofmidbrain  = varargin{i+3};
                    MRX.Par.tissue.rhotegmentumofmidbrain  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1medulla'} ))
                    MRX.Par.tissue.t1medulla  = varargin{i+1};
                    MRX.Par.tissue.t2medulla  = varargin{i+3};
                    MRX.Par.tissue.rhomedulla  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1medullarypyramids'} ))
                    MRX.Par.tissue.t1medullarypyramids  = varargin{i+1};
                    MRX.Par.tissue.t2medullarypyramids  = varargin{i+3};
                    MRX.Par.tissue.rhomedullarypyramids  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1inferiorolive'} ))
                    MRX.Par.tissue.t1inferiorolive  = varargin{i+1};
                    MRX.Par.tissue.t2inferiorolive  = varargin{i+3};
                    MRX.Par.tissue.rhoinferiorolive  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1superiorcolliculus'} ))
                    MRX.Par.tissue.t1superiorcolliculus  = varargin{i+1};
                    MRX.Par.tissue.t2superiorcolliculus  = varargin{i+3};
                    MRX.Par.tissue.rhosuperiorcolliculus  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1inferiorcolliculus'} ))
                    MRX.Par.tissue.t1inferiorcolliculus  = varargin{i+1};
                    MRX.Par.tissue.t2inferiorcolliculus  = varargin{i+3};
                    MRX.Par.tissue.rhoinferiorcolliculus  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1periacquaductalgreyouter'} ))
                    MRX.Par.tissue.t1periacquaductalgreyouter  = varargin{i+1};
                    MRX.Par.tissue.t2periacquaductalgreyouter  = varargin{i+3};
                    MRX.Par.tissue.rhoperiacquaductalgreyouter  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1periacquaductalgrey'} ))
                    MRX.Par.tissue.t1periacquaductalgrey  = varargin{i+1};
                    MRX.Par.tissue.t2periacquaductalgrey  = varargin{i+3};
                    MRX.Par.tissue.rhoperiacquaductalgrey  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1substantianigra'} ))
                    MRX.Par.tissue.t1substantianigra  = varargin{i+1};
                    MRX.Par.tissue.t2substantianigra  = varargin{i+3};
                    MRX.Par.tissue.rhosubstantianigra  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1pons'} ))
                    MRX.Par.tissue.t1pons  = varargin{i+1};
                    MRX.Par.tissue.t2pons  = varargin{i+3};
                    MRX.Par.tissue.rhopons  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1amygdala'} ))
                    MRX.Par.tissue.t1amygdala  = varargin{i+1};
                    MRX.Par.tissue.t2amygdala  = varargin{i+3};
                    MRX.Par.tissue.rhoamygdala  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1fornix'} ))
                    MRX.Par.tissue.t1fornix  = varargin{i+1};
                    MRX.Par.tissue.t2fornix  = varargin{i+3};
                    MRX.Par.tissue.rhofornix  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1putamen'} ))
                    MRX.Par.tissue.t1putamen  = varargin{i+1};
                    MRX.Par.tissue.t2putamen  = varargin{i+3};
                    MRX.Par.tissue.rhoputamen  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1thalamus'} ))
                    MRX.Par.tissue.t1thalamus  = varargin{i+1};
                    MRX.Par.tissue.t2thalamus  = varargin{i+3};
                    MRX.Par.tissue.rhothalamus  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1globuspallidus'} ))
                    MRX.Par.tissue.t1globuspallidus  = varargin{i+1};
                    MRX.Par.tissue.t2globuspallidus  = varargin{i+3};
                    MRX.Par.tissue.rhoglobuspallidus  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1caudate'} ))
                    MRX.Par.tissue.t1caudate  = varargin{i+1};
                    MRX.Par.tissue.t2caudate  = varargin{i+3};
                    MRX.Par.tissue.rhocaudate  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1thirdventricle'} ))
                    MRX.Par.tissue.t1thirdventricle  = varargin{i+1};
                    MRX.Par.tissue.t2thirdventricle  = varargin{i+3};
                    MRX.Par.tissue.rhothirdventricle  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1lateralventricle'} ))
                    MRX.Par.tissue.t1lateralventricle  = varargin{i+1};
                    MRX.Par.tissue.t2lateralventricle  = varargin{i+3};
                    MRX.Par.tissue.rholateralventricle  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1fourthventricle'} ))
                    MRX.Par.tissue.t1fourthventricle  = varargin{i+1};
                    MRX.Par.tissue.t2fourthventricle  = varargin{i+3};
                    MRX.Par.tissue.rhofourthventricle  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1brain'} ))
                    MRX.Par.tissue.t1brain  = varargin{i+1};
                    MRX.Par.tissue.t2brain  = varargin{i+3};
                    MRX.Par.tissue.rhobrain  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1mamillarybodies'} ))
                    MRX.Par.tissue.t1mamillarybodies  = varargin{i+1};
                    MRX.Par.tissue.t2mamillarybodies  = varargin{i+3};
                    MRX.Par.tissue.rhomamillarybodies  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1pinealgland'} ))
                    MRX.Par.tissue.t1pinealgland  = varargin{i+1};
                    MRX.Par.tissue.t2pinealgland  = varargin{i+3};
                    MRX.Par.tissue.rhopinealgland  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1cerebralaqueduct'} ))
                    MRX.Par.tissue.t1cerebralaqueduct  = varargin{i+1};
                    MRX.Par.tissue.t2cerebralaqueduct  = varargin{i+3};
                    MRX.Par.tissue.rhocerebralaqueduct  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1septumpellicudum'} ))
                    MRX.Par.tissue.t1septumpellicudum  = varargin{i+1};
                    MRX.Par.tissue.t2septumpellicudum  = varargin{i+3};
                    MRX.Par.tissue.rhoseptumpellicudum  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1lesn'} ))
                    MRX.Par.tissue.t1lesn  = varargin{i+1};
                    MRX.Par.tissue.t2lesn  = varargin{i+3};
                    MRX.Par.tissue.rholesn  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1cerebellumWM'} ))
                    MRX.Par.tissue.t1lesn  = varargin{i+1};
                    MRX.Par.tissue.t2lesn  = varargin{i+3};
                    MRX.Par.tissue.rholesn  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'t1cerebellumGM'} ))
                    MRX.Par.tissue.t1lesn  = varargin{i+1};
                    MRX.Par.tissue.t2lesn  = varargin{i+3};
                    MRX.Par.tissue.rholesn  = varargin{i+5};
                end
                if any( strcmpi( varargin{i}, {'simname'} ))
                    MRX.Par.scan.name  = [MRX.Par.scan.name varargin{i+1}];
                    
                end
                
                if any( strcmpi( varargin{i}, {'demo_gui'} )) % overwrite PARs from demo GUI
                    par = varargin{i+1};
                    names1 = fieldnames(par.contrast);
                    for n=1:length(names1)
                        if isfield( par.contrast, names1{n} )
                            MRX.Par.contrast.(names1{n}) = par.contrast.(names1{n});
                        end
                    end
                    names2 = fieldnames(par.scan);
                    for n=1:length(names2)
                        if isfield( par.scan, names2{n} )
                            MRX.Par.scan.(names2{n}) = par.scan.(names2{n});
                        end
                    end
                    % currently no tissue parameters modifiable by GUI
                end
            end
            
            % --------------------------------------------------------------------
            %   Calculate coil sensitivies (Biot-Savart)
            % --------------------------------------------------------------------
            sen = MRX.calculateCoilMaps; sen = single(sen);
            
            %             % --------------------------------------------------------------------
            %             % Compute standard deviation factor for noise addition (SNR way):
            %             % stddev = mean ROI signal / desired SNR
            %             % --------------------------------------------------------------------
            %             MRX.Par.scan.noisestd = MRX.computeNoiseStdDev( sen );
            %
            st = 0; % initialize waitbar
            
            % --------------------------------------------------------------------
            %   Produce MRXCAT phantom loops over heart phases & k-space segments
            % --------------------------------------------------------------------
            for t=MRX.Par.scan.phases % heart phases
                for s=1:MRX.Par.scan.segments % k-space segments
                    
                    % --------------------------------------------------------------------
                    %   Read data
                    % --------------------------------------------------------------------
                    xcat_no = t+((s-1)*MRX.Par.scan.frames);
                    data = MRX.readImgData(xcat_no); data = single(data);
                    
                    % ----------------------------------------------------------------
                    %   Map MR tissue properties
                    % ----------------------------------------------------------------
                    [img, msk] = MRX.mapTissueProps(data); img = single(img); msk = single(msk);
                    
                    % ----------------------------------------------------------------
                    %   Low Pass Filter (Blur)
                    % ----------------------------------------------------------------
                    [img, msk , KSPACE] = MRX.lowPassFilter(img, msk);
                    
                    save('kspace.mat', 'KSPACE')
                    
                    
                    % ------------------------------------------------------------------
                    % stddev = mean ROI signal / desired SNR -- ne signal
                    % should be scaled to the 12 bit range first
                    % --------------------------------------------------------------------
                    MRX.Par.scan.noisestd = MRX.computeNoiseStdDev(KSPACE, msk);
                    
                    % ----------------------------------------------------------------
                    %   Add coils
                    % ----------------------------------------------------------------
                    img       = MRX.multiplyCoilMaps(img,sen);
                    
                    % ----------------------------------------------------------------
                    %   Add noise
                    % ----------------------------------------------------------------
                    [KSPACEnoisy,noi] = MRX.addNoise(KSPACE);
                    
                    %-----------------------------------------------------------------
                    % Save data to Nifti file
                    %-----------------------------------------------------------------
                    
                    MRX.saveNifti(KSPACEnoisy, msk);
                    fprintf('\n---------------Saving is Done------------------------\n');
                    
                    % ----------------------------------------------------------------
                    %   Extract needed k-space segment
                    % ----------------------------------------------------------------
                    if s==1; MRX.Ksp = zeros(size(img)); end
                    MRX.extractSegment( s, img );
                    if s==MRX.Par.scan.segments %% FFT after k-space is filled
                        img = MRXCAT.k2i( MRX.Ksp, [1 2 3]);
                    end
                    
                end
                MRX.Ksp = [];
                
                % -----------------------------------------------------------------
                %   Regrid data for radial trajectory
                % -----------------------------------------------------------------
                if strcmpi(MRX.Par.scan.trajectory, 'Radial') || strcmpi(MRX.Par.scan.trajectory, 'GoldenAngle')
                    img = MRX.radialResample( img );
                end
                
                % -----------------------------------------------------------------
                %   Save data to .mat file
                % -----------------------------------------------------------------
                %                 MRX.saveImgData(img, msk, noi, sen, t);
                %                 st = MRX.textwaitbar( t, st, 'Writing MRXCAT output data'); % update waitbar
                
                
            end
            
            % --------------------------------------------------------------------
            %   Save Parameters for Recon
            % --------------------------------------------------------------------
            MRX.saveParsForRecon( img );
            
        end
        
        
        %=========================================================================
        % Calc noise std deviation based on desired SNR
        
        function stdev = computeNoiseStdDev( MRX, KSPACE, mskn )
            
            image = abs(ifft2c(KSPACE));
            %      scale to 12 bit range 0-4096
            inmin = min(image(find(image)));
            inmax = max(image(find(image)));
            image = rescale(image,0,4095,'InputMin',inmin,'InputMax',inmax);
            
            roi         = mskn==MRX.Par.act.white_matter_act;
            
            
            % find biggest connected mask = heart (if > 1 connected region)
            conn        = bwconncomp(roi);
            if conn.NumObjects > 1
                for k=1:conn.NumObjects
                    len(k)  = length(conn.PixelIdxList{k});
                end
                [~,maxk]    = max(len);
                idxs        = conn.PixelIdxList{maxk};
                roi         = zeros(size(roi));
                roi(idxs)   = 1;
            end
            roiall = roi.*image;
            smean = mean(roiall(roiall~=0));
            sosmean     = sqrt(sum(abs(smean).^2));
            stdev       = sosmean/MRX.Par.scan.snr;
            fprintf('adding noise w/ std dev : %f\n',stdev);
            clear image
        end
        
        
        %=========================================================================
        % Calc signal intensities using tissue and sequence pars
        function [img,msk] = mapTissueProps(MRX, data)
            
            act = cell2mat(struct2cell(MRX.Par.act));
            tis = fieldnames(MRX.Par.act);
            act = [act;1;18;19];
            tis = [tis; 'lesn_activity';'cerebellumWM_act';'cerebellumGM_act'];
            img = single(zeros(size(data)));
            msk = uint8(zeros(size(data)));
            
            tr=MRX.Par.scan.tr;
            ti=MRX.Par.scan.ti;
            t2p=MRX.Par.scan.t2p;
            a=MRX.Par.scan.flip;
            te=MRX.Par.scan.te;
            sequence=MRX.Par.scan.seq;
            ntis=1;
            
            for z1=1:320
                for x1=1:320
                    for y1=1:320
                        
                        z1
                        i=data(x1,y1,z1);
                        s=(tis(find(act==i)));
                        s1=s(1);
                        rho = 0;
                        % ----------------------------------------------------------------
                        %   Select tissue type
                        % ----------------------------------------------------------------
                        switch char(s1)
                            % ------------------------------------------------------------
                            %   Brain xcat
                            % ------------------------------------------------------------
                            
                            case {'eye_activity','skin_activity','lens_activity','sinus_activity','body_activity','salivary_activity'}
                                rho = MRX.Par.tissue.rhoeye;
                                t1  = MRX.Par.tissue.t1eye;
                                t2  = MRX.Par.tissue.t2eye;
                            case {'cortical_bone_activity'}
                                rho = MRX.Par.tissue.rhocorticalbone;
                                t1  = MRX.Par.tissue.t1corticalbone;
                                t2  = MRX.Par.tissue.t2corticalbone;
                            case {'bone_marrow_activity'}
                                rho = MRX.Par.tissue.rhobonemarrow;
                                t1  = MRX.Par.tissue.t1bonemarrow;
                                t2  = MRX.Par.tissue.t2bonemarrow;
                            case {'muscle_activity'}
                                rho = MRX.Par.tissue.rhomuscle;
                                t1  = MRX.Par.tissue.t1muscle;
                                t2  = MRX.Par.tissue.t2muscle;
                            case {'spinal_cord_activity'}
                                rho = MRX.Par.tissue.rhospinalcord;
                                t1  = MRX.Par.tissue.t1spinalcord;
                                t2  = MRX.Par.tissue.t2spinalcord;
                            case {'cartilage_activity'}
                                rho = MRX.Par.tissue.rhocartilage;
                                t1  = MRX.Par.tissue.t1cartilage;
                                t2  = MRX.Par.tissue.t2cartilage;
                            case {'spine_activity'}
                                rho = MRX.Par.tissue.rhospine;
                                t1  = MRX.Par.tissue.t1spine;
                                t2  = MRX.Par.tissue.t2spine;
                            case {'grey_matter_act','Hippocampus_act','Amygdala_act','cerebellumGM_act'}
                                t1 =  1433 + 18*randn(ntis,1);
                                t2 =  92 + 2*randn(ntis,1);
                                rho =  0.86 + 0*randn(ntis,1);
                            case {'white_matter_act','Corpus_Callosum_act','septum_pellicudum_act','Internal_capsule_act','Anterior_commissure_act','cerebellum_act','cerebellumWM_act','Cerebral_peduncles_act','Superior_cerebellar_peduncle_act','Middle_cerebellar_peduncle_act','Midbrain_act','Tegmentum_of_midbrain_act','Medulla_act' ,'Medullary_pyramids_act','Inferior_olive_act','Superior_colliculus_act','Inferior_colliculus_act','Periacquaductal_grey_outer_act','Periacquaductal_grey_act' ,'Substantia_nigra_act','Pons_act','Fornix_act'}
                                rho = 0.77 + 0*randn(ntis,1);
                                t1  = 999 + 27*randn(ntis,1);
                                t2  = 75 + 1*randn(ntis,1);
                            case {'Putamen_act'}
                                rho = 0.86 + 0*randn(ntis,1);
                                t1  = 1433 + 18*randn(ntis,1);
                                t2  = 81 + 3*randn(ntis,1);
                            case {'Thalamus_act'}
                                rho = 0.86 + 0*randn(ntis,1);
                                t1  = 1433 + 18*randn(ntis,1);
                                t2  = 95 + 1.4*randn(ntis,1);
                            case {'Globus_pallidus_act'}
                                rho = 0.86 + 0*randn(ntis,1);
                                t1  = 1433 + 18*randn(ntis,1);
                                t2  = 66 + 3*randn(ntis,1);
                            case {'Caudate_act'}
                                rho = 0.86 + 0*randn(ntis,1);
                                t1  = 1433 + 18*randn(ntis,1);
                                t2  = 111 + 7*randn(ntis,1);
                            case {'Third_ventricle_act','Lateral_ventricle_act','Fourth_ventricle_act','brain_activity','Mamillary_bodies_act','Pineal_gland_act' ,'Cerebral_aqueduct_act'}
                                rho = 1 + 0*randn(ntis,1);
                                t1  = 4200 + 263*randn(ntis,1);
                                t2  = 1500 + 7*randn(ntis,1);
                            case {'lesn_activity'}
                                rho = MRX.Par.tissue.rholesn;
                                t1  = MRX.Par.tissue.t1lesn;
                                t2  = MRX.Par.tissue.t2lesn;
                                
                            otherwise
                                rho = 0;
                                
                        end
                        img(x1,y1,z1) = steady_state_signal(rho,t1,t2,tr,te,a,sequence,ti,t2p);
                    end
                end
            end
            
            
            
            
            % ----------------------------------------------------------------
            %   Update tissue masks
            % ----------------------------------------------------------------
            msk = data;
            
            
            msk1=msk;
            save('msk1.mat','msk1')
        end
        
        %=========================================================================
        % Calc signal intensities using tissue and sequence pars
        function [img,msk] = mapTissueProps_orig(MRX, data)
            
            act = cell2mat(struct2cell(MRX.Par.act));
            tis = fieldnames(MRX.Par.act);
            act = [act;1;18;19]; % xcat
            tis = [tis; 'lesn_activity';'cerebellumWM_act';'cerebellumGM_act'];
            img = single(zeros(size(data)));
            msk = uint8(zeros(size(data)));
            
            for i=1:length(act)
                
                rho = 0;
                % ----------------------------------------------------------------
                %   Select tissue type
                % ----------------------------------------------------------------
                switch char(tis(i))
                    % ------------------------------------------------------------
                    %   Brain xcat
                    % ------------------------------------------------------------
                    
                    case {'eye_activity','skin_activity','lens_activity','sinus_activity','body_activity','salivary_activity'}
                        rho = MRX.Par.tissue.rhoeye;
                        t1  = MRX.Par.tissue.t1eye;
                        t2  = MRX.Par.tissue.t2eye;
                    case {'cortical_bone_activity'}
                        rho = MRX.Par.tissue.rhocorticalbone;
                        t1  = MRX.Par.tissue.t1corticalbone;
                        t2  = MRX.Par.tissue.t2corticalbone;
                    case {'bone_marrow_activity'}
                        rho = MRX.Par.tissue.rhobonemarrow;
                        t1  = MRX.Par.tissue.t1bonemarrow;
                        t2  = MRX.Par.tissue.t2bonemarrow;
                    case {'muscle_activity'}
                        rho = MRX.Par.tissue.rhomuscle;
                        t1  = MRX.Par.tissue.t1muscle;
                        t2  = MRX.Par.tissue.t2muscle;
                    case {'spinal_cord_activity'}
                        rho = MRX.Par.tissue.rhospinalcord;
                        t1  = MRX.Par.tissue.t1spinalcord;
                        t2  = MRX.Par.tissue.t2spinalcord;
                    case {'cartilage_activity'}
                        rho = MRX.Par.tissue.rhocartilage;
                        t1  = MRX.Par.tissue.t1cartilage;
                        t2  = MRX.Par.tissue.t2cartilage;
                    case {'spine_activity'}
                        rho = MRX.Par.tissue.rhospine;
                        t1  = MRX.Par.tissue.t1spine;
                        t2  = MRX.Par.tissue.t2spine;
                    case {'grey_matter_act','Hippocampus_act','Amygdala_act','cerebellumGM_act'}
                        rho = 1*MRX.Par.tissue.rhogreymatter;
                        t1  = MRX.Par.tissue.t1greymatter;
                        t2  = MRX.Par.tissue.t2greymatter;
                    case {'white_matter_act','septum_pellicudum_act','Internal_capsule_act','Anterior_commissure_act','cerebellum_act','cerebellumWM_act','Cerebral_peduncles_act','Superior_cerebellar_peduncle_act','Middle_cerebellar_peduncle_act','Midbrain_act','Tegmentum_of_midbrain_act','Medulla_act' ,'Medullary_pyramids_act','Inferior_olive_act','Superior_colliculus_act','Inferior_colliculus_act','Periacquaductal_grey_outer_act','Periacquaductal_grey_act' ,'Substantia_nigra_act','Pons_act','Fornix_act'}
                        rho = MRX.Par.tissue.rhowhitematter;
                        t1  = MRX.Par.tissue.t1whitematter;
                        t2  = MRX.Par.tissue.t2whitematter;
                    case {'Corpus_Callosum_act'}
                        rho = MRX.Par.tissue.rhocorpuscallosum;
                        t1  = MRX.Par.tissue.t1corpuscallosum;
                        t2  = MRX.Par.tissue.t2corpuscallosum;
                    case {'Putamen_act'}
                        rho = MRX.Par.tissue.rhoputamen;
                        t1  = MRX.Par.tissue.t1putamen;
                        t2  = MRX.Par.tissue.t2putamen;
                    case {'Thalamus_act'}
                        rho = MRX.Par.tissue.rhothalamus;
                        t1  = MRX.Par.tissue.t1thalamus;
                        t2  = MRX.Par.tissue.t2thalamus;
                    case {'Globus_pallidus_act'}
                        rho = MRX.Par.tissue.rhoglobuspallidus;
                        t1  = MRX.Par.tissue.t1globuspallidus;
                        t2  = MRX.Par.tissue.t2globuspallidus;
                    case {'Caudate_act'}
                        rho = MRX.Par.tissue.rhocaudate;
                        t1  = MRX.Par.tissue.t1caudate;
                        t2  = MRX.Par.tissue.t2caudate;
                    case {'Third_ventricle_act','Lateral_ventricle_act','Fourth_ventricle_act','brain_activity','Mamillary_bodies_act','Pineal_gland_act' ,'Cerebral_aqueduct_act'}
                        rho = MRX.Par.tissue.rhobrain;
                        t1  = MRX.Par.tissue.t1brain;
                        t2  = MRX.Par.tissue.t2brain;
                    case {'lesn_activity'}
                        rho = MRX.Par.tissue.rholesn;
                        t1  = MRX.Par.tissue.t1lesn;
                        t2  = MRX.Par.tissue.t2lesn;
                        
                    otherwise
                        rho = 0;
                end
                
                % ----------------------------------------------------------------
                %   Signal model (Balanced SSFP)
                % ----------------------------------------------------------------
                a   = MRX.Par.scan.flip;
                te  = MRX.Par.scan.te;
                tr=MRX.Par.scan.tr;
                ti=MRX.Par.scan.ti;
                t2p=MRX.Par.scan.t2p;
                sequence = MRX.Par.scan.seq ;
                sig = steady_state_signal(rho,t1,t2,tr,te,a,sequence,ti,t2p);
                
                % ----------------------------------------------------------------
                %   Update tissue compartment
                % ----------------------------------------------------------------
                img(find(data(:)==act(i))) = sig;
                
                % ----------------------------------------------------------------
                %   Update tissue masks
                % ----------------------------------------------------------------
                msk(find(data(:)==act(i))) = act(i);
                
            end
            msk1=msk;
            save('msk1.mat','msk1')
        end
        
        
        %=========================================================================
        % Extract k-space segment from fully sampled k-space
        function extractSegment(MRX, segm, img)
            nsegm    = MRX.Par.scan.segments;                   % no. of segments
            kyrange  = ceil((segm-1)*size(img,2)/nsegm)+1:min(ceil(segm*size(img,2)/nsegm),size(img,2));
            temp     = MRXCAT.i2k(img,[1 2 3]);                 % FFT
            MRX.Ksp(:,kyrange,:,:,:) = temp(:,kyrange,:,:,:);   % extract segment
        end
        
        
        %=========================================================================
        % Radial Trajectory Resampling + Calc
        function [img_rad,ksp] = radialResample( MRX, img )
            
            % standard radial or golden angle trajectory
            if strcmpi(MRX.Par.scan.trajectory, 'goldenAngle')
                ga = 1;     % golden angle
            else
                ga = 0;     % stamdard radial
            end
            
            % calculate radial trajectory
            samp = size(img,1);
            % Number of profiles acc. to radial Nyquist & undersampling factor.
            % Cf. MRXCAT.computeBoundingBox for larger FOV in radial than Cartesian.
            prof = round(size(img,1)*2/pi*1/MRX.Par.scan.undersample);  % # prof = 1/1.57*FOVx
            w    = MRXCAT_Brain.getRadWeights2D(samp,prof,0,ga,1);
            k    = MRXCAT_Brain.buildRadTraj2D(samp,prof,0,ga,1);
            % apply trajectory to image
            % Check if NUFFT toolbox (J.Fessler) and wrapper (M.Lustig) are present
            errormsg = ['NUFFT toolbox by J. Fessler and/or NUFFT wrapper ' ...
                'by M. Lustig not found in Matlab path. ' ...
                'Please download and install from ' ...
                ' http://web.eecs.umich.edu/~fessler/code/index.html' ...
                ' and http://www.eecs.berkeley.edu/~mlustig/Software.html'];
            if exist('NUFFT')<2 || exist('nufft')<2
                error(errormsg);
            else % in case of Matlab version < 2011b, exist is not case-sensitive
                [~,f0]=fileparts(which('nufft'));
                [~,f1]=fileparts(which('NUFFT'));
                % case-sensitive check
                if ~strcmp(f0,'nufft') || ~strcmp(f1,'NUFFT')
                    error(errormsg);
                end
            end
            
            %             e    = NUFFT(k,w,[0,0],[samp,samp]);        % Claudio Santelli version (no phase argin)
            e    = NUFFT(k,w,1,[0,0],[samp,samp],1);    % Miki Lustig wrapper (phase=1)
            
            % radially resample for all coil elements
            for k=1:MRX.Par.scan.coils
                ksp(:,:,:,k)     = e*double(img(:,:,:,k));
                img_rad(:,:,:,k) = e'*ksp(:,:,:,k);
            end
        end
        
    end
    
    
    
    methods ( Static )
        
        %=========================================================================
        % Calculate sampling weights for radial trajectory
        function w = getRadWeights2D(no_samples, no_profiles, alt_prof, gafl, normfl)
            %==========================================================================
            % Function which calculates analytically the sampling weights of a radial
            % trajectory.
            %
            % Inputs:
            % -------
            % no_samples:       Number of samples along each projection.
            % no_profiles:      Number of projections per frame and slice.
            % alt_prof:         Alternating profiles flag.
            % gafl:             Golden angle flag.
            % normfl:           Normalization flag.
            %
            % Outputs:
            % --------
            % w:                Radial Weights
            %
            % Function calls:    none
            % ---------------
            %
            % Claudio Santelli, 06/11/2012
            %==========================================================================
            
            % Check input
            %--------------------------------------------------------------------------
            error(nargchk(3,5, nargin))
            if nargin<5, normfl = false; end
            if nargin<4, gafl   = false; end
            
            % Initialize weight and build up kr vector (radial coordinate along spoke)
            %--------------------------------------------------------------------------
            w = []; kr = abs(-floor(no_samples/2):ceil(no_samples/2-1));
            
            if gafl
                % Calculate angle for every spoke and map it onto interval [0,pi]
                %----------------------------------------------------------------------
                phi      = (2*pi/(sqrt(5)+1))*(0:no_profiles-1)+pi/2;
                phi      = ((phi./pi)-floor(phi./pi))*pi;
                
                % Sort spokes according to their angles on the interval [0,pi],
                % calculate relative angular distances to neighboring spokes, and
                % finally, get corresponding total relative angles.
                %----------------------------------------------------------------------
                [phi, I] = sort(phi);
                dPhi1    = [phi(2:end) (phi(1)+pi)]-phi; % Left relative angular distance
                dPhi2    = circshift(dPhi1,[0 1]); % Right relative angular distance
                dPhi     = 0.5*(dPhi1+dPhi2);
                
                % Build up weighting matrix, d, where the i-th column corresponds to
                % the i-th spoke in the sorted order.
                %----------------------------------------------------------------------
                for i=1:no_profiles, w = [w, (dPhi(i)*kr).']; end
                
                % Modify zero point and bring columns of d back into the original order,
                % i.e. i-th column corresponds to the weighting of the i-th spoke.
                %----------------------------------------------------------------------
                w(kr==0,:) = pi/(4*no_profiles);
                w(:,I)     = w;
            else
                % Calculate weight for one spoke and modify zero-point (every spoke has
                % the same weighting).
                %----------------------------------------------------------------------
                w        = (pi/no_profiles)*kr;
                w(kr==0) = pi/(4*no_profiles);
                w        = repmat(w.', [1 no_profiles]);
                if alt_prof
                    w(:,2:2:end) = w(end:-1:1,2:2:end);
                end
            end
            
            if normfl, w = w./max(w(:)); end
            
        end
        
        
        %=========================================================================
        % Create radial trajectory for NUFFT
        function k = buildRadTraj2D(no_samples, no_profiles, alt_prof, gafl, normfl, dim_t, t_offset, dim_z, z_offset)
            %==========================================================================
            % Returns radial trajectory in a complex valued data array. The real part
            % corresponds to the x- and the imaginary part to y-component respectively.
            %
            % Inputs:
            % -------
            % no_samples:   Number of samples along each projection.
            % no_rofiles:   Number of projections per frame and slice.
            % alt_prof:     Alternating profiles flag.
            % gafl:         Golden angle flag.
            % normfl:       Coordinates normalization flag.
            % dim_t:        Number of time frames.
            % t_offset:     Profile offset between two adjacent time frames.
            % dim_z:        Number of slices in z-direction.
            % z_offset:     Profile offset between two adjacent slices.
            %
            % Outputs:
            % --------
            % k:            k-space coordinates.
            %
            % Function calls:    none
            % ---------------
            %
            % Claudio Santelli, 11/06/2012
            %==========================================================================
            
            % Check input
            %--------------------------------------------------------------------------
            error(nargchk(3,9, nargin))
            if nargin<9 || isempty(z_offset), z_offset = 0; end
            if nargin<8 || isempty(dim_z),    dim_z    = 1; end
            if nargin<7 || isempty(t_offset), t_offset = 0; end
            if nargin<6 || isempty(dim_t),    dim_t    = 1; end
            if nargin<5 || isempty(normfl),   normfl   = true; end
            if nargin<4 || isempty(gafl),     gafl     = false; end
            
            % Build trajectory
            %--------------------------------------------------------------------------
            k = [];
            
            % Initial spoke along ky-axis
            k0 = [zeros(1,no_samples); -floor(no_samples/2):ceil(no_samples/2-1)];
            
            % Angle increment
            if gafl
                goldenRatio = (sqrt(5)+1)/2;
                dPhi        = pi/goldenRatio;
            else
                dPhi = pi/no_profiles;
            end
            
            for z=1:dim_z
                for t=1:dim_t
                    for i=1:no_profiles
                        % Update rotation matrix
                        rot_angle = ((i-1)+(t-1)*t_offset+(z-1)*z_offset)*dPhi;
                        if alt_prof && ~mod(i,2)
                            rot_angle = rot_angle+pi;
                        end
                        R = [cos(rot_angle), -sin(rot_angle);
                            sin(rot_angle),  cos(rot_angle)];
                        % Rotate k0 vector accordingly
                        ktmp       = (R*k0).';
                        k(:,i,z,t) = ktmp(:,1)+1i*ktmp(:,2);
                    end
                end
            end
            
            if normfl, k = k./no_samples; end
            
        end
        
    end
    
end
