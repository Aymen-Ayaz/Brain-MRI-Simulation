%=========================================================================
%
% MRXCAT Basic superclass for all MRXCAT classes
%
% The MRXCAT superclass contains all common properties and methods of MRXCAT_CMR_PERF and MRXCAT_CMR_CINE.
%
% NOTE:     Do not call MRXCAT directly, but use MRXCAT_CMR_PERF or MRXCAT_CMR_CINE instead.
%
% PROPERTIES:
%
%           Data                        MRXCAT phantom data (usually empty)
%           Sensitivities               Coil sensitivities  (usually empty)
%           Mask                        XCAT object masks   (usually empty)
%           Par                         Parameter struct with arbitrary sublevels (scan, tissue, contrast/ etc.)
%           Ksp                         K-space data, used for non-Cartesian gridding and k-space segmentation
%           Filename                    Name of the XCAT *.bin file
%           Version                     MRXCAT software version
%
%
% METHODS:  (MRX refers to the MRXCAT instance)
%
%           readLogFile( MRX )          Reads parameters from XCAT log file using XCAT *.bin filename.
%                                       Make sure your log filename is the same as the *.bin name, without *_act_1.bin
%                                       ending and _log extension instead.
%           readImgData( MRX,t )        Reads mask data at time frame t from XCAT .bin file, where t is the time frame
%                                       in the file name, e.g. *_act_15.bin corresponds to t=15.
%  			addNoise( MRX, img )        Adds Gaussian noise to image img.
%                                       The standard deviation of the noise is given by Par.scan.noisestd, which is
%                                       calculated based on desired SNR / CNR.
%           computeBoundingBox( MRX )   Computes image cropping box in x,y,z based on relative values in Par.scan.bbox
%                                       To get the full XCAT image, use [0,1; 0,1; 0,1]
%           calculateCoilMaps( MRX )    Computes voxelwise coil sensitivity maps based on image geometry, and number,
%                                       location and size of coil elements.
%           coilCentres( MRX )          Places coil centres around the torso based on number and size of coil elements.
%           lowPassFilter( MRX, img, msk )
%                                       Applies a low-pass filter on image and mask (img, msk).
%                                       Low-pass filtering is controlled by the flag Par.scan.lowpass and filter
%                                       strength Par.scan.lowpass_str.
%           multiplyCoilMaps( MRX, img, sen )
%                                       Applies coil sensitivities (sen) to image (img) by replicating the image
%                                       Par.scan.coils times and multiplying by the sensitivities.
%           generateFilename( MRX, img )
%                                       Creates MRXCAT output filename (without extension). The filename contains
%                                       MRXCAT parameters from MRX fields (e.g. spatial resolution, # time frames)
%                                       and image parameters derived from img (e.g. matrix size).
%           saveImgData( MRX, img, msk, nois, sen, t )
%                                       Saves image (img), XCAT mask (msk), noise (nois) and sensitivity map (sen)
%										at time point t. For memory handling reasons, the MRXCAT phantom is
%										produced in a loop for each MRXCAT time frame t.
% 			saveParsForRecon( MRX, img )
% 										Saves MRXCAT parameters for recon. This includes the Par struct and
% 										additional image (img) geometry parameters.
%
% HIDDEN METHODS:
%
%			textwaitbar( MRX, t, str_len_in, title )
%										Displays a progress waitbar in the Matlab command window during MRXCAT
%										phantom production.
%
% STATIC, HIDDEN METHODS:
%
%           k2i( img, dims )            inverse FFT (k- -> image space) of img along dimensions dims
%           i2k( img, dims )            FFT (image -> k-space) of img along dimensions dims
%
%
% WEBSITE: 	http://www.biomed.ee.ethz.ch/mrxcat
%
%=========================================================================

%=========================================================================
% VERSION HISTORY:
%			140130LW OO IMPLEMENTATION - v1.0
% 			140321LW VERSION W/ COMMENTS - v1.2
%			150930LW ADDED MRXCAT_Showcase GUI - v1.3
%			170111LW BUGFIX coil calculus - v1.4
%
%
% AUTHORS:	Lukas Wissmann, Claudio Santelli, Sebastian Kozerke
% 			Institute for Biomedical Engineering, University and ETH Zurich
%
%=========================================================================

classdef MRXCAT < handle
    
    properties
        Data;
        Sensitivities;
        Mask;
        Par;
        Ksp;
        Filename;
        Version;
    end
    
    
    
    methods
        function MRX = MRXCAT()
            MRX.Version = 1.4;
        end
        
        
        %=========================================================================
        % Read XCAT log file generated during .bin file production
        function readLogFile( MRX )
            fname = MRX.Filename;
            % --------------------------------------------------------------------
            %   Update number of frames
            % --------------------------------------------------------------------
            MRX.Par.scan.frames_xcat = length(dir([fname(1:end-5),'*.bin']));
            
            % --------------------------------------------------------------------
            %   Trim file name
            % --------------------------------------------------------------------
            fname = [fname(1:end-9),'log'];
            
            % --------------------------------------------------------------------
            %   Open and read log file
            % --------------------------------------------------------------------
            fid = fopen(fname);
            if ( fid ~= -1 )
                while ( ~feof( fid ) )
                    str = fscanf( fid, '%s', 1 );
                    
                    if     strcmp(str, 'array_size'                     ) fscanf(fid,'%s',1); MRX.Par.scan.matrix               = fscanf( fid,'%d',1);
                        % brain
                    elseif strcmp(str, 'eye_activity'                   ) fscanf(fid,'%s',1); MRX.Par.act.eye_activity          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'lens_activity'                  ) fscanf(fid,'%s',1); MRX.Par.act.lens_activity         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'sinus_activity'                 ) fscanf(fid,'%s',1); MRX.Par.act.sinus_activity        = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'cortical_bone_activity'         ) fscanf(fid,'%s',1); MRX.Par.act.cortical_bone_activity= fscanf( fid,'%f',1);
                    elseif strcmp(str, 'bone_marrow_activity'           ) fscanf(fid,'%s',1); MRX.Par.act.bone_marrow_activity  = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'skin_activity'                  ) fscanf(fid,'%s',1); MRX.Par.act.skin_activity        = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'spine_activity'                 ) fscanf(fid,'%s',1); MRX.Par.act.spine_activity       = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'spinal_cord_activity'           ) fscanf(fid,'%s',1); MRX.Par.act.spinal_cord_activity          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'salivary_activity'              ) fscanf(fid,'%s',1); MRX.Par.act.salivary_activity        = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'body_activity'                  ) fscanf(fid,'%s',1); MRX.Par.act.body_activity         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'cartilage_activity'             ) fscanf(fid,'%s',1); MRX.Par.act.cartilage_activity         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'muscle_activity'                ) fscanf(fid,'%s',1); MRX.Par.act.muscle_activity       = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'grey_matter_act'                ) fscanf(fid,'%s',1); MRX.Par.act.grey_matter_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'white_matter_act'               ) fscanf(fid,'%s',1); MRX.Par.act.white_matter_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Internal_capsule_act'           ) fscanf(fid,'%s',1); MRX.Par.act.Internal_capsule_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Anterior_commissure_act'        ) fscanf(fid,'%s',1); MRX.Par.act.Anterior_commissure_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'cerebellum_act'                 ) fscanf(fid,'%s',1); MRX.Par.act.cerebellum_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Cerebral_peduncles_act'         ) fscanf(fid,'%s',1); MRX.Par.act.Cerebral_peduncles_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Superior_cerebellar_peduncle_act') fscanf(fid,'%s',1); MRX.Par.act.Superior_cerebellar_peduncle_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Middle_cerebellar_peduncle_act' ) fscanf(fid,'%s',1); MRX.Par.act.Middle_cerebellar_peduncle_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Corpus_Callosum_act'            ) fscanf(fid,'%s',1); MRX.Par.act.Corpus_Callosum_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Hippocampus_act'                ) fscanf(fid,'%s',1); MRX.Par.act.Hippocampus_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Midbrain_act'                   ) fscanf(fid,'%s',1); MRX.Par.act.Midbrain_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Tegmentum_of_midbrain_act'      ) fscanf(fid,'%s',1); MRX.Par.act.Tegmentum_of_midbrain_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Medulla_act'                    ) fscanf(fid,'%s',1); MRX.Par.act.Medulla_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Medullary_pyramids_act'         ) fscanf(fid,'%s',1); MRX.Par.act.Medullary_pyramids_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Inferior_olive_act'             ) fscanf(fid,'%s',1); MRX.Par.act.Inferior_olive_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Superior_colliculus_act'        ) fscanf(fid,'%s',1); MRX.Par.act.Superior_colliculus_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Inferior_colliculus_act'        ) fscanf(fid,'%s',1); MRX.Par.act.Inferior_colliculus_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Periacquaductal_grey_outer_act' ) fscanf(fid,'%s',1); MRX.Par.act.Periacquaductal_grey_outer_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Periacquaductal_grey_act'       ) fscanf(fid,'%s',1); MRX.Par.act.Periacquaductal_grey_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Substantia_nigra_act'           ) fscanf(fid,'%s',1); MRX.Par.act.Substantia_nigra_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Pons_act'                       ) fscanf(fid,'%s',1); MRX.Par.act.Pons_act         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Amygdala_act'                   ) fscanf(fid,'%s',1); MRX.Par.act.Amygdala_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Fornix_act'                     ) fscanf(fid,'%s',1); MRX.Par.act.Fornix_act         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Putamen_act'                    ) fscanf(fid,'%s',1); MRX.Par.act.Putamen_act         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Thalamus_act'                   ) fscanf(fid,'%s',1); MRX.Par.act.Thalamus_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Globus_pallidus_act'            ) fscanf(fid,'%s',1); MRX.Par.act.Globus_pallidus_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Caudate_act'                    ) fscanf(fid,'%s',1); MRX.Par.act.Caudate_act         = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Third_ventricle_act'            ) fscanf(fid,'%s',1); MRX.Par.act.Third_ventricle_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Lateral_ventricle_act'          ) fscanf(fid,'%s',1); MRX.Par.act.Lateral_ventricle_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Fourth_ventricle_act'           ) fscanf(fid,'%s',1); MRX.Par.act.Fourth_ventricle_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'brain_activity'                 ) fscanf(fid,'%s',1); MRX.Par.act.brain_activity          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Mamillary_bodies_act'           ) fscanf(fid,'%s',1); MRX.Par.act.Mamillary_bodies_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Pineal_gland_act'               ) fscanf(fid,'%s',1); MRX.Par.act.Pineal_gland_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'Cerebral_aqueduct_act'          ) fscanf(fid,'%s',1); MRX.Par.act.Cerebral_aqueduct_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'septum_pellicudum_act'          ) fscanf(fid,'%s',1); MRX.Par.act.septum_pellicudum_act          = fscanf( fid,'%f',1);
                    elseif strcmp(str, 'pixel'                          ) str = fscanf( fid, '%s', 1 );
                        if strcmp(str, 'width'                          ) fscanf(fid,'%s',1); MRX.Par.scan.rx_cm                = fscanf( fid,'%f',1); end
                    elseif strcmp(str, 'slice'                          ) str = fscanf( fid, '%s', 1 );
                        if strcmp(str, 'width'                          ) fscanf(fid,'%s',1); MRX.Par.scan.rz_cm                = fscanf( fid,'%f',1); end
                    elseif strcmp(str, 'starting'                       ) str = fscanf( fid, '%s', 1 ); str1 = fscanf( fid, '%s', 1 );
                        if strcmp(str, 'slice') && strcmp(str1, 'number') fscanf(fid,'%s',1); sl_start                          = fscanf( fid,'%f',1); end
                    elseif strcmp(str, 'ending'                         ) str = fscanf( fid, '%s', 1 ); str1 = fscanf( fid, '%s', 1 );
                        if strcmp(str, 'slice') && strcmp(str1, 'number') fscanf(fid,'%s',1); sl_end                            = fscanf( fid,'%f',1); end
                    elseif strcmp(str, '==>Total'                        ) str = fscanf( fid, '%s', 1 ); str1 = fscanf( fid, '%s', 1 );
                        if strcmp(str, 'Output') && strcmp(str1, 'Period') fscanf(fid,'%s',1); MRX.Par.scan.scan_dur            = fscanf( fid,'%f',1); end
                    elseif strcmp(str, 'beating'                        ) str = fscanf( fid, '%s', 1 ); str1 = fscanf( fid, '%s', 1 );
                        if strcmp(str, 'heart') && strcmp(str1, 'period') fscanf(fid,'%s',1); MRX.Par.scan.heartbeat_length     = fscanf( fid,'%f',1); end
                    elseif strcmp(str, 'Respiration motion and beating heart motions included');              MRX.Par.scan.resp = 1; %will be overwritten if Resp option set outside function
                    elseif strcmp(str, 'Respiration motion included only'                     );              MRX.Par.scan.resp = 1; %will be overwritten if Resp option set outside function
                    elseif strcmp(str, 'Beating heart motion included only'                   );              MRX.Par.scan.resp = 0; %will be overwritten if Resp option set outside function
                    end
                end
                fclose( fid );
                if( exist('sl_start', 'var') && exist('sl_end', 'var') ) % number of slices
                    MRX.Par.scan.slices = sl_end-sl_start+1;
                end
                
                fid = fopen(fname);
                if ( fid ~= -1 )
                    while ( ~feof( fid ) )
                        str = fgets(fid);
                        if     strfind(str, 'Respiration motion and beating heart motions included'); MRX.Par.scan.resp = 1; %will be overwritten if Resp option set outside function
                        elseif strfind(str, 'Respiration motion included only'                     ); MRX.Par.scan.resp = 1; %will be overwritten if Resp option set outside function
                        elseif strfind(str, 'Beating heart motion included only'                   ); MRX.Par.scan.resp = 0; %will be overwritten if Resp option set outside function
                        end
                    end
                end
                fclose( fid );
                
                fprintf('\nPhantom information:\n');
                fprintf('  matrix     :%4d\n',MRX.Par.scan.matrix);
                fprintf('  frames     :%4d\n',MRX.Par.scan.frames_xcat);
                fprintf('  resolution :   %1.2f x %1.2f x %1.2f cm3\n\n', MRX.Par.scan.rx_cm, MRX.Par.scan.rx_cm, MRX.Par.scan.rz_cm);
            else
                error('Cannot read XCAT log file. Aborting ... ');
            end
        end
        
        
        %=========================================================================
        % Read data from XCAT .bin file
        function img = readImgData(MRX,t)
            
            fname = MRX.Filename;
            
            % --------------------------------------------------------------------
            %   in case of more output than input frames, restart at first XCAT frame
            %   after last frame is reached.
            %   !!! Make sure breathing/heart rate fits with number of frames in XCAT parameter file !!!
            %   --> Only done for perfusion, but probably the 2nd condition <--
            %   --> (strcmpi(class(MRX), 'MRXCAT_CMR_PERF')) can be removed <--
            % --------------------------------------------------------------------
            if t > MRX.Par.scan.frames_xcat && strcmpi(class(MRX), 'MRXCAT_CMR_PERF')
                t = mod(t,MRX.Par.scan.frames_xcat);
            end
            
            % --------------------------------------------------------------------
            %   Trim file name
            %   Use only the 1st .bin file for breathhold perfusion (1 heart phase)
            % --------------------------------------------------------------------
            if MRX.Par.scan.resp || strcmpi(class(MRX), 'MRXCAT_CMR_CINE')
                fname = [fname(1:end-5),sprintf('%d',t),'.bin'];
            end
            
            % --------------------------------------------------------------------
            %   Open data file
            % --------------------------------------------------------------------
            fid = fopen(fname);
            
            if ( fid ~= -1 )
                img = fread(fid,'single');
                img = reshape(img,MRX.Par.scan.matrix,MRX.Par.scan.matrix,[]);
                
                [xdim,ydim,zdim] = MRX.computeBoundingBox;
                
                % crop image to bounding box
                img  = img(xdim,ydim,zdim);
                
                fclose( fid );
            else
                fprintf([fname,' cannot be read\n']);
            end
        end
        
        
        % Add Gaussian noise to image
        function [KSPACE, nois] = addNoise(MRX, KSPACE)
            stdev = MRX.Par.scan.noisestd;
            nois  = stdev*complex(randn(size(KSPACE)),randn(size(KSPACE)));
            KSPACE   = KSPACE + nois;
            
        end
        
        function saveNifti(MRX, KSPACE, msk)
            
            filename = [MRX.Par.scan.name,'_','Res_',num2str(MRX.Par.scan.lowpass_str(1)),'_',num2str(MRX.Par.scan.lowpass_str(2)),'_',num2str(MRX.Par.scan.lowpass_str(3))];
            image = abs(ifft2c(KSPACE));
            %%      scale to 12 bit range 0-4096
            inmin = min(image(find(image)));
            inmax = max(image(find(image)));
            image = rescale(image,0,4095,'InputMin',inmin,'InputMax',inmax);
            scale = (MRX.Par.scan.resizeimage / size(image,1));
            image = imresize(image,scale,'bicubic');
            msk = imresize(msk,scale,'nearest');
            
            image(msk==0)=0;
            niftiwrite(uint16((image)),['.\',filename,'_img.nii'],'Compressed',true)
            niftiwrite(uint16((msk)),['.\',filename,'_msk.nii'],'Compressed',true)
            info = niftiinfo(['.\',filename,'_img.nii','.gz']);
            infomsk = niftiinfo(['.\',filename,'_msk.nii','.gz']);
            
            info.PixelDimensions=[ MRX.Par.scan.lowpass_str(1)/scale  MRX.Par.scan.lowpass_str(2)/scale  MRX.Par.scan.lowpass_str(3)];
            info.SpaceUnits ='Millimeter';
            info.MultiplicativeScaling = 1;
            info.Transform = 'Qform';
            info.Qfactor = 1;
            % An affine3d object encapsulates a 3-D affine geometric transformation.
            A = [-0.566895215624018,-0.756577734457881,-0.308039938473936,0;0.474734319495109,7.85402026108167e-08,-0.873668119690399,0;5.34858966844201,-5.19094362154612,2.90631947255718,0;-80.4087371826172,202.598815917969,260.606781005859,1];
            tform = affine3d(A);
            info.Transform = tform;
            info.raw.scl_slope = 1;
            info.raw.qform_code = 1;
            info.raw.sform_code = 0;
            info.raw.qform_code_name = 'NIFTI_XFORM_SCANNER_ANAT';
            info.raw.xyzt_units = 2;
            info.raw.quatern_b = -0.129011020064354;
            info.raw.quatern_c = 0.549282431602478;
            info.raw.quatern_d = -0.695238173007965;
            info.raw.qoffset_x = -80.408737182617190;
            info.raw.qoffset_y = 2.025988159179688e+02;
            info.raw.qoffset_z = 2.606067810058594e+02;
            
            infomsk.PixelDimensions = info.PixelDimensions; infomsk.SpaceUnits = info.SpaceUnits;
            infomsk.MultiplicativeScaling = info.MultiplicativeScaling;
            infomsk.Transform=info.Transform;
            infomsk.Qfactor = info.Qfactor;
            infomsk.Transform = tform;
            infomsk.raw.scl_slope = 1;
            infomsk.raw.qform_code = 1;
            infomsk.raw.sform_code = 0;
            infomsk.raw.qform_code_name = 'NIFTI_XFORM_SCANNER_ANAT';
            infomsk.raw.xyzt_units = 2;
            infomsk.raw.quatern_b = -0.129011020064354;
            infomsk.raw.quatern_c = 0.549282431602478;
            infomsk.raw.quatern_d = -0.695238173007965;
            infomsk.raw.qoffset_x = -80.408737182617190;
            infomsk.raw.qoffset_y = 2.025988159179688e+02;
            infomsk.raw.qoffset_z = 2.606067810058594e+02;
            niftiwrite(uint16((image)),['.\',filename,'_img.nii'],info,'Compressed',true)
            niftiwrite(uint16((msk)),['.\',filename,'_msk.nii'],infomsk,'Compressed',true)
            
        end
        
        %=========================================================================
        % Compute bounding box based on MRX.Par.scan.bbox (0 < parameters < 1)
        function [xdim,ydim,zdim] = computeBoundingBox( MRX )
            
            xdim = fix((MRX.Par.scan.bbox(1)*MRX.Par.scan.matrix+1):MRX.Par.scan.bbox(4)*MRX.Par.scan.matrix);
            ydim = fix((MRX.Par.scan.bbox(2)*MRX.Par.scan.matrix+1):MRX.Par.scan.bbox(5)*MRX.Par.scan.matrix);
            zdim = fix((MRX.Par.scan.bbox(3)*MRX.Par.scan.slices+1):MRX.Par.scan.bbox(6)*MRX.Par.scan.slices);
            
            if( isempty(xdim) || isempty(ydim) || isempty(zdim) )
                error('Bounding box dimension empty! Increase BoundingBox size or XCAT size');
            end
            
            % adapt ydim to multiple of possible k-t factors
            f = factor(MRX.Par.scan.frames);
            if length(f) < 2 || f(1)~=2
                warning('Number of time points should be multiple of 2 and no prime number for k-t studies!');
            end
            if isfield(MRX.Par.scan, 'adaptPhaseEncDim') && MRX.Par.scan.adaptPhaseEncDim
                rdy   = rem(length(ydim),MRX.Par.scan.frames);
                ydim = ydim(1)+ceil(rdy/2):ydim(end)-floor(rdy/2);
            end
            
            if( length(xdim) ~= length(ydim) && ~strcmpi(MRX.Par.scan.trajectory, 'Cartesian') )
                % square FOV for radial resampling
                dly  = length(ydim)-length(xdim);
                xdim = xdim(1)-floor(dly/2):xdim(end)+ceil(dly/2);
                % adjust FOV for radial (radial Nyquist requirement)
                samp=round(length(xdim)*pi/2);
                xdimplus = (samp-length(xdim))/2;
                xdim = xdim(1)-floor(xdimplus):xdim(end)+floor(xdimplus);
                ydim = ydim(1)-floor(xdimplus):ydim(end)+floor(xdimplus);
                % !!! add check for negative indices here !!!
            end
        end
        
        
        %=========================================================================
        % Calc voxelwise coil maps
        function sen = calculateCoilMaps( MRX )
            
            % -----------------------------------------------------------------
            % Bounding box indices
            % -----------------------------------------------------------------
            [xdim,ydim,zdim] = MRX.computeBoundingBox;
            
            if MRX.Par.scan.coils>1
                
                nc = MRX.Par.scan.coils;
                rx = MRX.Par.scan.rx_cm*10; % [mm] voxel size
                rz = MRX.Par.scan.rz_cm*10; % [mm] slice width
                
                % -----------------------------------------------------------------
                % Coil centre locations, coil radius
                % -----------------------------------------------------------------
                [cc,R] = MRX.coilCentres; %nc coils, 450mm radius (??), 600 mm coil array length => replace/remove
                % -----------------------------------------------------------------
                % Define rotation
                % -----------------------------------------------------------------
                a = 0; %MRX.Par.scan.rotation(1); % (needs debugging for correct rotation!)
                b = 0; %MRX.Par.scan.rotation(2); % (needs debugging for correct rotation!)
                c = 0; %MRX.Par.scan.rotation(3); % verified!
                rotx = [1,0,0; 0,cos(a),-sin(a); 0,sin(a),cos(a)];
                roty = [cos(b),0,sin(b); 0,1,0; -sin(b),0,cos(b)];
                rotz = [cos(c),-sin(c),0; sin(c),cos(c),0; 0,0,1];
                invr = inv(rotx*rotz*roty);
                
                % -----------------------------------------------------------------
                % angles for integration
                % -----------------------------------------------------------------
                angles = 60;
                dtheta = 2*pi/angles;
                theta = -pi:dtheta:pi-dtheta;
                
                % -----------------------------------------------------------------
                % voxel coordinates with origin in image centre
                % -----------------------------------------------------------------
                x = 0:1:length(xdim)-1;
                x = x-x(end)/2;
                y = 0:1:length(ydim)-1;
                y = y-y(end)/2;
                z = 0:1:length(zdim)-1;
                z = z-z(end)/2;
                
                
                % -----------------------------------------------------------------
                % convert voxel coordinats to mm
                % -----------------------------------------------------------------
                x = x*rx;
                y = y*rx;
                z = z*rz;
                
                
                % -----------------------------------------------------------------
                % calculate sensitivity for each coil
                % -----------------------------------------------------------------
                
                sen = ones(numel(x),numel(y),numel(z),nc);
                
                % reset contrast concentrations
                MRX.Par.contrast.cm.lv=0;
                MRX.Par.contrast.cm.rv=0;
                MRX.Par.contrast.cm.la=0;
                MRX.Par.contrast.cm.ra=0;
                MRX.Par.contrast.ca.lv=0;
                MRX.Par.contrast.ca.rv=0;
                MRX.Par.contrast.ca.la=0;
                MRX.Par.contrast.ca.ra=0;
                % normalize sens
                msen    = mean(sen(sen~=0));
                sen     = 1/msen*(real(sen)+1i*imag(sen));
            else
                sen = ones(length(xdim),length(ydim),length(zdim),1);
            end
        end
        
        
        %=========================================================================
        % Calc coil centres on circles around body
        function [cc, rcoil] = coilCentres( MRX )
            % Generate coil centres on cylinder surface for sensitivity-encoded MR simulation
            %
            % Input:    MRX                         MRXCAT object
            %
            % Parameters used:
            %           MRX.Par.scan.coils          number of coils
            %           MRX.Par.scan.coildist       body radius = radius of coil array in ('AP','RL') plane
            %           MRX.Par.scan.coilsperrow    max. number of coils per coil array row / "ring" / circle
            %
            % Output:   cc [mm]         			coil centre xyz coordinates as [ncoils, 3] matrix
            %           rcoil [mm]					coil radius
            %
            
            
            nc        = MRX.Par.scan.coils;                    % No of coils
            rbody_mm  = MRX.Par.scan.coildist;                 % distance of coil centres from image centre
            c_per_row = MRX.Par.scan.coilsperrow;              % No of coils per full circle on cylinder surface
            % Set this parameter to max. 8
            
            % Define specific angular distributions (low #coils)
            % For >6 coils, equal distribution from 0 to 2 PI is assumed
            angles2 = [150, 210] * pi/180;                   % [rad] 2 anterior coil elements
            angles3 = [130, 180, 230] * pi/180;              % [rad] 3 anterior coil elements
            angles4 = [150, 210, 330, 30] * pi/180;          % [rad] 2 anterior, 2 posterior
            angles5 = [130, 180, 230, 330, 30] * pi/180;     % [rad] 3 anterior, 2 posterior
            angles6 = [130, 180, 230, 310, 0, 50] * pi/180;  % [rad] 3 anterior, 3 posterior
            
            % Determine number of coil rows ("rings") and number of coil elements per ring
            if rem(nc,c_per_row)
                compl_rows = max(floor((nc-c_per_row)/c_per_row),0);
                remc  = nc-compl_rows*c_per_row;
                % divide into 2 rows with as equal partition as possible
                if remc>c_per_row
                    remc = [floor(remc/2),ceil(remc/2)];
                end
            else
                compl_rows = nc/c_per_row;
                remc  = [];
            end
            rings = compl_rows + numel(remc);
            switch numel(remc)
                case 0
                    nc_ring(1:rings)    = c_per_row;
                case 1
                    nc_ring(1)          = remc;
                    nc_ring(2:rings)    = c_per_row;
                case 2
                    nc_ring(1)          = remc(1);
                    nc_ring(2:rings-1)  = c_per_row;
                    nc_ring(rings)        = remc(2);
            end
            
            % z coordinate based on coil element radius
            if any(nc_ring>6)
                minang = 2*pi/max(nc_ring);
            else
                minang = 50*pi/180;
            end
            rcoil = 1/2*rbody_mm*minang;
            z     = -(rings-1)*rcoil:2*rcoil:(rings-1)*rcoil;
            
            % x-y coordinates
            cc    = zeros(nc,3);
            lctr  = 0; % index counter
            for k=1:rings
                switch nc_ring(k)
                    case 1
                        [x{k}, y{k}] = pol2cart(0, rbody_mm);
                    case 2
                        [x{k}, y{k}] = pol2cart(angles2, rbody_mm);
                    case 3
                        [x{k}, y{k}] = pol2cart(angles3, rbody_mm);
                    case 4
                        [x{k}, y{k}] = pol2cart(angles4, rbody_mm);
                    case 5
                        [x{k}, y{k}] = pol2cart(angles5, rbody_mm);
                    case 6
                        [x{k}, y{k}] = pol2cart(angles6, rbody_mm);
                    otherwise % equally distribute on circle
                        c = nc_ring(k);
                        [x{k}, y{k}] = pol2cart(linspace(0,2*pi*(c-1)/c,c),rbody_mm);
                end
                % assign coil centre positions (x,y,z)
                cc(lctr+1:lctr+length(x{k}),1:3) = [x{k}', y{k}', z(k)*ones(length(x{k}),1)];
                lctr = lctr + length(x{k});
            end
            
        end
        
        
        %=========================================================================
        % Low-pass filter image and mask
        function [img, msk , KSPACE] = lowPassFilter( MRX, img, msk )
            
            if MRX.Par.scan.lowpass
                
                %========K-space Sampling with desired resolution
                
                resx=MRX.Par.scan.lowpass_str(1);    %field of view in mm
                resy=MRX.Par.scan.lowpass_str(2);
                resz=MRX.Par.scan.lowpass_str(3);    %resolution in mm
                
                XCATsize = size(img);
                %field of view for input high resolution XCAT
                XCATFOV = [MRX.Par.scan.rx_cm *10* XCATsize(1),MRX.Par.scan.rx_cm *10* XCATsize(2), MRX.Par.scan.rz_cm *10* XCATsize(3)];
                
                %=========% matrix dimention at each direction should be even for windowing the center with Turkey window
                nx = round(XCATFOV(1) / resx) ;
                if mod((nx),2) == 1
                    nx = nx + 1;
                end
                ny = round(XCATFOV(2) / resy) ;
                if mod((ny),2) == 1
                    ny = ny + 1;
                end
                nz = round(XCATFOV(3) / resz) ;
                
                %=========% make a TukeyWindow with desired resolution to sample from
                %the high resolution k-space derived from high resolution XCAT labels
                Lx = nx;
                Ly = ny;
                Lz = nz;
                t25x = tukeywin(Lx,0.5);
                t25y = tukeywin(Ly,0.5);
                [maskt25x,maskt25y]=meshgrid(t25x,t25y);
                Bt25xyz=maskt25x.*maskt25y;
                
                Bt25 = padarray(Bt25xyz,[(XCATsize(1)-(nx))/2,(XCATsize(2)-(ny))/2],'both');
                
                
                
                %instead of slice averaging we downsample in Z direction from the high
                %resolution image in z direction to desired slice number
                
                imglow = imresize3(img,[size(img,1),size(img,2),nz],'cubic');
                
                KSPACE = zeros (size(Bt25,1),size(Bt25,2),size(imglow,3));
                
                %=========loop over coils and slices
                for j=1:size(imglow,4)
                    for i=1:size(imglow,3)
                        
                        XCATKspace = fft2c(squeeze(imglow(:,:,i,j)));
                        KSPACEzeropad = Bt25 .* XCATKspace;
                        KSPACE(:,:,i) = KSPACEzeropad;
                        
                        
                        imglow(:,:,i,j) = ifft2c(KSPACEzeropad);
                        
                    end
                end
                img = imglow;
                msk = imresize3(msk,[size(msk,1) size(msk,2) size(KSPACE,3)],'nearest');
                clear A FOV resx resy ns imgks img2;
                
                
            end
            
        end
        
        %=========================================================================
        % Replicate image and multiply w/ coil sensitivity maps
        
        
        function img = multiplyCoilMaps( MRX, img, sen)
            sen = ones(size(img));
            img = repmat(img,[1 1 1 MRX.Par.scan.coils]);
            img = img.*sen;
        end
        
        
        %=========================================================================
        % Create file name for saving MRXCAT phantom data
        function filename = generateFilename(MRX, img)
            % determine short filename
            [p,f,~] =   fileparts(MRX.Filename);
            fname   =   [p,filesep,f(1:4)];
            
            % create different strings to be concatenated in a filename string
            if MRX.Par.scan.resp==0
                bhstr = '_bh';
            else
                bhstr = '_fb';
            end
            resolstr = ['_',num2str(MRX.Par.scan.rx_cm*10),'x',num2str(MRX.Par.scan.rx_cm*10),'x',num2str(MRX.Par.scan.rz_cm*10),'mm'];
            matstr = ['_',num2str(size(img,1)),'x',num2str(size(img,2)),'x',num2str(size(img,3)),'x',num2str(MRX.Par.scan.frames)];
            if MRX.Par.scan.coils>1
                matstr = [matstr,'x',num2str(MRX.Par.scan.coils)];
            end
            snrstr      = ['_snr',num2str(MRX.Par.scan.snr)];
            fastr       = ['_fa',num2str(round(MRX.Par.scan.flip*180/pi))];
            % case separation between perf and other phantom types
            if strcmpi(class(MRX), 'MRXCAT_CMR_PERF') % perfusion
                dosestr     = ['_dose',num2str(MRX.Par.contrast.dose)];
                tshiftstr   = ['_tshift',num2str(MRX.Par.contrast.tshift)];
                if MRX.Par.contrast.rs == 2
                    restStress = 'Stress';
                elseif MRX.Par.contrast.rs == 1
                    restStress = 'Rest';
                else
                    restStress = '';
                end
                filename = [fname,restStress,resolstr,matstr,snrstr,dosestr,fastr,tshiftstr,bhstr];
            else % cine, ...
                filename = [fname,resolstr,matstr,snrstr,fastr,bhstr];
            end
        end
        
        
        %=========================================================================
        % Save image data (cpx), mask (msk), noise (noi) and sensitivities (sen)
        function saveImgData( MRX, img, msk, nois, sen, t)
            
            % --------------------------------------------------------------------
            %   Generate filename and append extensions
            % --------------------------------------------------------------------
            fname = MRX.generateFilename( img );
            fimg = [fname,'.cpx'];
            fmsk = [fname,'.msk'];
            fnoi = [fname,'.noi'];
            fsen = [fname,'.sen'];
            
            % --------------------------------------------------------------------
            %   Write/append data file
            % --------------------------------------------------------------------
            if t==1
                fprintf('\nOutput file information:\n');
                fprintf('  matrix :%4d x%4d x%4d\n',size(img,1),size(img,2),size(img,3));
                fprintf('  frames :%4d\n',MRX.Par.scan.frames);
                fprintf('  coils  :%4d\n',MRX.Par.scan.coils);
                
                fidimg = fopen(fimg,'w');
                fidmsk = fopen(fmsk,'w');
                fidnoi = fopen(fnoi,'w');
                fidsen = fopen(fsen,'w');
            else
                fidimg = fopen(fimg,'a');
                fidmsk = fopen(fmsk,'a');
                fidnoi = fopen(fnoi,'a');
                fidsen = fopen(fsen,'a');
            end
            
            if ( fidimg ~= -1 )
                dim = size(img);
                if length(dim) < 3, dim(3) = 1; end
                img  = [real(img),imag(img)];
                img  = reshape(img,dim(1),dim(2),2,dim(3),MRX.Par.scan.coils);
                img  = permute(img,[3 1 2 4 5]);
                nois = [real(nois),imag(nois)];
                nois = reshape(nois,dim(1),dim(2),2,dim(3),MRX.Par.scan.coils);
                nois = permute(nois,[3 1 2 4 5]);
                sen  = [real(sen),imag(sen)];
                sen  = reshape(sen,dim(1),dim(2),2,dim(3),MRX.Par.scan.coils);
                sen  = permute(sen,[3 1 2 4 5]);
                
                fwrite(fidimg,single(img),'single');
                fwrite(fidmsk,uint8(msk),'uint8');
                fwrite(fidnoi,single(nois),'single');
                fwrite(fidsen,single(sen),'single');
                
                fclose(fidimg);
                fclose(fidmsk);
                fclose(fidnoi);
                fclose(fidsen);
            else
                fprintf([fimg,' cannot be written\n']);
            end
        end
        
        
        %=========================================================================
        % Save MRXCAT parameters for recon
        function saveParsForRecon( MRX, img )
            % append recon parameters to parameter struct and save it to file
            
            % --------------------------------------------------------------------
            %   Generate filename
            % --------------------------------------------------------------------
            fname   = MRX.generateFilename( img);
            
            fpar = [fname,'_par.mat'];
            xdim = size(img,1);
            ydim = size(img,2);
            zdim = size(img,3);
            
            MRX.Par.acq_ovs        = ones(1,4);
            MRX.Par.acq_matrix     = [xdim ydim zdim 1];        % ACQ matrix
            MRX.Par.acq_ne0        = xdim;
            MRX.Par.acq_ne1        = ydim;
            MRX.Par.acq_ne2        = zdim;
            MRX.Par.rec_matrix     = MRX.Par.acq_matrix(1:3);   % RECON marix
            MRX.Par.rec_file_name  = [fname,'.cpx'];
            MRX.Par.proto_array    = [10 11 7 1 0 0 100 0];     % [kt-Factor, ky_train, kz_train, kt-PCA flag, plug in
            %  training, phase recon, #PCA components, comp-based flag]
            MRX.Par.ncoils         = MRX.Par.scan.coils;        % duplicate for ReconKtData.m
            
            % save par struct
            Par = MRX.Par;
            save( fpar, 'Par');
        end
        
    end
    
    
    
    methods ( Hidden )
        %=========================================================================
        % Text waitbar (replaces graphical waitbar)
        %=========================================================================
        function strlen = textwaitbar( MRX, t, str_len_in, title )
            % MRX:			MRXCAT instance
            % T:         	between 0 and # MRXCAT frames
            % STR_LEN_IN:   length of previous waitbar
            % TITLE:        optional
            %
            % HOWTO USE:    E.g. in a for loop, write
            %
            %               strlen = 0;
            %               for k=1:N
            %                   strlen = textwaitbar(MRX, k, strlen, 'Calculating task name');
            %                   .... (tasks WITHOUT command line output)
            %               end
            
            % L.Wissmann, 12 Apr 2013
            
            frac    = t/MRX.Par.scan.frames;
            
            len     = 40;
            numblk  = round(frac*len);
            strlen  = 0;
            
            if nargin > 1 && ~isempty(str_len_in) && str_len_in > 0
                for k=1:str_len_in
                    fprintf('\b');
                end
            else
                fprintf('\n------------------------------------------\n');
            end
            if nargin > 2
                fprintf('%s\n\n',title);
                strlen = strlen + length(title)+2;
            else
                stdstr = 'Computing ...\n\n';
                fprintf(stdstr);
                strlen = strlen + length(stdstr)-2;
            end
            fprintf('[');
            strlen = strlen + 1;
            for k=1:numblk
                fprintf('o');
                strlen = strlen+1;
            end
            for k=1:len-numblk
                fprintf(' ');
                strlen = strlen+1;
            end
            fprintf(']');
            strlen = strlen + 1;
            
            if frac >= 1
                fprintf('\n\nDone.\n------------------------------------------\n');
            end
            
        end
        %=========================================================================
        
    end
    
    
    
    methods ( Hidden, Static )
        %=========================================================================
        % FFT (image -> k-space)
        %=========================================================================
        function img = i2k(img, dims)
            
            dim_img = size(img);
            
            if nargin < 2
                factor  = prod(dim_img);
                img = 1/sqrt(factor)*fftshift(fftn(ifftshift( img )));
            else
                for dim = dims
                    if size(img,dim)>1
                        img = 1/sqrt(dim_img(dim))*fftshift(fft(ifftshift( img, dim ),[],dim),dim);
                    end
                end
            end
            
        end
        
        
        %=========================================================================
        % inverse FFT (k- -> image space)
        %=========================================================================
        function img = k2i(img, dims)
            
            dim_img = size(img);
            
            if nargin < 2
                factor  = prod(dim_img);
                img = sqrt(factor)*fftshift( ifftn(ifftshift(img)) );
            else
                for dim = dims
                    if size(img,dim)>1
                        img = sqrt(dim_img(dim))*fftshift( ifft(ifftshift(img,dim),[],dim), dim);
                    end
                end
            end
            
        end
        
    end
    
end
