%read_pillar_bin Read PILLAR binary data.
%   [META, DATA] = read_pillar_bin(FILENAME) reads binary data from the
%   specified file by the string FILENAME. The return value META stores
%   file meta information, parameters used to generate this file. and DATA
%   stores the raw coordinates, drift-corrected coordinates and
%   corresponding deflections if existed.
function [file_header, data] = read_pillar_bin(fname)
    file_header = struct;
    data = struct;
    
    %% open the file
    fileID = fopen(fname);
    readmode = 'b';
    
    %% read header
    fileversion1 = 589450;  %file version control
    fileversion2 = 589451;  %since v1.1.3 if grid restoration applied
    file_header.version = fread(fileID,1,'int',readmode);
    file_header.npillars = 0;
    file_header.nframes = 0;
    file_header.suc_load = false;
    if(file_header.version==fileversion1 || file_header.version==fileversion2)
        file_header.npillars    = fread(fileID,1,'int',readmode);
        file_header.nframes     = fread(fileID,1,'int',readmode);
        file_header.pixel_size  = fread(fileID,1,'double',readmode);
        file_header.diameter    = fread(fileID,1,'double',readmode);
        file_header.spacing     = fread(fileID,1,'double',readmode);
        file_header.oblique     = fread(fileID,1,'double',readmode);
        file_header.grid_angle  = fread(fileID,1,'double',readmode);
        file_header.sigma_PSF   = fread(fileID,1,'double',readmode);
        file_header.catch_radius= fread(fileID,1,'double',readmode);
        file_header.kernel_w    = fread(fileID,1,'int',readmode);
        file_header.box_limit   = fread(fileID,1,'int',readmode);
        file_header.dark_pillar     = fread(fileID,1,'int8',readmode);
        file_header.use_minimum_std = fread(fileID,1,'int8',readmode);
        file_header.use_metric_CG   = fread(fileID,1,'int8',readmode);
        file_header.use_enhancer    = fread(fileID,1,'int8',readmode);  
        if(file_header.version==fileversion2)
            file_header.use_fft     = fread(fileID,1,'int8',readmode);  
        else
            file_header.use_fft = 0;
        end
        file_header.enhance_loaded_suc = true;
        if(file_header.use_enhancer>0)
            file_header.enhancer.normalized_cc = fread(fileID,1,'int8',readmode);
            file_header.enhancer.use_gauss_psf = fread(fileID,1,'int8',readmode);
            file_header.enhancer.gauss_psf_dark= fread(fileID,1,'int8',readmode);
            file_header.enhancer.gauss_psf_radius= fread(fileID,1,'int',readmode);
            file_header.enhancer.gauss_psf_sigma = fread(fileID,1,'double',readmode);
            file_header.enhancer.psf_w = fread(fileID,1,'int',readmode);
            file_header.enhancer.psf_h = fread(fileID,1,'int',readmode);        
            if(file_header.enhancer.psf_w>0 && file_header.enhancer.psf_h>0)
                file_header.enhancer.PSF = fread(fileID,[file_header.enhancer.psf_w, file_header.enhancer.psf_h],'double',readmode);            
            else
                file_header.enhance_loaded_suc = false;
            end
        end
        if(file_header.use_fft>0)
            file_header.FFT.off_center_radius               = fread(fileID,1,'int',readmode);
            file_header.FFT.center_radius                   = fread(fileID,1,'int',readmode);
            file_header.FFT.start_off_center_radius         = fread(fileID,1,'int',readmode);
            file_header.FFT.end_off_center_center_radius    = fread(fileID,1,'int',readmode);
            file_header.FFT.num_points                      = fread(fileID,1,'int',readmode);
            file_header.FFT.points                          = fread(fileID,[2,file_header.FFT.num_points],'int',readmode);
        end
        if(file_header.enhance_loaded_suc)
            startXY= fread(fileID,[2,file_header.npillars],'int',readmode); 
            file_header.startX = startXY(1,:);
            file_header.startY = startXY(2,:);
            file_header.suc_load = true;
        end
    elseif(mod(file_header.version,2)==0)
        file_header.npillars = version/2;
        file_header.nframes = fread(fileID,1,'int',readmode);
        file_header.suc_load = true;
    end

    %% check whether the file header is loaded successfully.
    if(file_header.suc_load && file_header.npillars>0 && file_header.nframes>0)
        %% read out the tracks
        ncols = file_header.npillars*2;
        nframes = file_header.nframes;
        A = fread(fileID,[ncols, nframes],'double',readmode);
        B = fread(fileID,[ncols, nframes],'double',readmode);
        if(file_header.version==fileversion2)
            %% deflections 
            data.DX=B(1:2:ncols,:)';
            data.DY=B(2:2:ncols,:)';
            %% drift corrected data
            B = fread(fileID,[ncols, nframes],'double',readmode);
        end

        %% raw data
        data.tracksX=A(1:2:ncols,:)';
        data.tracksY=A(2:2:ncols,:)';

        %% drift correct data
        if(numel(B)>0)
            data.drift_correctX=B(1:2:ncols,:)';
            data.drift_correctY=B(2:2:ncols,:)';
            % drift information
            driftXY = fread(fileID,[2, nframes],'double',readmode);
            data.driftX = driftXY(1,:);
            data.driftY = driftXY(2,:);
            % reference pillars flags 
            reference_flags = fread(fileID,file_header.npillars,'int8',readmode);
            data.reference = (reference_flags>0);
        end
    else
        file_header.suc_load = false;
        disp('loading file failed.');
    end

    %% colse the file handle
    fclose(fileID);
end