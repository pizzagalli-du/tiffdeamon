%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts the TIFF files of the 2PM microscope (IRB specific) to AVI     %
% Automatically finds the contrast of the channels                        %
% Automatically adjusts the channels and creates a MIP projection         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

load('start_path');
TIFF_PATH = uigetdir(start_path);
start_path = TIFF_PATH;
save('start_path.mat', 'start_path');

lookup_table = [0,1,0; 1,0,0; 0,0,1; 0.7,0.7,0.7];
lookup_names= {'GFP', 'RFP', 'PB', 'FR'};
lookup_idx = zeros(1,4);

files = dir([TIFF_PATH, '/', '*.tif']);

strT = 'Time Time';
strC = '] _C';
strZ = ' Z';
strC_name_start = 'PMT [';

Z = 0;
C = 0;
T = 0;

% first obtain dataset size
for file = files'
    currfn = file.name;
    parts = strsplit(currfn,'_');
    
    strT_idx = strfind(currfn, strT);
    strC_idx = strfind(currfn, strC);
    strZ_idx = strfind(currfn, strZ);
    strC_name_start_idx = strfind(currfn, strC_name_start);
    
    if(~isempty(strT_idx))
        curr_t = str2num(currfn(strT_idx+(9:12)));
    else
        curr_t = 0;
    end
    %TODO: solo se seguito da 4 numeri
    if(~isempty(strZ_idx))
        curr_z = str2num(currfn((strZ_idx+2):(strZ_idx+5)));
    else
        curr_z = 0;
    end
    
    if(~isempty(strC_idx))
        curr_c = str2num(currfn(strC_idx+(4:5)));
    else
        curr_x = 0;
    end
    
    if(~isempty(strT_idx))
        curr_t = str2num(currfn(strT_idx+(9:12)));
    else
        curr_t = 0;
    end
    
    if(~isempty(strC_name_start_idx))
        curr_c_name = currfn((strC_name_start_idx+5):strC_idx-1);
        lookup_idx(curr_c+1) = find(strcmp(lookup_names, curr_c_name));
    else
        curr_c_name = 'GFP';
    end
    
    Z = max(Z, curr_z);
    C = max(C, curr_c);
    T = max(T, curr_t);
end

Z = Z+1;
C = C+1;
T = T+1;
Itemp = imread([TIFF_PATH,'/',currfn]);
H = size(Itemp,1);
W = size(Itemp,2);

% Allocate memory to store data
zstack = zeros(H,W,Z,C,T,'uint16');

% Read the entire dataset
h = waitbar(0, 'Reading dataset');
totfiles = numel(files');
    
countfiles = 0;
for file = files'
    currfn = file.name;
    parts = strsplit(currfn,'_');
    
    strT_idx = strfind(currfn, strT);
    strC_idx = strfind(currfn, strC);
    strZ_idx = strfind(currfn, strZ);
    strC_name_start_idx = strfind(currfn, strC_name_start);
    
    if(~isempty(strT_idx))
        curr_t = str2num(currfn(strT_idx+(9:12)));
    else
        curr_t = 0;
    end
    
    if(~isempty(strZ_idx))
        curr_z = str2num(currfn((strZ_idx+2):(strZ_idx+5)));
    else
        curr_z = 0;
    end
    
    if(~isempty(strC_idx))
        curr_c = str2num(currfn(strC_idx+(4:5)));
    else
        curr_x = 0;
    end
    
    if(~isempty(strT_idx))
        curr_t = str2num(currfn(strT_idx+(9:12)));
    else
        curr_t = 0;
    end
    
    if(~isempty(strC_name_start_idx))
        curr_c_name = currfn((strC_name_start_idx+5):strC_idx-1);
    else
        curr_c_name = 'GFP';
    end
    
    curr_t = curr_t + 1;
    curr_z = curr_z + 1;
    curr_c = curr_c + 1;
        
    if (curr_t <= T) && (curr_z <= Z) && (curr_c <= C)
        Itemp = imread([TIFF_PATH,'/',currfn]);
        if(size(Itemp) == size(zstack(:,:,curr_z,curr_t)))
            zstack(:,:,curr_z,curr_c,curr_t) = Itemp;
        else
            errordlg(['Error in ', file.name]);
            return;
        end
    else
        disp 'Warning: there are files not included in the size of the datset';
    end
    countfiles = countfiles + 1;
    if(totfiles > 0)
        waitbar(countfiles / totfiles, h);
    end
end

% creating PNG snapshot
% exporting MP4 preview
%TODO: MIP for each channel
waitbar(0, h, 'Writing video');
[folder_parent,folder_name,~] = fileparts(TIFF_PATH);
outputVideo = VideoWriter([folder_parent,'/',folder_name,'_preview.mp4'], 'MPEG-4');
outputVideo.FrameRate = 5;
open(outputVideo)

I_temp_rgb = zeros(H,W,3);
for curr_t = 1:T
    MIP_RGB = zeros(H,W,3);
    for curr_c = 1:C
        u = zstack(:,:,:,curr_c,curr_t);
        u = imadjust(u(:));
        v = reshape(u, [H,W,Z]);
        curr_lookup = lookup_table(lookup_idx(curr_c), :);
        for curr_z = 1:Z
            I_temp = double(v(:,:,curr_z));
            for curr_rgb = 1:3
                I_temp_rgb_ch = I_temp * curr_lookup(curr_rgb);
                I_temp_rgb(:,:,curr_rgb) = I_temp_rgb_ch;
            end
        end
        MIP_RGB = max(MIP_RGB, I_temp_rgb);
    end
    writeVideo(outputVideo, mat2gray(MIP_RGB));
    if(totfiles > 0)
        waitbar(curr_t / T, h);
    end
end
close(outputVideo);
close(h);
implay([folder_parent,'/',folder_name,'_preview.mp4']);