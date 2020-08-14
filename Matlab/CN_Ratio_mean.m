% -------------------------------------------------------------------------
%   This script accompanies the manuscript                                 
%   Simon et al., (2020) Developmental Cell                                
%   Repository available on https://github.com/therealkatlab               
%   Please consult READ_ME for more information                            
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% This script is for use with EmbryoPipeline_*.cpipe CellProfiler pipelines
% Calculates cytoplasmic : nuclear (C:N) ratio of mean ERK-KTR-mClover 
% fluorescence intensity after background subtration.
%
% Inputs: Segmentation of pre-nuclei from CellProfiler as *.tiff images
%         Mean intensity measurements of background, nuclei and cytoplasm 
%           from CellProfiler as *.mat file
%
% Output: 16-bit images with raw C:N values (*10000) for import into Imaris 
%         for further downstream analysis
% -------------------------------------------------------------------------

%% Image structure and inputs for each embryo

t = 24 ; % total time points
z = 30 ; % total z slices
totalimages = t * z ;

path = 'Y:\LAB DATA CONFOCAL\Claire\ERK-KTR Clover\100319_ERK_KTR_hom_mKate2_E4.0\Image_analysis\Cell_profiler\p7_A'
% path with *.mat and *.tiff files output from CellProfiler

%% Matrix for all images

Times = ones(1,totalimages) ;
c = 1 ;
for a = 1:t
    for b = 1:z
        Times(c) = a ;
        c = c+1 ;
    end
end

Zs = ones(1, totalimages) ;
c = 1 ;
for a = 1:t
    for b = 1:z
        Zs(c) = b ;
        c = c+1 ;
    end
end
        
%% Calculate C:N for each z-slice and time-point

for ii = 1
    X{1} = load([path,'\DefaultOUT.mat']) ; % CellProfiler *.mat output

    mkdir([path,'\OutputCNmeanratio']) ; % Output folder directory to write C:N images to

   for a = 1:totalimages
        Img = imread([path,'\img_',sprintf('%09d',a),'_mKate_000_prenuclei.tiff']) ; % CellProfiler *.tiff segmentation
        Bg_ERKKTR(1) = min([X{1}.handles.Measurements.Background.Intensity_MeanIntensity_ERKKTR{1,a}]) ;
        % Background intensity for the image measured in CellProfiler
        FinalImg = zeros(size(Img));
        for b = 1:length(X{1}.handles.Measurements.nuclei.Intensity_MeanIntensity_ERKKTR{1,a}) 
            % for each Z plane get the mean nuclear intensity for each time point
                Ratio = (X{1}.handles.Measurements.cytoring.Intensity_MeanIntensity_ERKKTR{1,a}(b)-Bg_ERKKTR(1))...
                ./(X{1}.handles.Measurements.nuclei.Intensity_MeanIntensity_ERKKTR{1,a}(b)-Bg_ERKKTR(1)) ; 
                % calculates the ratio between the cytoplasm and nuclei for each time point and each Z plane
          
            FinalImg(Img==b) = Ratio*10000; % this should not excede 18-bit range of image
           
       end
       imwrite(uint16(FinalImg),[path,'\OutputCNmeanratio\img_',sprintf('%09d',Times(a)),'_CNRatioERK_',sprintf('%03d',Zs(a)),'.TIF']);
       
   end
    
    
end