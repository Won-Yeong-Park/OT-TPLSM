%% Simple stack
clear all
close all
clc
fclose all;
for channel= 1:2

clearvars -except channel
close all
clc
fclose all;

% pwd
step=1;
step_ratio=step/0.325;
height_max=850;%raw image상 아
height_min=338;%raw image상 위

scan_direction=2;%1:plus,2:minus
sampling_method=1;%1: undersampling, 2: upsampling, 3: noting
re_skew=1;%%1:yes,2:no

if channel==1
mkdir C1
files = dir('*C1*.tif');
max_=35;%Local normalization value
min_=0;%Local normalization value
cd C1
avg=2;
save('Build_Parameter_C1.mat','max_','min_','step_ratio','height_max','height_min','scan_direction','sampling_method','avg')

elseif channel==2
mkdir C2
files = dir('*C2*.tif');
max_=125;
min_=85;
cd C2
avg=10;
save('Build_Parameter_C2.mat','max_','min_','step_ratio','height_max','height_min','scan_direction','sampling_method','avg')

else %microsphere measurment mode
mkdir out
files = dir('*im*.tif');
scan_direction=1;
sampling_method=3;
avg=1;
cd out
save('Build_Parameter.mat','step_ratio','height_max','height_min','scan_direction','sampling_method','avg')
end



n_files = size(files,1);
image_3d_comp = zeros(height_max-height_min,2060,round(n_files/step));

%% image read along scanning direction
  tic 
if scan_direction==1
 
parfor i_image=1:n_files-1
    temp_image = imread([files(n_files-i_image).folder,'\',files(n_files-i_image).name]);   
    image_3d_comp(:,:,i_image) = temp_image(height_min:height_max-1,:);
end

else
    
parfor i_image=1:round(n_files/step)-1
    temp_image = imread([files(i_image).folder,'\',files(i_image).name]);
    image_3d_comp(:,:,i_image) = temp_image(height_min:height_max-1,:);
end

end
    toc
clear temp_image


%% sampling method
tic
if sampling_method==1
    image_3d_comp_=imresize3(image_3d_comp, [(height_max-height_min)/(step_ratio*sqrt(2)) 2060/step_ratio n_files]);%undersampling   
elseif sampling_method==2
    image_3d_comp_=imresize3(image_3d_comp, [(height_max-height_min)/sqrt(2) 2060 n_files*step_ratio]);%upsampling
else
    image_3d_comp_=image_3d_comp;  
end

clear image_3d_comp %memory save
pause(8) %time for menoey stable
%% matrix re-skew
n_file = size(image_3d_comp_);
if re_skew==1
image_3d_tilt = zeros(n_file(1),n_file(2),n_file(3)+n_file(1)-1);

% re-skew
    for i=1:n_file(3)
%         disp(i)
        for j=1:n_file(1)
             
             image_3d_tilt(end-j+1,:,i+j-1) = image_3d_comp_(end-j+1,:,i);

        end
    end
    
toc

clear image_3d_comp_%memory save


imageA_volume=zeros(size(image_3d_tilt,2),size(image_3d_tilt,3),avg);

else
imageA_volume=zeros(size(image_3d_comp_,2),size(image_3d_comp_,3),avg);
end
    

for i=1:n_file(1)-avg+1

    if channel==1%fluorescence save
    disp(i)
       if re_skew==1
    for i_avg=1:avg
        imageA_volume(:,:,i_avg)=squeeze(image_3d_tilt(i+i_avg-1,:,:));
    end
    
       else
        for i_avg=1:avg
        imageA_volume(:,:,i_avg)=squeeze(image_3d_comp_(i+i_avg-1,:,:));
        end
       end
    imageA_volume_mean=mean(imageA_volume,3);
    
    normalization=real(localnormalize(double(imageA_volume_mean),400,400));
    normalization=normalization-min(normalization(:));
    normalization=normalization-min_;
    normalization=normalization./(max_-min_); 
    normalization=uint16(normalization.*2^16);
    normalization=uint16(imageA_volume_mean);

    normalization_sharp=imsharpen(normalization,'Radius',2,'Amount',1);
    save=ind2gray(normalization_sharp,gray(65536));
    imageoutput=sprintf('%s%d%s','image_',i,'.tif');
    imwrite(save,imageoutput)

    elseif channel==2 %SHG save
    disp(i)        
        if re_skew==1
    for i_avg=1:avg
        imageA_volume(:,:,i_avg)=squeeze(image_3d_tilt(i+i_avg-1,:,:));
    end
    
        else
        for i_avg=1:avg
        imageA_volume(:,:,i_avg)=squeeze(image_3d_comp_(i+i_avg-1,:,:));
        end
       end
        
        imageA_volume_mean=mean(imageA_volume,3);
        
        imageA_volume_mean=(imageA_volume_mean-min_);
        imageA_volume_mean=imageA_volume_mean.*(1/(max_-min_));
        imageA_volume_mean=imageA_volume_mean.*65536;
        SHG=uint16(imageA_volume_mean);
    % SHG save

    SHG_column_offsets = median(SHG,2);
    SHG_column_offsets = SHG_column_offsets - min(SHG_column_offsets);
	SHG_new_im = bsxfun(@minus,SHG,SHG_column_offsets);
    SHG_new_im =imgaussfilt(SHG_new_im,2);
    
    save=ind2gray(SHG_new_im,gray(65536));
    imageoutput=sprintf('%s%d%s','image_',i,'.tif');
    imwrite(squeeze(save),imageoutput)
    
    else
    save=squeeze((image_3d_tilt(i,:,:)));
    save=uint16(save);
    imageoutput=sprintf('%s%d%s','image_',i,'.tif');
    imwrite(squeeze(save),imageoutput)
    end
    
        
end
    
end