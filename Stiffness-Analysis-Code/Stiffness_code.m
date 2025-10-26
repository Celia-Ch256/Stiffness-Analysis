% Stiffness_code.m----Modified version of Chacon et al. 2014 JCB code(Stiffness_code_JCB.m)
%
% Purpose: Calculate inter-kinetochore distance and fit the spring constant based on the separation data 
% and saves that info to file in order:inter-kinetochore distance, spring constant
%
% Requires: image processing toolbox, statistics toolbox.
%
% use on low signal : noise images at your own risk!
%
% with our setup, the background is well detected with a gaussian
% filter of stdev = 5. your mileage may vary depending on imaging
% conditions.
%
% this script assumes that each time series is saved as its own tif in numbered order,
% each TIFF contains one color channel, and the suffix of each time series is '1.tif' where
% 1 is the time series number (e.g. 1, 2, ... 999)
%
% Fill out all the relevant info under user info!
%
% Note: At the moment, this code shows you the fits over scatter plots of the image
% data. If the fits do not look right, THEY PROBABLY AREN'T and your image
% is probably too poor of SNR to get useful data.
%


%(1) clear data
clc;
clear;

%(2) user info
pix_convert = 9.498;  % nm/pixel  
total_frames = 224; % n frames in movie  
time_step = 0.022;% sampling interval(s)
pole_color = 1; % color channel
file_prefix = 'data '; % file_prefix
output_file_name = strcat(file_prefix,'Stiffness.xlsx');

x = 0.0;  % save the coordinates of manually selected image centers
y = 0.0;  % save the coordinates of manually selected image centers

%(3) generate a mask based on the first image
fullpath = mfilename('fullpath');    
[path,name]=fileparts(fullpath);  % retrieve runtime path
my_dir = path;
filename = [my_dir,'\',file_prefix,num2str((1)),'.tif'];

I=single(imread(filename,1));% read an image
I = mat2gray(I);

figure(20);
imshow(I);   
set(gcf,'WindowState','maximized');  

[c,r]=ginput(12); % define ROI by selecting a 12-point polygon

Bw= roipoly(I,c,r); % generate a binary mask based on selected region

imwrite(Bw,'mask1.jpg');% export the mask

%(4) generate seed coordinates for centromere centers based on the first image

%(4-1) filter and subtract background
Ss = fspecial('gaussian', 64,7);
filter_Ss = imfilter(I, Ss, 'symmetric');

Ls = fspecial('gaussian', 256, 256);
filter_Ls = imfilter(I, Ls, 'symmetric');
    
I = filter_Ss - filter_Ls;
%(4-2) contrast stretching
imshow(imadjust(1*I),[],'InitialMagnification',400); text(5,5,'click on the center of each pole','Color','white');
set(gcf,'WindowState','maximized');

%(4-3) seed point generation
% manually select the center point and press Enter to confirm; 
% this will be used as a seed point
[x, y, ~] = impixel(); close;

% inter-kinetochore distance calculation(i.e. the variable "dis_mearsure")
xw_t = [x]';
yw_t = [y]';
dis_measure=(sqrt((xw_t(1)-xw_t(2))^2 + (yw_t(1)-yw_t(2))^2))*pix_convert;

% For save distance and coordination list
dis = []; 
coopx = [];
coopy = [];
coopx1 = [];
coopy1 = [];

%(5) extract centromere center point distances from 100 images over 2 seconds

for i0 = 1:100 %total_frames  

    this_frame=single(imread(filename,i0)); % read an image
    this_frame = (this_frame - min(this_frame(:))) / (max(this_frame(:)) - min(this_frame(:)));
    this_frame=100*this_frame;
    
    %(5-2) filter and subtract background
    Ss = fspecial('gaussian', 64,5);
    filter_Ss = imfilter(this_frame, Ss, 'symmetric');
    
    Ls = fspecial('gaussian', 256, 256);
    filter_Ls = imfilter(this_frame, Ls, 'symmetric');
    
    this_frame = filter_Ss - filter_Ls;
    this_frame=immultiply(this_frame,Bw);

    figure(20);
    imshow(0.05*(this_frame));
    title(num2str((i0-1))); 
     
    %(5-3) estimate centromere center positions using a 2D Gaussian fit
    [x1, y1, ~] = fit1Gaussian(this_frame(:,:,pole_color),x(1),y(1),6);
    [x2, y2, ~] = fit1Gaussian(this_frame(:,:,pole_color),x(2),y(2),6);
    
    x_t = [x1 x2]';
    y_t = [y1 y2]';
  
    x_t=x_t*pix_convert; y_t=y_t*pix_convert;% per pixel size 

    %(5-4) distance calculation
    distance = (sqrt((x_t(1)-x_t(2))^2 + (y_t(1)-y_t(2))^2));
    
    dis(i0) = distance;    

    coopx(i0) = x1;
    coopx1(i0) = x2;
    
    coopy(i0) = y1;
    coopy1(i0) = y2;    
    
end

%(6) calculate the residuals of centromere coordinates

%(6-1) calculate coordinate shifts of centromere centers
dis_t= diff(dis);
dis_t = dis_t';

%(6-2)  generate the time sequence

l2=length(dis_t);  % calculate the number of time intervals
X=[];
X(1)=0;
for i=2:l2
    X(i,1)=X(i-1,1)+time_step;% 22ms time intervals
end

%(6-3) build regression model
% solve: dis_t = ¦Á+¦Â*X
% where ¦Á is the slope and ¦Â is the intercept
% b = [¦Á,¦Â]'
b = regress(dis_t, [ones(length(X),1) X]); %Calculate drift corrected residuals

% calculate Yhat
Yhat = [ones(length(X),1) X]*b;
% calculate residuals
resid2 = dis_t - Yhat;
% normalize the residual values
residual=resid2./sqrt(2);

%(7) calculate MSD
%(7-1) for save  data
m = size(residual,1);  % number of residuals
tau_array=[];          % time sequence
msd_array=[];          % MSD sequence
std_msd_array=[];      % MSD variance sequence

%(7-2) core calculation formula
for time_sc = 0:m-1  % How many time points for rough MSD vs time plot
    newarray = [];
    test_end = 1;
    counter = 0;

    while (test_end+time_sc) <= m
        counter = counter + 1;               
        newarray(counter) = sum(residual(test_end:test_end+time_sc));
        test_end = test_end + 1;
    end    
    
    tau_array(time_sc+1) = (time_sc+1)*time_step;   % time sequence
    msd_array(time_sc+1) = (mean(newarray.^2));     % MSD sequence
    std_msd_array(time_sc+1)=(std(newarray.^2));    % MSD variance sequence
end;


%(8) calculate the spring constant from MSD fitted based on confined diffusion motion

%(8-1) define the fitting data used in the regression model
F_Fit_tau_array = tau_array(:); % time sequence
F_Fit_msd_array = msd_array(:); % MSD sequence
F_xdata = F_Fit_tau_array;
F_ydata = F_Fit_msd_array;
%(8-2) define the MSD fitting equation: MSD(t) = ¦Ò^2*(1-exp(-t/tau)) +4Dt + v^2t +b¦È
F_fun = @(F_para,F_xdata)(F_para(1)*(1-exp(-F_xdata/F_para(2)))+4*F_para(3)*F_xdata+F_para(4)^2.*F_xdata.^2+F_para(5));

%(8-3) set initial values for model parameters in the following order: ¦Ò^2, tau,D,v,b¦È
% Note: adjust ¦Ò^2 initial value based on data to improve fitting accuracy
F_para0=[mean(msd_array(1,5:90))-14,0.1,0.01,0.1,mean(msd_array(1,1:5))];

%(8-4) set constraints and perform parameter fitting
options = optimset('MaxFunEvals',2000,'MaxIter',10000);
F_lb=[0,0,0,0,0];                     % lower bound
F_ub=[10000,10000,10000,10000,10000]; % upper bound
F_para = lsqcurvefit(F_fun,F_para0,F_xdata,F_ydata,F_ub,F_lb,options);

%(8-5) calculate RMSE
F_ydataPred = F_fun(F_para,F_xdata);
F_rmse = sqrt(mean((F_ydataPred-F_ydata).^2));

%(8-6) calculate R2
F_r2 = 1 - (sum((F_ydataPred - F_ydata).^2) / sum((F_ydata - mean(F_ydata)).^2));

%(8-7) calculate spring constant(i.e. the variable "k2")
KbT=1000000000000*1000000000*1.38E-23*(273+37);
k2=(KbT/(F_para(1)))*1000 ;  

%(8-8) visual analysis
F_times = linspace(F_xdata(1),F_xdata(end));
figure(3),plot(F_xdata,F_ydata,'ko',F_times,F_fun(F_para,F_times),'b-');
text(0,0,"¦Ò2=" + num2str(F_para(1))+";"+"RMSE="+ num2str(F_rmse)+";"+ "R2="+num2str(F_r2)+";"+"Intercept="+num2str(F_para(5))+";"+"spring constant="+num2str(k2)); 
legend('Data','Fitted exponential');
title('Data and Fitted Curve');

%(8-9) data output
output_file = [my_dir,'\',output_file_name];

if F_para(5) < 50
    if  F_r2 <1
        xlswrite(output_file, {'inter-kinetochore distance'},1,'C1');
        xlswrite(output_file, mean(dis_measure)/1000,1,'C2');
        xlswrite(output_file, {'spring constant'},1,'D1');
        xlswrite(output_file, k2,1,'D2');
    end
end


%% Fit one gaussian in 2D
function [x ,y, fit] = fit1Gaussian(image_to_fit,x_coord,y_coord, offset)
x_coord = round(x_coord);
y_coord = round(y_coord);

image_size = size(image_to_fit);

% while loop makes sure the offsets don't go beyond the image
good_size = 0;
while (good_size == 0)
    if x_coord - offset < 1 || x_coord + offset > image_size(2) || y_coord - offset < 1 || y_coord + offset > image_size(1)
        offset = offset - 1;
    else
        good_size = 1;
        offset;
    end
end

% put the image data into a format gmdistribution can use
xy_data = []; 
counter = 1;
for image_column = x_coord - offset : x_coord + offset
    for image_row = y_coord - offset : y_coord + offset
        pixel_intensity = image_to_fit(image_row,image_column);
        
        for this_intensity = 1 : pixel_intensity
            xy_data(counter,1) = image_column;
            xy_data(counter,2) = image_row;
            counter = counter + 1;
        end
    end
end

fit = gmdistribution.fit(xy_data,1);
x = fit.mu(1);
y = fit.mu(2);

end


