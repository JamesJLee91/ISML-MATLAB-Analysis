%% Author: James J Lee

%Assign Savefile names
testname = 'Cylinder_2_Instron_NEW'; 
file = 'Test3.steps.tracking.csv'; %Force file name
file2 = 'Cylinder_2_Strain.csv'; %Strain file name
n = 2; %Test number
decimate = 'Y';

%File directory
directory = 'F:\My Drive\Oxford 2022-24\Caltech Chainmail\Instron\Material\Raw Data';
rawdata = fullfile(directory,file);
rawdata2 = fullfile(directory,file2);
data = readmatrix(rawdata); %load force data
data2 = readmatrix(rawdata2); %load strain data
testdata = readmatrix('F:\My Drive\Oxford 2022-24\Caltech Chainmail\Instron\Material\Pre_Test_Analysis_Characterisation_Cylinder');

% Raw Data and Parameters

testnumbers = testdata(:,1);
m = find(testnumbers==n); %Test identifier
r_S = (testdata(m,5)/2)*1e-3; %radius
A_S = r_S^2 * pi; %area
t_S = testdata(m,11)*1e-3; %thickness

%Raw Data

%Force Time
time = data(:,2);
dt = time(2) - time(1);
%Strain Time
time2 = data2(:,1);

%Force Data
force = -data(:,9)*1e3; %in N (initially kN)

%Strain Data
straindata = -data2(:,2)*1e-2; % percentage to decimal




%ONLY when decimated in image processing
if decimate == 'Y'
    strain = interp(straindata,10); %increase sampling rate
    time_s = interp(time2,10);
elseif decimate == 'N'
    strain = straindata;
    time_s = time2;
end

%%
%Isolate Strain Data
%ZOOM in and press 'Enter' to select the point
% figure;
% title('Choose strain end time')
% plot(time_s,strain)
% zoom on;
% waitfor(gcf,'CurrentCharacter',char(13))
% zoom reset
% zoom off
% grid minor
% % xlim([1e-3,2e-3])
% [x2]=ginput(1);
% %Choose where strain ends
% j=1;
% while time_s(j) < x2(1)       %find the equivalent index
%     j = j+1;
% end
% close gcf %automatically close figure once selection is done
% strainEndTime = time_s(j);
strainEndTime = time_s(end);

%Isolate Force Data
figure;
title('Choose force end time')
zoom on
plot(time,force)
zoom on;
waitfor(gcf,'CurrentCharacter',char(13))
zoom reset
zoom off
grid minor
% xlim([1e-3,2e-3])
[x3]=ginput(1);
%Choose where strain begins and ends
k=1;
while time(k) < x3(1)       %find the equivalent index
    k = k+1;
end
close gcf %automatically close figure once selection is done
forceEndTime = time(k);


%%

%Find the time difference between two end times
timeDiff = strainEndTime - forceEndTime;
timeDiffIndex = ceil(abs(timeDiff)/dt); %convert to index

%If strain time is later, shift strain, if force time is later, shift
%force. Choose the time that isn't shifted.
if timeDiff > 0
    strain_cut = strain(timeDiffIndex:end); %strain is shifted
    time_cut = time_s(1:end-timeDiffIndex+1); %
    time_N = time(1:k); %use force time
    if length(strain) == length(force)
        interp_strain = strain;
    else 
        interp_strain = interp1(time_cut,strain_cut,time_N); %match the data lengths
    end
    strain_N = interp_strain;
    force_N = force(k-length(strain_N)+1:k); %match the beginning of force to strain
  
elseif timeDiff < 0 
    force_cut = force(1:k);
    time_cut = time(1:k);
    force_N = force_cut(timeDiffIndex:end); %force is shifted 'backwards'
    time_N = time(1:k-timeDiffIndex+1); %time is readjusted to the shifted force
    if length(strain) == length(force_cut) %check if the data lengths are same
        interp_strain = strain;
    else
        interp_strain = interp1(time_s,strain,time_N); %if NOT match the data lengths
    end
    strain_N = interp_strain(end-length(force_N)+1:end); %match the beginning of strain to force
      
end

%Engineering Stress on sample
stress_N = force_N/A_S; 
%Strain Rate 
srate=[]; 
for i=1:length(strain_N)-1
    srate(i) = (strain_N(i+1)-strain_N(i))/(time_N(i+1)-time_N(i));
end
%True parameters (assume volume conservation)
t_stress = stress_N.*(1+strain_N);
t_strain = log(1+strain_N);

%% Plotting
figure(3)
subplot(2,2,1)
hold on
grid on
box on
yyaxis left;
plot(time,force,'LineWidth',1.3)
ylabel('Force (N)')
yyaxis right;
plot(time_s,strain,'LineWidth',1.3)
ylabel('Strain')
ylim([0,Inf]);

xlabel('Time (s)')
title('Raw data time sync')
legend('force data','strain data')

subplot(2,2,2)
box on
plot(time_N,stress_N*1e-6,'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Stress (MPa)')
title('Stress History')
grid on

subplot(2,2,3)
box on
plot(time_N,strain_N,'LineWidth',1.3)
ylim([0,Inf])
xlabel('Time (s)')
ylabel('Strain')
title('Strain History')
grid on

subplot(2,2,4)
box on
plot(time_N(1:end-1),srate,'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Strain Rate')
% ylim([0,15]);
grid on 

%Engineering stress-strain=================================================
figure(4)
box on
plot(strain_N,stress_N*1e-6,'LineWidth',1.3);
xlim([0,Inf]);
ylim([0,Inf]);
xlabel('Engineering Strain')
ylabel('Engineering Stress (MPa)')  
title('Engineering Stress-Strain')
grid on

%True stress-strain========================================================
figure(5)
box on
plot(t_strain,t_stress*1e-6,'LineWidth',1.3);
xlim([0,Inf]);
ylim([0,Inf]);
xlabel('True Strain')
ylabel('True Stress (MPa)')  
title('True Stress-Strain')
grid on

%% SAVE FILES============================================================== 

%Save images
destdirectory = 'F:\My Drive\Oxford 2022-24\Caltech Chainmail\Instron\Material\Analysed\';
destdirectory2 ='F:\My Drive\Oxford 2022-24\Caltech Chainmail\Instron\Material\Figures\';
savedir = fullfile(destdirectory,testname);
mkdir(savedir); 
saveas(figure(3), fullfile(savedir, 'parameters.png'));
saveas(figure(4), fullfile(savedir, 'eng_stress_strain.png'));
saveas(figure(5), fullfile(savedir, 'true_stress_strain.png'));

%For csv file
filename = fullfile(savedir, [testname,'.xlsx']);
%Headers
col_header = {'Force_Time_Raw (s)','Force_Raw (kN)','Strain_Time_Raw (s)','Strain_Raw (%)','Time (s)','Stress (Pa)','Strain','True_Stress (Pa)','True_Strain'};
xlswrite(filename,col_header,'Sheet1','A1');

%Raw Data
forcedata = [];
forcedata(:,1) = time;
forcedata(:,2) = force;
straindata = [];
straindata(:,1) = time_s;
straindata(:,2) = strain;
xlswrite(filename,forcedata,'Sheet1','A2');
xlswrite(filename,straindata,'Sheet1','C2');


%Analysed Data
stress_strain = [];
stress_strain(:,1) = time_N;
stress_strain(:,2) = stress_N;
stress_strain(:,3) = strain_N;
stress_strain(:,4) = t_stress;
stress_strain(:,5) = t_strain;
xlswrite(filename,stress_strain,'Sheet1','E2');


fprintf('SAVED!')
