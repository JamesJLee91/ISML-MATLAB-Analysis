%% Low Rate Instron Analysis 2023 
%% Author: James J Lee
%Assign Savefile names=====================================================
testname = 'Cylinder_1_Instron'; 
file = 'Test2.steps.tracking.csv'; %Force file name
file2 = 'Cylinder_1_Strain.csv'; %Strain file name
n = 1; %Test number

%File directory============================================================
directory = 'F:\My Drive\Oxford 2022-24\Caltech Chainmail\Instron\Material\Raw Data';
rawdata = fullfile(directory,file);
rawdata2 = fullfile(directory,file2);
data = readmatrix(rawdata); %load force data
data2 = readmatrix(rawdata2); %load strain data
testdata = readmatrix('F:\My Drive\Oxford 2022-24\Caltech Chainmail\Instron\Material\Pre_Test_Analysis_Characterisation_Cylinder');

%%
%Raw Data
time = data(:,1);
time = time-time(1); %zero
mdt = diff(time); %check time increments
force = -data(:,9)*1e3; %in N (initially kN)
force = force-force(1); %zero

%Sample Parameters=========================================================
testnumbers = testdata(:,1);
m = find(testnumbers==n); 
r_S = (testdata(m,5)/2)*1e-3; %radius
A_S = r_S^2 * pi; %area
t_S = testdata(m,11)*1e-3; %thickness

%Strain Data
straindata = -data2(:,2)/1000;
straindata = straindata-straindata(1);
strain = straindata;

%%ONLY when decimated in image processing
strain = interp(straindata,10); %increase sampling rate

%Time for the strain data
dt = 1e-6;
time2 = zeros(length(strain),1);
for i = 1:length(strain)
    time2(i) = dt*(i-1);
end

%Match data lengths

t_lag = time(end)-time2(end); %time lag between force and strain data
m = round(t_lag/mdt(1)); %convert in terms of index
if m>0 %if force data is longer
    force_n = force(abs(m)+1:end); %adjusted force
    strain_n = strain;
    time_n = time2;
else %if strain data is longer
    force_n = force;
    strain_n = strain(abs(m)+1:end); %adjusted strain
    time_n = time2(abs(m)+1:end);
end

%isolate strain data
figure(1)
plot(time_n,strain_n)
grid minor
% xlim([1e-3,2e-3])
[x1,y1]=ginput(1);
%Beginning of incident pulse
i=1;
while time_n(i) < x1(1)       %find the equivalent index
    i = i+1;
end
close gcf %automatically close figure once selection is done

%isolate force data
figure(5)
plot(force_n)
grid minor
[x2,y2]=ginput(1);
%Beginning of incident pulse
j=1;
while time_n(j) < x2(1)       %find the equivalent index
    j = j+1;
end
close gcf %automatically close figure once selection is done

%Isolated
strain2 = strain_n(i:end); 
force2 = force_n(j:end);

%Make both equal lengths
if length(strain2)<length(force2)
    k = length(strain2);
else
    k = length(force2);
end

%New adjusted data
strain3 = strain2(1:k);
force3 = force2(1:k);
time3 = time_n(1:k);

%Parameter Calculations

%Engineering Stress on sample
stress = force3/A_S; 
%Strain Rate 
srate=[]; 
for i=1:length(strain3)-1
    srate(i) = (strain3(i+1)-strain3(i))/(time3(i+1)-time3(i));
end
%True parameters
t_stress = stress.*(1+strain3);
t_strain = log(1+strain3);

%Plotting================================================================

figure(2)
subplot(2,2,1)
hold on
grid on
box on
plot(time3,force3,'LineWidth',1.3)
plot(time3,strain3,'LineWidth',1.3)
xlabel('Time (s)')
title('Raw data time sync')
legend('force data','strain data')

subplot(2,2,2)
box on
plot(time3,stress*1e-6,'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Stress (MPa)')
title('Stress History')
grid on

subplot(2,2,3)
box on
plot(time3,strain3,'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Strain')
title('Strain History')
grid on

subplot(2,2,4)
box on
plot(time3(1:end-1),srate,'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Strain Rate')
% ylim([0,15]);
grid on 

%Engineering stress-strain=================================================
figure(3)
box on
plot(strain3,stress*1e-6,'LineWidth',1.3);
xlabel('Engineering Strain')
ylabel('Engineering Stress (MPa)')  
title('Engineering Stress-Strain')
grid on

%True stress-strain========================================================
figure(4)
box on
plot(t_strain,t_stress*1e-6,'LineWidth',1.3);
xlabel('True Strain')
ylabel('True Stress (MPa)')  
title('True Stress-Strain')
grid on

% SAVE FILES============================================================== 

%Save images
destdirectory = 'C:\Users\cieli\Google Drive\Oxford 2022-23\Caltech Chainmail\Instron\Material\Analysed\';
savedir = fullfile(destdirectory,testname);
mkdir(savedir); 
saveas(figure(2), fullfile(savedir, 'plots.png'));
saveas(figure(3), fullfile(savedir, 'stress_strain.png'));
saveas(figure(4), fullfile(savedir, 'true_stress_strain.png'));

%For csv file
filename = fullfile(savedir, [testname,'.xls']);

%Headers
col_header = {'Force_Time_Raw','Force_Raw','Strain_Time_Raw','Strain_Raw','Time','Stress','Strain','True_Stress','True_Strain'};
xlswrite(filename,col_header,'Sheet1','A1');

%Raw Data
forcedata = [];
forcedata(:,1) = time;
forcedata(:,2) = force;
straindata = [];
straindata(:,1) = time2;
straindata(:,2) = strain;
xlswrite(filename,forcedata,'Sheet1','A2');
xlswrite(filename,straindata,'Sheet1','C2');

%Analysed Data
stress_strain = [];
stress_strain(:,1) = time3;
stress_strain(:,2) = stress;
stress_strain(:,3) = strain3;
stress_strain(:,4) = t_stress;
stress_strain(:,5) = t_strain;
xlswrite(filename,stress_strain,'Sheet1','E2');

