%% Low Rate Instron Analysis Non-Uniform Samples 2023 
%% Author: James J Lee

%Assign Savefile names=====================================================
testname = 'SCR_12_Instron_REV'; 
file = 'Test14.steps.tracking.csv'; %Force file name
n = 14; %Test number

%File directory============================================================
directory = 'C:\Users\cieli\Google Drive\Oxford 2022-23\Caltech Chainmail\Instron\Structure\Raw Data';
rawdata = fullfile(directory,file);
data = readmatrix(rawdata); %load force data
testdata = readmatrix('C:\Users\cieli\Google Drive\Oxford 2022-23\Caltech Chainmail\Instron\Structure\Pre_Test_Analysis');

%Raw Data
time = data(:,1);
time = time-time(1); %zero
mdt = diff(time); %check time increments
force = -data(:,9)*1e3; %in N (initially kN)
force = force-force(1); %zero

%Sample Parameters=========================================================
r_S = (testdata(n,5)/2)*1e-3; %radius
A_S = r_S^2 * pi; %area
t_S = testdata(n,5)*1e-3; %thickness

%Strain from the Instron data
disp = -data(:,10)*1e-3; %in m (initially mm)
time2 = data(:,1);
strain =[];
for m = 1:length(disp)
    eps=(disp(m)-disp(1))/t_S;
    strain(m,1)=eps;
end

%isolate force
figure(1)
plot(time,force)
grid minor
% xlim([1e-3,2e-3])
[x3,y3]=ginput(2);
%Beginning of pulse
i=1;
while time(i) < x3(1)       %find the equivalent index
    i = i+1;
end
%End of pulse
j=1;
while time(j) < x3(2)
    j = j+1;
end
close gcf %automatically close figure once selection is done

force2 = force(i:j);

%isolate strain
figure(6)
plot(time2,strain)
grid minor
% xlim([1e-3,2e-3])
[x4,y4]=ginput(1);
%Beginning of pulse
p=1;
while time2(p) < x4(1)       %find the equivalent index
    p = p+1;
end
close gcf %automatically close figure once selection is done

strain2 = strain(p:end);
time2_n = time2(1:length(strain2));

%match the lengths
m = length(force) - length(strain2);
if m>0
    force3 = force(1:length(strain2));
    strain3 = strain2;
    time3 = time2_n;
else
    force3 = force;
    strain3 = strain2(1:length(force));
    time3 = time2_n(1:length(strain3));
end

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

%% Plotting================================================================

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

%% SAVE FILES============================================================== 

%Save images
destdirectory = 'C:\Users\cieli\Google Drive\Oxford 2022-23\Caltech Chainmail\Instron\Structure\Analysed\';
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
