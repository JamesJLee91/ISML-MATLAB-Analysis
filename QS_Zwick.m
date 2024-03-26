%Load data


testname = 'Batch1_NAM4_QS_B2_2';
file = 'Batch1_NAM4_QS_B2_2';
n = 73;

directory = 'F:\My Drive\Oxford 2022-24\NAM 3-4\Test Data';
rawdata = fullfile(directory,file);
data = readmatrix(rawdata); %load scope data
testdata = readmatrix('F:\My Drive\Oxford 2022-24\NAM 3-4\NAM Test Records');

%Save directory
destdirectory = 'F:\My Drive\Oxford 2022-24\NAM 3-4\Analysed\';
savedir = fullfile(destdirectory,testname);
savedir2 = fullfile('F:\My Drive\Oxford 2022-24\NAM 3-4\Figures\NAM4\QS\');
mkdir(savedir); 


%Locate the relevant row 
testnumbers = testdata(:,1);
m = find(testnumbers==n);
%Get data
x = data(:,3)*1e-3;
force = data(:,2);
time = data(:,1);

x2 = (-1.484e-21*force.^4) + (1.092e-16*force.^3) + (-2.994e-12*force.^2) + (6.064e-8*force) + (8.574e-6) ; %compliance
displacement = x - x2; %sample displacement

%Sample Parameters
L_S = mean(testdata(m,3:6))*1e-3; %sample thickness
R_B = (12.3e-3)/2; %compression radius
A_B = pi*(R_B^2); %compression area

% %%
% hold on
% plot(x)
% plot(x2)
% plot(displacement)
% legend('total','rig','sample');


%isolate data
figure(1)
plot(time,force)
%xlim([1e-3,2e-3])
[x1,y1]=ginput(2);
%Beginning of incident pulse
i=1;
while time(i) < x1(1)       %find the equivalent index
    i = i+1;
end
% Manually selecting points when signal is not clear
% Beginning of reflected pulse
j=1;
while time(j) < x1(2)
    j = j+1;
end
close gcf %automatically close figure once selection is done


% Parameter Calculations

f = force(i:j);
x = displacement(i:j);
ntime = time(i:j);

stress = f/A_B;
strain =[];
for q = 1:length(x)
    eps = (x(q)-x(1))/L_S;
    strain(q,1)=eps;
end

%zero 
stress = stress-stress(1);
strain = strain-strain(1);
ntime = ntime-ntime(1);

%Plotting
figure(2)
subplot(2,1,1)
hold on
grid on
box on
plot(ntime,stress*1e-6,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Stress (MPa)')

subplot(2,1,2)
hold on
grid on
box on
plot(ntime-10,strain,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Strain')

figure(3)
hold on
grid on
box on
plot(strain,stress*1e-6,'LineWidth',1.5)
xlabel('Strain')
ylabel('Stress (MPa)')


%% SAVE FILES==============================================================

%Save images
saveas(figure(2), fullfile(savedir, 'plots.png'));
saveas(figure(3), fullfile(savedir, 'stress_strain.png'));
% saveas(figure(4), fullfile(savedir, 'true_stress_strain.png'));

%For csv file
filename = fullfile(savedir, [testname,'.xlsx']);
filename2 = fullfile(savedir2,[testname,'.xlsx']);

%Headers
col_header = {'Time_Raw','Standard_Force','Standard_Travel','Time','Stress','Strain'};
xlswrite(filename,col_header,'Sheet1','A1');
xlswrite(filename2,col_header,'Sheet1','A1');

%Raw Data
rawdata = [];
rawdata(:,1) = time;
rawdata(:,2) = force;
rawdata(:,3) = x2;
xlswrite(filename,rawdata,'Sheet1','A2');
xlswrite(filename2,rawdata,'Sheet1','A2');

%Analysed Data
finaldata = [];
finaldata(:,1) = ntime;
finaldata(:,2) = stress;
finaldata(:,3) = strain;
xlswrite(filename,finaldata,'Sheet1','D2');
xlswrite(filename2,finaldata,'Sheet1','D2')

fprintf('SAVED!')