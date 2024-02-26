% Executing MRXCAT_Brain for variable tissue parameters, variable sequence
% paramters and SNR levels

clc
clear

%----------------------------
% path and filename for XCAT subjects or subjects in bin and log file
% format
%----------------------------

filepath = 'your path to subject';
subject = 'name of subject';
fileID = fopen([date,'_exp.txt'],'w');
c=1; % no of contrasts

finfo = dir (filepath);
filename = finfo(3).name;

pars = generate(filepath,filename,subject,c);
simpar = struct(subject,pars);

fprintf(' subject was simulated.....')
save('simpar.mat')
fclose(fileID);