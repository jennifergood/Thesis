clc
clear all
close all

%Move files to Matlab drive/JennFolder data for Symplectic Check and STM
deg = 10;

%Create files for Symplectic Checks and STMs from Fortran
% copyfile fort.41 Symp1.txt
% copyfile fort.42 Symp2.txt
% copyfile fort.43 Symp3.txt
% copyfile fort.44 Symp4.txt
% 
% copyfile fort.53 STM2.txt
% copyfile fort.54 STM3.txt
% copyfile fort.59 STM4.txt

name = dir('STM*.txt');
for i = 2:length(name)
    Source = sprintf('Degrees.%d/STM*',deg);
    Destination = 'C:\Users\jgood3383\MATLAB Drive\Thesis\Jenn\Data\STMDATA';
    copyfile(Source,Destination)
end
