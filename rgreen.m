function [Gss Gds Gcl]=rgreen(f)

%DMM 01/2011
%
%Read Green function files








%STRIKE SLIP

fid=fopen([f '.ss']);
%first 12 lines are headers, 13th line contains info on file structure,
%keep it!
for k=1:12
    fstruc=fgets(fid);
end
fstruc=str2num(fstruc);
%Define structure sizes for reading Green functions
N1=fstruc(1);   %No. of radial distances
N2=fstruc(4);   %No of depth slices
N=N1*N2;
M=10;
Gss=fscanf(fid,'%f%f%f%f%f%f%f%f%f%f',[M N]);
Gss=Gss';
fclose(fid);


%DIP SLIP

fid=fopen([f '.ds']);
%first 12 lines are headers, 13th line contains info on file structure,
%keep it!
for k=1:12
    fstruc=fgets(fid);
end
fstruc=str2num(fstruc);
%Define structure sizes for reading Green functions
Gds=fscanf(fid,'%f%f%f%f%f%f%f%f%f%f',[M N]);
Gds=Gds';
fclose(fid);


%CLVD

fid=fopen([f '.cl']);
%first 12 lines are headers, 13th line contains info on file structure,
%keep it!
for k=1:12
    fstruc=fgets(fid);
end
fstruc=str2num(fstruc);
%Define structure sizes for reading Green functions
M=7;  %Displacements+strains+tilts (3 less for CLVD)
Gcl=fscanf(fid,'%f%f%f%f%f%f%f%f%f%f',[M N]);
Gcl=Gcl';
fclose(fid);