%workfolder = input('Enter folder with .bin files: \n','s')
addpath(cd)
workfolder = uigetdir
cd(workfolder)

files = dir('*.bin');
fileID = fopen('samplerates.txt', 'a');
for file = files'
    
    if contains(file.name, 'Wfm') ~= 1;

        [I,t,s] = RTOReadBin(file.name);
        dim = size(I);
        for i=1:dim(2)

            srate = num2sip(1/s.SignalResolution);
            srate = append(regexprep(srate, '\s+', ''),'Sa');
            
            disp(append(file.name,'_{wf',num2str(i),'}'))
            fname  = append(file.name(1:end-4),'_{wf',num2str(i),'}_[',srate,']_{Y}.csv');
            writematrix(I(:,i),fname); 
            fnameX = append(file.name(1:end-4),'_{wf',num2str(i),'}_[',srate,']_{X}.csv');
            writematrix(t',fnameX); 
            
            fprintf(fileID,fname);
            fprintf(fileID,'\n');
            fprintf(fileID,'s.SignalResolution = %d\n', s.SignalResolution);
            fprintf(fileID,'s.TimeScale = %d\n----------\n', s.TimeScale);
            
        end
    end
end
fclose(fileID);
