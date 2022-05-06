files = dir('*.bin');
fileID = fopen('samplerates.txt', 'a');
for file = files'
    
    if contains(file.name, 'Wfm') ~= 1;

        [I,t,s] = RTOReadBin(file.name);
        dim = size(I);
        for i=1:dim(2)

            srate = num2sip(1/s.SignalResolution);
            srate = append(regexprep(srate, '\s+', ''),'Sa');
            
            fname = append(file.name,'__',srate,' ',num2str(i),'.csv');
            disp(fname)
            writematrix(I(:,i),fname); 
            
            fprintf(fileID,fname);
            fprintf(fileID,'\n');
            fprintf(fileID,'s.SignalResolution = %d\n', s.SignalResolution);
            fprintf(fileID,'s.TimeScale = %d\n----------\n', s.TimeScale);
            
        end
    end
end
fclose(fileID);
