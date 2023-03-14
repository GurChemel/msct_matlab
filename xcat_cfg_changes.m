clear;clc;
fid=fopen('230202/original_xcat2.cfg');

data_bin=fread(fid,'uint32');
disp('OLD')
disp(datetime(data_bin(2),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS'))
disp(datetime(data_bin(9),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS'))
disp(datetime(data_bin(end),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS'))

% data_bin(9)=data_bin(end)+60*60*24*365;
data_bin(end)=data_bin(end)+2*60*60*24*365;


fileID = fopen('230202/newest_xcat2.cfg','w');
fwrite(fileID,data_bin,'uint32');
fclose(fileID);

fclose(fid);

fid=fopen('230202/newest_xcat2.cfg');
disp('NEW')
data_bin=fread(fid,'uint32');
disp(datetime(data_bin(2),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS'))
disp(datetime(data_bin(9),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS'))
disp(datetime(data_bin(end),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS'))
fclose(fid);
