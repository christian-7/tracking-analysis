
base='A647_EGF_10ms_1500mW_COT_Au__'; % select the base name of the folder

for i=[1,2,5,6,7,9,10,12,14,15];

path=['Z:\Christian-Sieben\data_HTP\2016-01-07_blinking_test_A647_EGF_complex\locResults_blinking\' base, num2str(i)];

filename_peaks=[base, num2str(i), '_MMStack_locResults_DC_Merged'];
 
cd(path)

[res]=CG_tracking(filename_peaks); % X, Y, T, ID

cd('Z:\Christian-Sieben\data_HTP\2016-01-07_blinking_test_A647_EGF_complex\locResults_blinking\Tracking')

if  i==1;
    
    filename=['All_tracks_merged' base '.txt'];
    dlmwrite(filename, res);
    
else
    
    filename=['All_tracks_merged' base '.txt'];
    
    preRes=dlmread(filename);
    res(:,4)=res(:,4)+max(preRes(:,4));
    
    new_plus_old=cat(1,preRes,res);
    dlmwrite(filename, new_plus_old);

end

end

fprintf('\n -- Finished Processing --\n')

CG_blink_analysis