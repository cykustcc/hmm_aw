function gen_data_n(dim, len, data, datafile)
% The file C,a,b for mass transportation problem are generated using matlab storing in a binary fi#
    fid = fopen(datafile,'w');
    fwrite(fid,[dim,len],'int');
    fwrite(fid,data,'float');
    fclose(fid);
end