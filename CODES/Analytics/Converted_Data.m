% This function takes the raw text output from a simulation and converts it
% into a useable matlab format with some meta information asssinged to it

raw_dir = '..\Data\Raw\';
conv_dir = '..\Data\Converted\';

dd = dir( raw_dir );
for ii = find( ~[dd.isdir] )
    fo = fopen( horzcat( raw_dir, dd(ii).name ), 'r');
    t_ct = 0;
    while  t_ct < 1250 & ~feof(fo)  % stop after a certain number of times steps
        t_ct = t_ct + 1;
        
        fgetl(fo);
        t = fscanf( fo, '%i\n', 1);
        
        for jj = 1 : 3
            fgetl( fo );
        end
        
        lims = fscanf( fo, '%f %f\n', [ 2 3] )';
        
        fgetl( fo );
        
        A = fscanf(fo, '%f %f %f %f %f %f\n',[ 6 7999 ])';
        if (mod( t_ct,100)-1) == 0
            fn = dd(ii).name;
            % Each file is saved with a unique index that maps back to a
            % timestep in a file, I probably should save file pointers
            save( horzcat( conv_dir, num2str(ii),'_',num2str(t_ct),'.mat'), 'A','fn','lims','t','t_ct');
        end
        if t_ct> 12501; end
        disp( [ii t_ct'])
    end
    
    fclose(fo)
    
end