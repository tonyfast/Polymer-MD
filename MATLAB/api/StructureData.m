function StructureData( filename, data, writeall )
% Convert Matlab structure/record array to an HDF5 file.
% filename - an .h5 file name
% data - a structure array containing aggregate and spatial data.
% writeall -  if writeall is true then all the spatial data will be
% written.  false (default)
if ~exist( 'writeall','var')
    writeall = false;
end
    
    for ii = 1 : numel( data )

        if isfield( data{ii}, 'name')
            dsetname = horzcat('/', data{ii}.name);
        else
            dsetname = horzcat( '/DATASET_',num2str( round( 10000*rand(1))) );
        end
        h5CreateGroup( filename, dsetname );
        
        if isfield( data{ii}, 'aggregate')
            h5writecompound( filename, ...
                sprintf( '%s/%s',dsetname,'aggregate'), ...
                data{ii}.aggregate );
        end
        
        if isfield( data{ii}, 'spatial')
            if writeall
                h5writecompound( filename, ...
                    sprintf( '%s/%s',dsetname,'spatial'), ...
                    data{ii}.spatial );
            else
                h5writecompound( filename, ...
                    sprintf( '%s/%s',dsetname,'spatial'), ...
                    structfun( @(x)0, data{ii}.spatial,'UniformOutput',false) );
            end
        end
        
        if isfield( data{ii}, 'link')
            if ischar( data{ii}.link )
                data{ii}.link = cellstr(data{ii}.link);
            end
            h5writeatt( filename,  sprintf( '%s/%s',dsetname,'spatial'),'link', strjoin(data{ii}.link,';') );
        end
        
        
        if isfield( data{ii}, 'image')
            if ischar( data{ii}.image )
                data{ii}.image = cellstr(data{ii}.image);
            end
            h5writeatt( filename,  sprintf( '%s/%s',dsetname,'spatial'),'image', strjoin(data{ii}.image,';') );
        end
        
        
        
        
    end
    
end