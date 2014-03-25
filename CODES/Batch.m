%% Create HDF5 files without large datasets contained within

pub_dir = 'https://dl.dropboxusercontent.com/u/22455492/Polymer/';
top_dir = '~/Dropbox/Public/Polymer/';
directories = strsplit(genpath(top_dir),':');

dictname = 'polyer-md-sim';

for ii = 1 : numel( directories )
    if ~strcmp( directories{ii},  top_dir );
        data_files = dir( directories{ii} );
        for ff = 1 : numel( data_files )
            if ~data_files(ff).isdir
                disp( sprintf('Converting file : %s with %i bytes\n',data_files(ff).name,data_files(ff).bytes));
                fn = fullfile(directories{ii}, data_files(ff).name );
                DATA = StructureMD( fn );
                h5out = fullfile( 'DATA', horzcat( data_files(ff).name, '.h5' ) );
                for kk = unique([ 1 : 10 : numel( DATA), numel( DATA) ])
                    % Plot some of the data
                    DATA{kk}.image = VisualizePolymer( DATA{kk}, false);
                    DATA{kk}.link =  sprintf( '%s/%s', regexprep( directories{ii}, top_dir, pub_dir ), data_files(ff).name );
                    
                    for qq = 1 : numel( DATA{kk}.image )
                        s = flickrGetImage( 'flickr.photos.search','Medium 640','text', regexprep(fliplr(strtok(fliplr(DATA{kk}.image{qq}),'/')),'.png',''));
                        DATA{kk}.image{qq} = char(s);
                    end
                    % Add data and image locations before saving the structured
                    % data.
                    
                end
                
                
                
                StructureData( h5out, DATA, false );
                
                ml = createDataset( h5out );
                ml.dict = dictname;
                AttachDictionary( ml, [], dictname );
                PublishDataset(ml);
                
                
            end
        end
    end
end

return

%%

