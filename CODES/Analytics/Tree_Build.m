% Build and save the kd-trees from the converted data

conv_data_loc = '../Data/Converted/';
tree_data_loc = '../Data/Tree/';

ff = dir( horzcat( conv_data_loc,'*.mat') );

for ii = 30 : numel( ff )
    %     load( horzcat( conv_data_loc, ff(ii).name),  );
    load( horzcat( conv_data_loc, ff(ii).name), 'A', 'lims');
    
    xx = diff( lims, [], 2 );
    % Minimum image convention since the data is periodic
    An = Periodic_Data_For_Tree( A(:,4:6), lims );
    
    tr2 = kdtree_build( An );
    % Store tr in matlab memory for saving
    tr = kdtree_io_to_mat( tr2 );
    kdtree_delete( tr2 ); % delete tree
    pause(1.5)  % This gives the computer time to clean things up in RAM
    save( horzcat( tree_data_loc, ff(ii).name), 'tr');
    disp(ii)
end