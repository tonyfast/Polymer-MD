return
conv_data_loc = '../Data/Converted/';
dist_data_loc = '../Data/Distance/';

ff = dir( horzcat( conv_data_loc,'*.mat') );

for ii = 1 : numel( ff )
    %     load( horzcat( conv_data_loc, ff(ii).name),  );
    load( horzcat( conv_data_loc, ff(ii).name), 'A', 'lims', 'N' );
    
    xx = diff( lims, [], 2 );
    
    dX = zeros( 1, (N^2 - N)./2,3);
    [ dX(:,:,1) dX(:,:,2) dX(:,:,3) ] = ...
        deal( squareform( bsxfun( @minus, A(:,4) , A(:,4)' ) ), ...
        squareform( bsxfun( @minus, A(:,5) , A(:,5)' ) ), ...
        squareform( bsxfun( @minus, A(:,6) , A(:,6)' ) ) ) ;
    
    % make the distances periodic
    dX(:) = dX + bsxfun( @gt, abs(dX), permute( xx./2, [ 2 3 1])) .*  bsxfun( @times, sign(dX), -1.*permute( xx, [ 2 3 1]));
    
    
    save( horzcat( dist_data_loc, ff(ii).name), 'dX','xx');
    disp(ii)
end