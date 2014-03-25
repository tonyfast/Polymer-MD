dist_data_loc = '../Data/Distance/';
conv_data_loc = '../Data/Converted/';

ff = dir( horzcat( dist_data_loc,'*.mat') );

ll = [Inf*ones(3,1), -Inf*ones(3,1)];

for ii = 1 : numel(ff)
    
    load( horzcat( conv_data_loc, ff(ii).name), 'lims', 'N' );
    load( horzcat( dist_data_loc, ff(ii).name), 'xx' );
    ll = [ min( lims(:,1), ll(:,1)) max( lims(:,2), ll(:,2)) ];
end
%%

dx = [ 1 1 1 ]*2
xs = [ 0 : dx(1) : max(ll(1,:)) ]; xs = horzcat( -1 * fliplr( xs( 2 : end ) ),xs );
ys = [ 0 : dx(2) : max(ll(2,:)) ]; ys = horzcat( -1 * fliplr( ys( 2 : end ) ),ys );
zs = [ 0 : dx(3) : max(ll(3,:)) ]; zs = horzcat( -1 * fliplr( zs( 2 : end ) ),zs );

ss = [ numel(xs), numel(ys), numel(zs) ];
ms = [max(xs);max(ys);max(zs)];
%%
F = zeros( ss );
for ii = numel(ff)
    
    load( horzcat( conv_data_loc, ff(ii).name), 'A','lims', 'N' );
    load( horzcat( dist_data_loc, ff(ii).name), 'dX', 'xx' );
    
    
    dX(:,:,1) = round2( dX(:,:,1), dx(1))./dx(1);
    dX(:,:,2) = round2( dX(:,:,2), dx(2))./dx(2);
    dX(:,:,3) = round2( dX(:,:,3), dx(3))./dx(3);
    
    B = find(all( bsxfun( @lt, abs(dX), permute( fix( ss./2 ), [ 1 3 2])),3));
    
    
    I = sub2ind( ss, dX(1,B,1) + ms(1)./dx(1)+1, dX(1,B,2) + ms(2)./dx(2)+1, dX(1,B,3) + ms(3)./dx(3)+1 );    
    F(:) = accumarray( I', ones( size(I') ), [  prod( ss ) 1], @sum, 0 );
%     
    I = sub2ind( ss, -1*dX(1,B,1) + ms(1)./dx(1)+1, -1*dX(1,B,2) + ms(2)./dx(2)+1, -1*dX(1,B,3) + ms(3)./dx(3)+1 );    
    F(:) = F(:) + accumarray( I', ones( size(I') ), [  prod( ss ) 1], @sum, 0 );
    return

end
    
    