raw_data_loc = '../Data/Raw/pe400tilt6l/';
conv_data_loc = '../Data/Converted/';
ff = dir( raw_data_loc );

 
for ii = find( ~[ff.isdir] );
    fo = fopen( horzcat( raw_data_loc, ff(ii).name) );
 s2 = ['s'];   
   
    
    ct = 0;
    while numel( s2 ) ~= 6
        ct = ct + 1;
        s = fgetl( fo );
        s2 = sscanf( s, '%f' );
        if ct == 2
            T = s2;
        end
        if ct == 4
            N = s2;  % the number of atoms
        end
        
        if ismember( ct, 6 : 8)
            lims( ct - 5,:) = s2;  % size of the simulation
        end
        
        if ct == 10
            break
        end
    end
    
    A = horzcat( s2, fscanf( fo, '%f', [6 N-1]) )';
    
    size(A)
    fclose( fo )
    
    save( horzcat( conv_data_loc, ff(ii).name, '.mat'), 'T','N','lims','A')
end

return

%%

for qq = 5
%     subplot(1,1,qq)
b = ismember( A(:,2), qq);

A2 = A(b,:);

% subplot(3,3,[1 2 4 5 7 8]);
scatter3( A2(:,4), A2(:,6), A2(:,5), 50, A2(:,2), 'filled', 'MarkerEdgeColor', 'k'); axis equal
% subplot(3,3,[3:3:9]); 
% hist( A(b,2), unique( A(b,2) ) );
% plotbrowser

colormap( rand( 100, 3 ));
%%

C = sortrows([A2(:,1), A2(:, [ 4 6 5 3])],1);

if numel( qq ) == 1
P{qq} = C(:,2:4);
end
C(:,1) = C(:,1) - min(C(:,1)) + 1;

C1 = C(:,1);
% S = sparse( C1, circshift( C1, [ 1 0]), 1, max(C1), max(C1));
% S(:) = S + sparse( C1, circshift( C1, [ -1 0]), 1, max(C1), max(C1));
% S(:) = S + eye( max(C1));

S = diag( ones(1, sum(b)) ) + diag( ones(1, sum(b)-1), 1 );

D = Dmatrix( C( : , 2 : 4 ));

if 1
gplot3( S & D < 20, C(:,2:4), 'g:o');
hp = findobj( gcf, 'Color','g','LineStyle',':');
set( hp, 'Color','k','Linestyle','-','MarkerFaceColor','b','Markersize',8)
set( [hp], 'Linewidth',3);



hold on
gplot3( S & D > 20, C(:,2:4), 'c:');
hnp = findobj( gcf, 'Color','c','LineStyle',':');
set( hnp, 'Color','r','Linestyle','--')
% set( hp, 'Color','k','Linestyle','-','MarkerFaceColor','b')
set( [hnp hp], 'Linewidth',2);


bstart = (C(:,end) == 2 );
he = plot3( C(bstart, 2), C(bstart, 3), C(bstart, 4), 'pc','Markersize',30,'Markerfacecolor','m')
hold off

legend( [hp hnp he], 'Internal Bonds','Periodic Bonds','Terminal of Chain')
end

grid on
box on
axis equal


end


save('somedescriptiveariable.mat','P')

%% %

% scatter3( A(:,4), A(:,5), A(:,6), 20, A(:,2), 'filled')
return
%%
C = 1:20;
b = ismember( A(:,2), [ C]);
D = Dmatrix( A( b, 4 : 6 ), [],lims );
[S2,I] = sort( D ,2);
S = S2(:, 2:3);
% ez3dplot( S);

d = pdist( S );
z = linkage( d );
c = cluster( z, 'cutoff',mean(S(:)),'criterion','distance');

subplot(1,2,1);
scatter( S(:,1), S(:,2), 50, c, 'filled', 'MarkerEdgeColor', 'k'); 
subplot(1,2,2);
[y,x] = hist( c , unique(c) );
plot( x, y,'o-');

st = find( c == x( y <= 2), 1, 'first');

nc = sum( b );
p = zeros( nc, 1);
d = zeros( 1, nc);

p(1) = st;
for ii = 2 : sum( b )
    i = find( ~ismember( I( p( ii-1),:), p ), 1, 'first');
    p(ii) = I( p( ii-1),i);
end

[j] = find( b );
p2 = j( p );

S = sparse( repmat( p, 2, 1), ...
    [ circshift( p, [ 1 0 ]); circshift( p, [ -1 0 ])], ...
    1, nc, nc ) + eye(nc);

clf
G = S;
S( S & D > 10) = 0;
gplot3( G, A(b,4:6), 'r:');
hp = findobj( gcf, 'Color','r','LineStyle',':');
hold on
G = S;
S( S & D < 10) = 0;
gplot3( G, A(b,4:6), 'k-o');
hnp = findobj( gcf, 'Color','k','LineStyle','-','Marker','o');

legend( [hp,hnp],'Periodic Connections','NonPeriodic Connections');
grid on
box on
hold off

set( [hp, hnp], 'Linewidth',2);
set( hnp, 'Markerfacecolor','b')