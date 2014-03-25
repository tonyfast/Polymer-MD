% Analyze_Stats

% Analyze Stats
conv_loc = '../Data/Converted/';
tree_loc = '../Data/Tree/';
part_loc = '../Data/Partitions/';
stats_loc = '../Data/Statistics/';

dd = dir( horzcat( stats_loc, '*.mat') );

bnds = repmat( bsxfun( @times, ones(3,1), [-Inf Inf]), [ 1 1 5]);

for mult = [  1 ]
    for ii = 1 : numel( dd );
        load( horzcat( stats_loc, dd(ii).name), horzcat('X2_',num2str(mult)) );
        load( horzcat( part_loc, dd(ii).name), horzcat('X_',num2str(mult)) );
        
        eval(horzcat('X = ',horzcat('X_',num2str(mult))));
        
        bndstmp = [cellfun( @min, X)', cellfun( @max, X)'];
        bnds(:,:,mult) = [ max( bnds(:,1,mult), bndstmp(:,1) ), min( bnds(:,2,mult), bndstmp(:,2) )];
        
        clear( horzcat('*_',num2str(mult) ) );
    end
end
%%

bndsgrid = fix( bnds  ) - 2 * sign(bnds);
bnds_idx = arrayfun( @(x,y) -1* min(abs(x),abs(y)) : min(abs(x),abs(y)), bndsgrid(:,1), bndsgrid(:,2), UO, false );

%%

for mult = [ 1]
    for ii = 1:numel( dd );
        load( horzcat( conv_loc, dd(ii).name), 'lims', 't_ct', 'fn', 'A' );
        AA{ii} = A;
        if ii == 1 fn_bank{1} = fn; fn_indx = 1; else
            fn_indx = find( cellfun( @(x)strcmp( fn, x ), fn_bank) );
            if numel( fn_indx ) == 0
                fn_bank{ numel( fn_bank ) + 1 } = fn;
                fn_indx = numel( fn_bank ) + 1;
            end
        end
        info( ii,: ) = [ fn_indx t_ct];
        load( horzcat( stats_loc, dd(ii).name), horzcat('F2_',num2str(mult)), horzcat('X2_',num2str(mult)) );
        load( horzcat( part_loc, dd(ii).name), horzcat('X_',num2str(mult)) );
        
%         eval( horzcat( 'X = cellfun( @(x,y) x.* mean( diff( y ) ),', ...
%             'X2_',num2str(mult),', X_',num2str(mult), ' ,UO, false);' )  );
        eval( horzcat( 'X2 = ', 'X_',num2str(mult) ,';'));        
        eval( horzcat( 'X = ', 'X2_',num2str(mult) ,';'));        
        eval( horzcat( 'F = ', 'F2_',num2str(mult) ,';'));
        
        X = cellfun( @(x)x(1:(end-1)), X, UO, false);  % This is cause I saved the data funny

        bndstmp = [cellfun( @min, X)', cellfun( @max, X)'];
%         bnds = [ min( bnds(:,1), bndstmp(:,1) ), max( bnds(:,2), bndstmp(:,2) )];
        


        B = cellfun( @(x,y) x > bnds( y,1, mult ) & x < bnds( y,2, mult ), X, {1,2,3} , UO, false);
        
        F2 = F( find(B{1}), find(B{2}), find(B{3}) );
        [R SS T] = ndgrid( X{1}( B{1} ), X{2}( B{2} ), X{3}( B{3} ) );
        [R2 S2 T2] = ndgrid( bnds_idx{1}, bnds_idx{2}, bnds_idx{3} );
        
        Fi = interp3( SS, R,  T, F2, S2, R2, T2 ); 
        if any( isnan(Fi(:)) ) return; end
        if ii == 1;
            P = zeros( numel(dd), numel(Fi) );
        end
        P(ii,:) = Fi(:);
        clear( horzcat('*_',num2str(mult) ) );
    end
    break
end


s = cellfun( @(x)getOutput(@strtok,1,x,'_'),fn_bank,UO,false);
is1 = find(cellfun( @(x)numel(strfind( x, 'tilt6l'))>0, fn_bank))
s = cellfun( @(x)str2num( x(end-[2:-1:0])), s )
i2 = find( s < 100 );
i = find( ismember( info(:,1),i2) );

info(i,1) = 1;
info(i,2) = s(i2);

%%
[ SS, sid] = sortrows( info , [1 2]);
SS = [SS,sid]
M = [ circshift( SS(:,:), [1 0])];
b = SS(:,1) == circshift( SS(:,1), [1 0]);

G = zeros( size(SS,1) );

G(sub2ind( size(G),SS(b,3), M(b,3))) = 1;
G(sub2ind( size(G),M(b,3), SS(b,3))) = 1;

%% polymer plot
load PCA_Polymer.mat info P G

nu = find(info(:,2)>Inf );
vals = setdiff( 1:size(G,1), nu);
P2 = P( vals, any( ~isnan( P ),1) );
P2(:) = bsxfun(@rdivide, P2, max(P2,[],2) );

[U SS V ] = pca( P2, 10 );
U2 = U * SS;
subplot(2,3,[ 1 2 4 5])
gplot3( G(vals, vals), U2(:,1:3) )
hold on
scatter3(  U2(:,1),U2(:,2), U2(:,3), 50,info(vals,1), 'filled' )
text(U2(:,1),U2(:,2), U2(:,3), num2str(info(vals,2)))

hold off

% dcm = datacursormode;
% set(dcm,'UpdateFcn',@(x,y)view_molecules(x,y,U2, AA))
return
%%

subplot(1,2,1);
pcolor( Fi(:,:,bnds_idx{3}==0))
colorbar
axis equal;
subplot(1,2,2);
pcolor( F(:,:,X{3}==0))
colorbar
axis equal;