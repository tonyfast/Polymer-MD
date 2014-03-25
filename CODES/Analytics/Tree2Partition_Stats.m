% Partition the spatial domain and compute the statistics for each
% partition

conv_data_loc = '../Data/Converted/';
tree_data_loc = '../Data/tree/';

part_data_loc = '../Data/Partitions/';
stats_data_loc = '../Data/Statistics/';

ff = dir( horzcat( conv_data_loc,'*.mat') );
ll = [Inf*ones(3,1), Inf*ones(3,1)];
for qq = 1 : numel(ff)
    
    load( horzcat( conv_data_loc, ff(qq).name), 'lims');
    ll = [ min( abs(lims(:,1)), ll(:,1)) min( abs(lims(:,2)), ll(:,2)) ];
end

ff = dir( horzcat( tree_data_loc,'*.mat') );



%%


N = 7999; % NUmber of particles
mulv = [ 1 2 3 .25 .5];  % size of the partitions
for mult = [ 1 ];;
    
    for qq = 108:numel( ff )
        
        load( horzcat( tree_data_loc, ff(qq).name), 'tr');
        load( horzcat( conv_data_loc, ff(qq).name), 'lims','A' );
        
        
        V = ceil(diff(lims,[],2))*mult + abs(mod(ceil(diff(lims,[],2)*mult),2)+1);
        if all(V<500) % I forget what this flag was for
            for zz = 1 : 3;
                X{zz} = linspace( lims(zz,1), lims(zz,2), (V(zz)-1)./2 + 1);
                dx(zz,:) = bsxfun( @times, mean(diff(X{zz})), [-1 1]./2);
            end
            
            %     X{1} =
            tr = kdtree_io_from_mat( tr );
            
            [ R S T ] = ndgrid( X{1},X{2}, X{3} );
            
            % Partition the space and index the points in each partition
            I = arrayfun( @(x,y,z) getOutput(@kdtree_range_query, 1, tr , [x + dx(1,:); y + dx(2,:); z + dx(3,:)] )', R(:), S(:), T(:), UO, false );
            I = cellfun( @(x)mod(x,N) + (mod(x,N)==0)*N, I, UO, false );
            
            
            kdtree_delete( tr );
            pause(.5)
            B = zeros( size(R));
            B(:) = cellfun( @numel, I);
            
            F2 = fftshift(ifftn( fftn( B ) .* conj( fftn( B) ) ));
            a = zeros(1,3);
            [ a(1) a(2) a(3)] = ind2sub( size(B),getOutput(@max,2,F2(:)))
            
            for zz = 1 : 3;
                X2{zz} = [ 0 : numel(X{zz})]*dx(zz,2)*2 - (a(zz)-1)*dx(zz,2)*2;
            end
            
            assignin( 'base', horzcat('I_',num2str(find(  mult==mulv ))), I )
            assignin( 'base', horzcat('X_',num2str(find(  mult==mulv ))), X )
            
            if exist( horzcat( part_data_loc, ff(qq).name), 'file')
                save( horzcat( part_data_loc, ff(qq).name), horzcat('I_',num2str(find(  mult==mulv ))), horzcat('X_',num2str(find(  mult==mulv ))),'-append');
            else
                save( horzcat( part_data_loc, ff(qq).name), horzcat('I_',num2str(find(  mult==mulv ))), horzcat('X_',num2str(find(  mult==mulv ))));
            end
            clear( horzcat('I_',num2str(find(  mult==mulv ))), horzcat('X_',num2str(find(  mult==mulv ))) );
            
            assignin( 'base', horzcat('F2_',num2str(find(  mult==mulv ))), F2 )
            assignin( 'base', horzcat('X2_',num2str(find(  mult==mulv ))), X2 )
            if exist( horzcat( stats_data_loc, ff(qq).name), 'file')
                save( horzcat( stats_data_loc, ff(qq).name), horzcat('F2_',num2str(find(  mult==mulv ))), horzcat('X2_',num2str(find(  mult==mulv ))),'-append');
            else
                save( horzcat( stats_data_loc, ff(qq).name), horzcat('F2_',num2str(find(  mult==mulv ))), horzcat('X2_',num2str(find(  mult==mulv ))));
            end
            clear( horzcat('F2_',num2str(find(  mult==mulv ))), horzcat('X2_',num2str(find(  mult==mulv ))) );
        end
    end
end