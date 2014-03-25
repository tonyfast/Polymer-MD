function An = Periodic_Data_For_Tree( A, lims )
%%
d = size( A,2);

if d == 2
    [ X Y] = meshgrid( [ -1 : 1], [ -1 : 1] );
    shifts = [X(:) Y(:)];
elseif d ==3
    [ X Y Z] = meshgrid( [ -1 : 1], [ -1 : 1], [-1 :1] );
    shifts = [X(:) Y(:) Z(:)]; 
end

n = size(A,1);
An = zeros( size(shifts,1)*n, 3);


dx = diff(lims,[],2);

for ii = 1 : size(shifts,1)
    An( [1:n] + (ii-1)*n,:) = bsxfun( @plus, A,shifts(ii,:) .*dx(:)' );
end