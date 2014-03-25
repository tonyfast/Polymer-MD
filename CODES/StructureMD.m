function DATA = StructureMD( filename );

% Pull out timestep
% Number of Atoms
% SImulation Boundaries
% Atom Index
% Chain ID
% Atom Type
% Position X,Y,Z

fo = fopen( filename, 'r' );

ct = 0;
while ~feof( fo )
    
    ct = ct + 1;
    DATA{ct}.name = horzcat( fliplr(strtok(fliplr(filename),'/')) , sprintf('_%i',ct ) );
    s = fgetl(fo);
    DATA{ct}.aggregate.step = fscanf( fo,'%i\n',[1 1]);
    s = fgetl(fo);
    DATA{ct}.aggregate.natom = fscanf( fo,'%i\n',[1 1]);
    s = fgetl(fo);
    temp = fscanf( fo,'%f %f\n',[2 1]);
    
    [ DATA{ct}.aggregate.Xlo, DATA{ct}.aggregate.Xhi ] = deal( temp(1),temp(2) );
    temp = fscanf( fo,'%f %f\n',[2 1]);
    [ DATA{ct}.aggregate.Ylo, DATA{ct}.aggregate.Yhi ] = deal( temp(1),temp(2) );
    temp = fscanf( fo,'%f %f\n',[2 1]);
    [ DATA{ct}.aggregate.Zlo, DATA{ct}.aggregate.Zhi ] = deal( temp(1),temp(2) );
    
    DATA{ct}.aggregate.Volume = abs(DATA{ct}.aggregate.Xhi - DATA{ct}.aggregate.Xlo) .* ...
        abs(DATA{ct}.aggregate.Yhi - DATA{ct}.aggregate.Ylo) .* ...
        abs(DATA{ct}.aggregate.Zhi - DATA{ct}.aggregate.Zlo);
    
    s = fgetl(fo);
    temp = fscanf( fo,'%i %i %i %f %f %f\n', [ 6 DATA{ct}.aggregate.natom])';
    [ DATA{ct}.spatial.index,...
        DATA{ct}.spatial.chainid,...
        DATA{ct}.spatial.atomtrait,...
        DATA{ct}.spatial.X,...
        DATA{ct}.spatial.Y,...
        DATA{ct}.spatial.Z ] = deal( temp(:,1), temp(:,2), temp(:,3), temp(:,4), temp(:,5), temp(:,6) );
    
    if log10(DATA{ct}.aggregate.Volume)>6 | mod( DATA{ct}.aggregate.step, 1e5 )~=0
        % Record data for small values and over 10 time steps
        ct = ct - 1;
    end
end