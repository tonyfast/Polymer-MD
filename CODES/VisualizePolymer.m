function names = VisualizePolymer( data, issave )

maxID = max( data.spatial.chainid )
co = cbrewer( 'qual', 'Paired', maxID);

figure(1);
for ii = 1 : maxID
    ii
    b = data.spatial.chainid == ii;
    h(ii) = plot3(data.spatial.X(b), data.spatial.Y(b), data.spatial.Z(b) );
    if ii == 1
        hold on;
    end
    
    set( h( ii ), 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',co(ii,:),...
        'Markersize',12, 'LineStyle','none');
end
 hold off;


    
    
    grid on;
   
    axis equal;
    axis tight
    ezlabel;
    
for jj = 1 : 4
    
    
    if jj == 4 view(3); 
    else
        v = [0 0 0]; v(jj) = 1;
        view( v );
    end
    
    names{jj} = fullfile('FIGURES', sprintf('%s_%i.png',data.name,jj) ) ;
    if issave
        saveas(1, names{jj}); 
    end
end
figure(gcf)
