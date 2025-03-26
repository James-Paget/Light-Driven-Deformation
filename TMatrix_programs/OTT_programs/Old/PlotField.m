function [X, Y, Zcol] = PlotField(field, z_plane, plotSize, plotN, shouldLogarithm, shouldPlot)
    % Plots the magnitude of the field in the z=z_plane plane and returns
    % the X,Y,Z matrices used to make the surface. plotSize sets the x,y
    % limits and plotN gives the number of sampled points on each axis.

    linX = linspace(-plotSize, plotSize, plotN);
    linY = linspace(-plotSize, plotSize, plotN);
    Zcol = zeros(plotN, plotN);
    [X,Y] = meshgrid(linX, linY);
    
    positions = zeros(plotN^2,3);
    i=1;
    for x = linX
        for y = linY
            positions(i,:) = [x,y,z_plane];
            i = i+1;
    
        end
    
    end
    
    EMfield = field.emFieldXyz(positions.').';
    Zcol(:,:) = vecnorm(reshape(EMfield, [plotN, plotN,3]),2,3);
    
    if shouldLogarithm
        Zcol(:,:) = log(Zcol(:,:));
    end
    
    if shouldPlot
        figure
        surf(X, Y, Zcol)
        shading interp
        xlabel("x")
        ylabel("y")
        zlabel("Field")
        % View can set it to look down from above for a 2D plot.
        view(0,90)
    end
end