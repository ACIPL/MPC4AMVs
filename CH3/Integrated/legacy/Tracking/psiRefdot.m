function [ psiRdot ] = psiRefdot( xRdot,xRddot,yRdot,yRddot )

psiRdot=(xRdot*yRddot-yRdot*xRddot)/(xRdot^2+yRdot^2);


end

