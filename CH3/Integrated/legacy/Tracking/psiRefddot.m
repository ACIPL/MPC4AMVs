function [ psiRddot ] = psiRefddot( xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot )

psiRddot=(1/(xRdot^2+yRdot^2)^2)*((xRddot*yRddot+xRdot*yRdddot-yRddot*xRddot-yRdot*xRdddot)*(xRdot^2+yRdot^2)-(2*xRdot*xRddot+2*yRdot*yRddot)*(xRdot*yRddot-yRdot*xRddot));


end

