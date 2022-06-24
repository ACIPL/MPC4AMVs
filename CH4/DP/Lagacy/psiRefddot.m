function [ psiRddot ] = psiRefddot( xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot )

% xRdot=xRdot(t);
% xRddot=xRddot(t);
% xRdddot=xRdddot(t);
% 
% yRdot=yRdot(t);
% yRddot=yRddot(t);
% yRdddot=yRdddot(t);

% psiRddot=(1/(xRdot^2+yRdot^2)^2)*((xRddot*yRddot+xRdot*yRdddot-yRddot*xRddot-yRdot*xRdddot)*(xRdot^2+yRdot^2)-(2*xRdot*xRddot+2*yRdot*yRddot)*(xRdot*yRddot-yRdot*xRddot));

psiRddot=0;

end

