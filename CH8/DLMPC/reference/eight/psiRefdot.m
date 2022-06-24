function [ psiRdot ] = psiRefdot( xRdot,xRddot,yRdot,yRddot )

% xRdot=xRdot(t);
% xRddot=xRddot(t);
% yRdot=yRdot(t);
% yRddot=yRddot(t);

psiRdot=(xRdot*yRddot-yRdot*xRddot)/(xRdot^2+yRdot^2);


end

