function [ psiR ] = psiRef( xRdot,yRdot )

% yRdot=yRdot(t);
% xRdot=xRdot(t);

psiR=atan2(yRdot,xRdot);

end

