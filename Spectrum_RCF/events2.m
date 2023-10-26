function [value,isterminal,direction] = events2(x,E)

% Locate the position when energy passes through zero and stop integration.

global epaisseur_totale

value(1) = E;
isterminal(1) = 1;   % Stop the integration
direction(1) = 0;   % all directions

value(2) = x-epaisseur_totale;
isterminal(2) = 1;   % Stop the integration
direction(2) = 0;   % all directions


