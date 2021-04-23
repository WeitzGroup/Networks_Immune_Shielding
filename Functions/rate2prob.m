% Function converting rate to probability per time
% Function written by Joey Leung
function h = rate2prob(rate,Dt)

h = 1-exp(-rate*Dt);

end