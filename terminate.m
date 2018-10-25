function [look, stop ,direction] = terminate (t , state , JD)

r_alt = norm(state(1:3));
look = r_alt - 100;
stop = 1;
direction = -1;

end
