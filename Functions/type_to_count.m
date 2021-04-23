%converts a vector with values 0, 1, 2 or 3 (standing for S, E, I or R) to
%the number of S, E, I and R
function r = type_to_count(G_type)
    %vector to count the number of S, E, I and R individuals
    G_count = zeros(4,1);
    for i = 1:4
        G_count(i) = sum(G_type == i-1);
    end
    r =  G_count;
end