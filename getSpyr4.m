function [imSpyr, pind, pyr] = getSpyr4(im, S, O)
%GETSPYR wrapper function for buildSFspyr
%splits positive and negative channels

N = size(im,1);
[spyr, pind] = buildSFpyr(im, S, O-1);
z = 0;
for r = 1:size(pind, 1);
    zNext = prod(pind(r,:));

    f = spyr(z+1:z+zNext);
    n = ((N+2)/pind(r,1))^2;
    spyr(z+1:z+zNext) = f/n;        

    z = z+zNext;
end
imSpyr = spyr;

% divide channel amplitudes
pyr = [S O N];

end
