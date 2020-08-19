function [ coords ] = V2DerivFilterCoords(pind)
%V2DERIVFILTERCOORDS generate coordinate functions for spyr derivatives


coordPosX = nan(1,sum(prod(pind,2)));
coordPosY = nan(1,sum(prod(pind,2)));
coordScale = nan(1,sum(prod(pind,2)));
coordOri = nan(1,sum(prod(pind,2)));
coordRes = zeros(1,sum(prod(pind,2)));

idx = cumsum(prod(pind,2));
idx = [0; idx];
oriNum = sum(sum(pind == pind(1,1),1))/2-1;
oriSet = 90:-180/oriNum:-89;
n = 1;
s = (length(pind)-2)/oriNum;
N = pind(1,1);
for z = 1:length(pind)
        
        %set x and y positional coords
        
        S = N/pind(z,1);
        dspan = ceil(S/2):S:N-ceil(S/2)+1;
        [Dc, Dr] = meshgrid(dspan, dspan);

        coordPosX(idx(z)+1:idx(z+1)) = Dc(:);
        coordPosY(idx(z)+1:idx(z+1)) = Dr(:);

        %do ori
        coordOri(idx(z)+1:idx(z+1)) = oriSet(n)*ones(size(Dc));
        
        %do scale
        coordScale(idx(z)+1:idx(z+1)) = s*ones(size(Dc));
        
        if z == 1 || z == length(pind);
            %residual bands
            coordRes(idx(z)+1:idx(z+1)) = ones(pind(z,1), pind(z,2));
            coordRes(idx(z)+1:idx(z+1)) = ones(pind(z,1), pind(z,2));
        else
            n = n + 1;
            if n == oriNum+1
                n = 1;
                s = s-1;
            end
        end
end




coords = [coordPosX; coordPosY; coordScale; coordOri; coordRes];
end

