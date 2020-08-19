%% SPYR V3: 4x4


%init important parameters
%addpath('C:\Users\oleskiw\Google Drive\wrx\Program\Toolbox\matlabPyrTools');
addpath('/arc/1.2/p3/ruidiazp/Documents/MATLAB/matlabPyrTools2/PyrTools/V2DerivFilterSparsity');

%addpath('C:\Users\oleskiw\Google Drive\wrx\Program\Toolbox\matlabPyrTools\MEX')

S = 4;
O = 4;
N = 64;

M = 2^11;

%% This scrambles the image: 1 - permute, 0 - normal
FLAG_SCRAMBLE_IMAGE = 0;

%%
%demo pyramid
imgFull = imread('./testImages/bark01.jpg');
imgFull = hsv2rgb(rgb2hsv(imgFull).*repmat(reshape([1,0,1],[1 1 3]), [size(imgFull,1) size(imgFull,2), 1]));
imgFull = squeeze(imgFull(:,:,1));

%image permutation
if FLAG_SCRAMBLE_IMAGE == 1
    %permute pixels
    sx = size(imgFull);
    imgFull = reshape(imgFull(randperm(numel(imgFull))), sx);
else
    %do nothing
end

%This cuts many pictures and gets their spyr representation
spyrDataMatrix = nan(M,3);
for m = 1:M
    x0 = [randi(size(imgFull,1)-(N+1)) randi(size(imgFull,2)-(N+1))];
    imgPatch = imgFull(x0(1):x0(1)+N-1, x0(2):x0(2)+N-1);
    imgPatch = imgPatch - 0.5;

    [spyrTest, pind, pyr] = getSpyr4(1*imgPatch, S, O);
    spyrDataMatrix(m,1:length(spyrTest)) = spyrTest;
end
%rectify
spyrDataMatrix = abs(spyrDataMatrix);


%figure('position', [21   477   649   757]);
figure(3)
subplot(2,1,1)
imagesc(imgPatch);
colormap(gray);
axis equal off;

subplot(2,1,2);
spyrTest = spyrDataMatrix(M,:);
spyrTest2 = spyrTest/max(abs(spyrTest(:)));
spyrTest3 = spyrTest2-mean(spyrTest2);
spyrDisp4(spyrTest3,pind,pyr)
%spyrDisp4(spyrTest/max(abs(spyrTest(:))), pind, pyr);

%% construct a spyr coordinate matrix
% each column is a spyr element, each row a dimension
spyrCoords = V2DerivFilterCoords(pind);
coordPosX = spyrCoords(1,:); 
coordPosY = spyrCoords(2,:); 
coordScale = spyrCoords(3,:);
coordOri = spyrCoords(4,:);
coordRes = spyrCoords(5,:);

% figure('position', [24          76        1644         966]);
% subplot(2,3,1);
% spyrDisp4(coordPosX/(max(abs(coordPosX(:)))), pind, pyr);
% title('X');
% subplot(2,3,2);
% spyrDisp4(coordPosY/(max(abs(coordPosY(:)))), pind, pyr);
% title('Y');
% subplot(2,3,3);
% spyrDisp4(coordScale/(max(abs(coordScale(:)))), pind, pyr);
% title('Scale');
% subplot(2,3,4);
% spyrDisp4(coordOri/(max(abs(coordOri(:)))), pind, pyr);
% title('Ori');
% subplot(2,3,5);
% spyrDisp4(coordRes/(max(abs(coordRes(:)))), pind, pyr);
% title('Residual');
% 
% 
% spyrCoords = [coordPosX; coordPosY; coordScale; coordOri; coordRes];


clf;

filterMatrix = [];
filterSupportMatrix = [];

fDimA = [0 45 90 135];
fDimB = [90 45 0 135];

fCount = 0;
filterCoordA = nan(1,numel(fDimA)*numel(fDimB));
filterCoordB = nan(1,numel(fDimA)*numel(fDimB));
for i = 1:length(fDimA)
    for j = 1:length(fDimB);
        fCount = fCount + 1;
        filterCoordA(fCount) = fDimA(i);
        filterCoordB(fCount) = fDimB(j);

        %where the derivative is to be evaluated
        gMu = [32 32 2 fDimB(j)];

        %direction of the derivative, normalized;
        gDirection = [cos(deg2rad(fDimA(i))) sin(deg2rad(fDimA(i))) 0 0];
        gDirection = gDirection / norm(gDirection);

        %std weights off/on direction vector
        gScale = .01;
        freq = 2;
        sigma(1) = gScale';
        sigma(2:4) = 0.2*gScale' * [1 1 1];

        %'size' of each dimension
        dimScale = ([N, N, 5, -180]);

        d = V2DerivFilter(spyrCoords, gMu, gDirection, dimScale, sigma, freq);
        filterMatrix(fCount,:) = d.op;
        filterSupportMatrix(fCount,:) = d.window;
    end
end


%% shuffle matrices
filterShuffleMatrix = filterMatrix;
for i = 1:size(filterMatrix,1)
    w = filterSupportMatrix(i,:);
    f = filterMatrix(i,:);
    fs = f(w>0);
    filterShuffleMatrix(i,w>0) = fs(randperm(length(fs)));
end

%% multiply filter matrix with data matrix
responseNormal = spyrDataMatrix * filterMatrix';
responseShuffle = spyrDataMatrix * filterShuffleMatrix';

%normalize by global response maximum
grm = max([responseNormal(:); responseShuffle(:)]);
responseNormal = responseNormal / grm;
responseShuffle = responseShuffle / grm;









