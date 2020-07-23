close all;
clear all;
clc;

%init important parameters
addpath('C:\Users\oleskiw\Google Drive\wrx\Program\Toolbox\matlabPyrTools');
addpath('C:\Users\oleskiw\Google Drive\wrx\Program\Toolbox\matlabPyrTools\MEX')

S = 4;
O = 4;
N = 64;

M = 2^11;

FLAG_SCRAMBLE_IMAGE = false;

%%

%demo pyramid
imgFull = imread('.\testImages\bark01.jpg');
imgFull = hsv2rgb(rgb2hsv(imgFull).*repmat(reshape([1,0,1],[1 1 3]), [size(imgFull,1) size(imgFull,2), 1]));
imgFull = squeeze(imgFull(:,:,1));
if FLAG_SCRAMBLE_IMAGE
    %permute pixels
    sx = size(imgFull);
    imgFull = reshape(imgFull(randperm(numel(imgFull))), sx);
else
    %do nothing
end
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


figure('position', [21   477   649   757]);
subplot(2,1,1)
imagesc(imgPatch);
colormap(gray);
axis equal off;

subplot(2,1,2);
spyrTest = spyrDataMatrix(M,:);
spyrDisp4(spyrTest/max(abs(spyrTest(:))), pind, pyr);

%%
% construct a spyr coordinate matrix
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


%%
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


%% display a filter for diagnostics

%lookSet = [1 2 3 4];
lookSet = [3 7 11 15];

for i = 1:length(lookSet)
    clf;
    n = lookSet(i);
    f = filterMatrix(n,:);
    w = filterSupportMatrix(n,:);
    w(w>0) = 1;
    fs = filterShuffleMatrix(n,:);

    subplot(2,2,1);
    spyrDisp4(f/max(f(:)), pind, pyr);
    title('Filter')
    subplot(2,2,2);
    spyrDisp4(w/max(w(:)), pind, pyr);
    title('shuffling window');
    subplot(2,2,3);
    spyrDisp4(fs/max(fs(:)), pind, pyr);
    title('Shuffled Filter');


    %compute histograms
    respBins = linspace(0,1,30);
    nh = histc(responseNormal(:,n), respBins);
    sh = histc(responseShuffle(:,n), respBins);

    subplot(2,2,4); hold on;
    plot(respBins, nh/sum(nh), 'b-');
    plot(respBins, sh/sum(sh), 'r-');
    title('Response Histograms');
    xlabel('Normalized Response');
    ylabel('Frequency');
    legend({'Norm.', 'Shuff.'});
    
    drawnow();
    pause();  
end

    % 
    % %plot
    % gWindow = d.window;
    % gSinusoid = d.diff;
    % gDerivOperator = d.op;
    % subplot(3,1,1);
    % spyrDisp4(gWindow/max(gWindow(:)), pind, pyr);
    % subplot(3,1,2);
    % spyrDisp4(gSinusoid/(max(gSinusoid(:))), pind, pyr);
    % subplot(3,1,3);
    % spyrDisp4(gDerivOperator/(max(gDerivOperator(:))), pind, pyr);


