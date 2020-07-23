close all;
clear all;
clc;

%init important parameters
addpath('C:\Users\oleskiw\Google Drive\wrx\Program\Toolbox\matlabPyrTools');
addpath('C:\Users\oleskiw\Google Drive\wrx\Program\Toolbox\matlabPyrTools\MEX')

S = 5;
O = 6;
N = 128;

%%

%demo pyramid
imgTest = imread('C:\Users\Oleskiw\Google Drive\wrx\Program\textureShorts\coefficientBundle\images\bark01.jpg');
imgTest = hsv2rgb(rgb2hsv(imgTest).*repmat(reshape([1,0,1],[1 1 3]), [size(imgTest,1) size(imgTest,2), 1]));
imgTest = squeeze(imgTest(:,:,1));

rng(1337);
x0 = randi(N,1,2);

imgTest = imgTest(x0(1):x0(1)+N-1, x0(2):x0(2)+N-1);
imgTest = imgTest - 0.5;


[spyrTest, pind, pyr] = getSpyr4(1*imgTest, S, O);

figure('position', [21   477   649   757]);
subplot(2,1,1)
imagesc(imgTest);
colormap(gray);
axis equal off;

subplot(2,1,2);
spyrDisp4(spyrTest/max(abs(spyrTest(:))), pind, pyr);

%%
% construct a spyr coordinate matrix
% each column is a spyr element, each row a dimension

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

figure('position', [24          76        1644         966]);
subplot(2,3,1);
spyrDisp4(coordPosX/(max(abs(coordPosX(:)))), pind, pyr);
title('X');
subplot(2,3,2);
spyrDisp4(coordPosY/(max(abs(coordPosY(:)))), pind, pyr);
title('Y');
subplot(2,3,3);
spyrDisp4(coordScale/(max(abs(coordScale(:)))), pind, pyr);
title('Scale');
subplot(2,3,4);
spyrDisp4(coordOri/(max(abs(coordOri(:)))), pind, pyr);
title('Ori');
subplot(2,3,5);
spyrDisp4(coordRes/(max(abs(coordRes(:)))), pind, pyr);
title('Residual');


spyrCoords = [coordPosX; coordPosY; coordScale; coordOri; coordRes];


%%
clf;
%where the derivative is to be evaluated
gMu = [64 64 3.5 30];

%direction of the derivative, normalized;
gDirection = [0 0 1 0];
gDirection = gDirection / norm(gDirection);

%std weights off/on direction vector
gScale = .025;
gFrequency = 2;
sigmaOn = gScale';
sigmaOff = 0.33*gScale' * [1 1 1];

%'size' of each dimension
dimScale = ([N, N, 5, -180]);

%create eigenvectors and values for covariance matrix
eVec = [gDirection' null(gDirection)];
eVec = eVec .* repmat(dimScale', [1 size(eVec,2)]);
eVal = diag([sigmaOn, sigmaOff]);
gSigma = eVec * eVal * eVec';

%transform ori coords to enforce periodoc domain in ori
spyrCoords2 = spyrCoords(1:4,:);
spyrCoords2A = spyrCoords2;
spyrCoords2A(4,:) = spyrCoords(4,:) - -180;
spyrCoords2B = spyrCoords2;
spyrCoords2B(4,:) = spyrCoords(4,:) - 0;
spyrCoords2C = spyrCoords2;
spyrCoords2C(4,:) = spyrCoords(4,:) - 180;


% construct gaussian window
gWindow = mvnpdf(spyrCoords2A', gMu, gSigma) + mvnpdf(spyrCoords2B', gMu, gSigma) + mvnpdf(spyrCoords2C', gMu, gSigma);
gWindow(gWindow < 1e-8) = 0;

%project coordinates onto 
derivativeDirection = gDirection .* dimScale;
coordProj = (spyrCoords(1:4,:)' - repmat(gMu, [length(spyrCoords), 1])) * derivativeDirection';
gSinusoid = sin(gFrequency*coordProj*((2*pi)) / norm(derivativeDirection).^2);
gDerivOperator = gWindow .* gSinusoid;

%normalize operator to have unit connection strengths
gDerivOperator = gDerivOperator / sum(abs(gDerivOperator(:)));

%plot
subplot(3,1,1);
spyrDisp4(gWindow/max(gWindow(:)), pind, pyr);
subplot(3,1,2);
spyrDisp4(gSinusoid/(max(gSinusoid(:))), pind, pyr);
subplot(3,1,3);
spyrDisp4(gDerivOperator/(max(gDerivOperator(:))), pind, pyr);


