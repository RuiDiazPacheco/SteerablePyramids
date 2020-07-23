function [ deriv ] = V2DerivFilter(spyrCoords, mu, direction, scale, sigma, f)
    gMu = mu;
    gDirection = direction;
    gFrequency = f;
    dimScale = scale;
    sigmaOn = sigma(1);
    sigmaOff = sigma(2:4);

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
    gWindow(gWindow < 1e-5) = 0;

    %project coordinates onto 
    derivativeDirection = gDirection .* dimScale;
    coordProj = (spyrCoords(1:4,:)' - repmat(gMu, [length(spyrCoords), 1])) * derivativeDirection';
    gSinusoid = sin(gFrequency*coordProj*((2*pi)) / norm(derivativeDirection).^2);
    gDerivOperator = gWindow .* gSinusoid;

    %normalize operator to have unit connection strengths
    deriv.op = gDerivOperator / sum(abs(gDerivOperator(:)));
    deriv.window = gWindow;
    deriv.diff = gSinusoid;
end