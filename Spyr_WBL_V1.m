
%% Spyr_WBL_v1
% Fourth iteration of Spyr_Hist

%% Obtaining Data 
%If one wishes to create a new set of Spyrs one must run code in
%Spyr_Hist_V1


%% Alocate space and Load File

%File Location & Filters
ImgDir = 'Dwdd/images/';
load('Dwdd/lists/splits.mat');
filts = 'sp3Filters';

%Number of pictures
N = input('How many pictures? ');


%Number of frames per picture
NN = input('How many sections per picture? ');

%Preallocating Space
%Master_Matrix = zeros(40960,N*NN);
Master_MatrixOG = zeros(40960,N*NN);

%Prealocating space for pyramids
smFSPyr1_4 = zeros(40960,NN);
smFSPyr1_4OG = zeros(40960,NN); 
%%
%Load Mat File
tic
load('SPYR_2000pic_20sects.mat');
toc
%~4min


%% Basis Set: Create a set of randomly generated Filters
% Space Set [1,2,3,4]
%1st element is irrelevant 
%2nd element will determine ORIENTATION
%3rd element will determine SPACE FREQUENCY - this should be used to color
%the graphs later

tic
freqN = 1:6;

%4th element will determine the phase

%Aperture must remain the same 
BasisSet = input('Pick number of filters: ');
BasisSet1 = ones(BasisSet,1);
ap = 6.5;

%Create Random orientation
spSet = 2*pi*rand(BasisSet,1);

%Create Random Spatial Frequence
spFreq = (pi*1/6) + ((pi*4/6)*rand(BasisSet,1));

%Create Random phase set
phaseSet = 2*pi*BasisSet1;

tic
filterSet =  sPyrDer([1 0 100 pi], [1 0 0 100], [1 100 100], ap);

    for i = 1:length(spFreq);
        filterSet = [filterSet sPyrDer([1 spSet(i) spFreq(i) phaseSet(i)], [1 0 0 100], [1 100 100], ap)];
    end
toc
% 14 MIN 

% Permutating filters
% Extracting the Indexes of relevant pixels

OG_circ = filterSet(:,1);
Der_circ = filterSet(:,1)~=0;
 
Der_circ_ind = find(Der_circ);
Der_circ_ind0 = find(Der_circ==0);
figure(11),spyrDisp(Der_circ);

% Filter Permutation
tic

  
 In_perm = randperm(size(Der_circ_ind,1));
 Der_circ_indP = Der_circ_ind(In_perm);
 
 filterSet_p = filterSet;
 
 for r = 1:size(filterSet,2)
     filterSet_p(Der_circ_ind,r) = filterSet(Der_circ_indP,r);
 end
toc

%1 sec
%% Transpose
tic
Master_MatrixT = Master_Matrix.';
toc
%20 secs

%% Passing Spyrs through Mult Filters - Multiplication 
%Matrix multiplication
tic
Rsp_MxA1 = Master_MatrixT * filterSet;
toc
% 70 secs
tic
Rsp_MxA2 = Master_MatrixT * filterSet_p;
toc
% 80 secs
%% Rectification

%Full Wave rectification  
tic
Rsp_MxA1_FWR = abs(Rsp_MxA1);
Rsp_MxA2_FWR = abs(Rsp_MxA2);
toc

%Half Wave Rectification
tic
Rsp_MxA1_HWR1 = Rsp_MxA1;
Rsp_MxA2_HWR1 = Rsp_MxA2;

for i =1:size(Rsp_MxA1_HWR1,2)
    Rsp_MxA1_HWR{i} = Rsp_MxA1_HWR1(Rsp_MxA1_HWR1(:,i)>0,i); 
end
    
 for i =1:size(Rsp_MxA2_HWR1,2)
     Rsp_MxA2_HWR{i} = Rsp_MxA2_HWR1(Rsp_MxA2_HWR1(:,i)>0,i); 
 end
toc    
 

%Squaring rectification
ticsteam
Rsp_MxA1_SQRD = Rsp_MxA1.^2;
Rsp_MxA2_SQRD = Rsp_MxA2.^2;
toc
 %11 secs
%% Histogram data and response size normalization

%Normalization - here we make the greatest response = 1. It is also
%importnt to mention that the first and last responses were removed since
%they are slightly different.
MaxA1 = round(max(max(Rsp_MxA1_FWR(:,2:(size(Rsp_MxA1_FWR,2)-1)))));
MaxA2 = round(max(max(Rsp_MxA1_SQRD(:,2:(size(Rsp_MxA1_SQRD,2)-1)))));
%MaxA3 = round(max(max(Rsp_MxA1_HWR(:,2:(size(Rsp_MxA1_HWR,2)-1)))));

Rsp_MxA1_FWR_N = zeros(size(Rsp_MxA1_FWR));
Rsp_MxA2_FWR_N = zeros(size(Rsp_MxA1_FWR));
Rsp_MxA1_SQRD_N = zeros(size(Rsp_MxA1_FWR));
Rsp_MxA2_SQRD_N = zeros(size(Rsp_MxA1_FWR));


%These two numbers are important and must be kept the same for all
%normalizations
nbinz = 300;
SampleSz = size(Rsp_MxA1_FWR,1);
tic
%FWR
Rsp_MxA1_FWR_N(:,:) = Rsp_MxA1_FWR(:,:)./MaxA1;
Rsp_MxA2_FWR_N(:,:) = Rsp_MxA2_FWR(:,:)./MaxA1;

%SQRD
Rsp_MxA1_SQRD_N(:,:) = Rsp_MxA1_SQRD(:,:)./MaxA2;
Rsp_MxA2_SQRD_N(:,:) = Rsp_MxA2_SQRD(:,:)./MaxA2;
 
% HWR
% for i=1:size(Rsp_MxA1_HWR,2)
%     Rsp_MxA1_HWR_N{i} = Rsp_MxA1_HWR{i}/MaxA1;
%     Rsp_MxA2_HWR_N{i} = Rsp_MxA2_HWR{i}/MaxA1;
% end

toc

%22 secs
%% Obtaining Prob Distributions

% from the histogram function, we normalize and extract Values 
%This method is faster
tic
ProbDistA1_FWR = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
ProbDistA2_FWR = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
ProbDistA1_HWR = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
ProbDistA2_HWR = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
ProbDistA1_SQRD = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
ProbDistA2_SQRD = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
toc

tic
for n = 2:size(Rsp_MxA1_FWR_N,2)
    hprdA1 = histogram(Rsp_MxA1_FWR_N(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA1_FWR(:,n) = hprdA1.Values;
    hprdA2 = histogram(Rsp_MxA2_FWR_N(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA2_FWR(:,n) = hprdA2.Values;
end
toc 
%~40 sec

tic
for n = 2:size(Rsp_MxA1_FWR_N,2)
    hprdA1 = histogram(Rsp_MxA1_SQRD_N(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA1_SQRD(:,n) = hprdA1.Values;
    hprdA2 = histogram(Rsp_MxA2_SQRD_N(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA2_SQRD(:,n) = hprdA2.Values;
end
toc
%~50 sec

%% Check integration is correct 
qq = zeros(size(filterSet,2),6);
tic
for t = 1:size(ProbDistA1_FWR,2)
    qq(t,1) = sum(ProbDistA1_FWR(:,t));
end 

for t = 1:size(ProbDistA1_FWR,2)
    qq(t,2) = sum(ProbDistA2_FWR(:,t));
end 



for t = 1:size(ProbDistA1_FWR,2)
    qq(t,5) = sum(ProbDistA1_SQRD(:,t));
end 

for t = 1:size(ProbDistA1_FWR,2)
    qq(t,6) = sum(ProbDistA2_SQRD(:,t));
end 
toc

%% FITTING WBL to PROBABILITY DISTRIBUTIONS

%This won't accept any zeros so let's first add a constant 
ixf = randi(size(filterSet,2),(size(Rsp_MxA1_FWR_N,2)-1),1);
xbinz = linspace(0,1,nbinz);

% Get all A and B values from Histogram
WBLA1_FWR = zeros(size(Rsp_MxA1_FWR_N,2),2);
WBLA2_FWR = zeros(size(Rsp_MxA2_FWR_N,2),2);
WBLA1_SQRD = zeros(size(Rsp_MxA1_FWR_N,2),2);
WBLA2_SQRD = zeros(size(Rsp_MxA2_FWR_N,2),2);
tic
for p = 2:size(Rsp_MxA1_FWR_N,2)
    try 
        hfitT1 = fitdist(Rsp_MxA1_FWR_N(:,p),'weibull');
    catch me
        continue
    end
    WBLA1_FWR(p,1) = hfitT1.A;
    WBLA1_FWR(p,2) = hfitT1.B;
end

for p = 2:size(Rsp_MxA2_FWR_N,2)
    try 
        hfitT1 = fitdist(Rsp_MxA2_FWR_N(:,p),'weibull');
    catch me
        continue
    end
    WBLA2_FWR(p,1) = hfitT1.A;
    WBLA2_FWR(p,2) = hfitT1.B;
end

for p = 2:size(Rsp_MxA1_FWR_N,2)
    try 
        hfitT1 = fitdist(Rsp_MxA1_SQRD_N(:,p),'weibull');
    catch me
        continue
    end
    WBLA1_SQRD(p,1) = hfitT1.A;
    WBLA1_SQRD(p,2) = hfitT1.B;
end

for p = 2:size(Rsp_MxA1_FWR_N,2)
    try 
        hfitT1 = fitdist(Rsp_MxA2_SQRD_N(:,p),'weibull');
    catch me
        continue
    end
    WBLA2_SQRD(p,1) = hfitT1.A;
    WBLA2_SQRD(p,2) = hfitT1.B;
end
toc

%~1min
%% CREATE weibull plots
nbinz2 = 300;
xbinz2 = linspace(0,1,nbinz2);

tic
wblplotsA1_FWR = zeros(nbinz2,size(Rsp_MxA1_FWR_N,2));
wblplotsA2_FWR = zeros(nbinz2,size(Rsp_MxA1_FWR_N,2));
wblplotsA1_SQRD = zeros(nbinz2,size(Rsp_MxA1_SQRD_N,2));
wblplotsA2_SQRD = zeros(nbinz2,size(Rsp_MxA1_SQRD_N,2));

for p = 1:size(Rsp_MxA1_FWR_N,2)
    wblplotsA1_FWR(:,p)  = wblpdf(xbinz2,WBLA1_FWR(p,1),WBLA1_FWR(p,2));
    wblplotsA2_FWR(:,p)  = wblpdf(xbinz2,WBLA2_FWR(p,1),WBLA2_FWR(p,2));
    wblplotsA1_SQRD(:,p)  = wblpdf(xbinz2,WBLA1_SQRD(p,1),WBLA1_SQRD(p,2));
    wblplotsA2_SQRD(:,p)  = wblpdf(xbinz2,WBLA2_SQRD(p,1),WBLA2_SQRD(p,2));
end
toc
%~1sec

%% Average response Per filter
tic
Rsp_MxA1_FWR_N_Avg = mean(Rsp_MxA1_FWR_N(:,:));
Rsp_MxA2_FWR_N_Avg = mean(Rsp_MxA2_FWR_N(:,:));

Rsp_MxA1_SQRD_N_Avg = mean(Rsp_MxA1_SQRD_N(:,:));
Rsp_MxA2_SQRD_N_Avg = mean(Rsp_MxA2_SQRD_N(:,:));

%This should give us a matrix of 20x3502
toc

%~ 3 secs
%% FWR Scatter Plots with Orientation, SpFreq, and phase

ORI_set = zeros(size(Rsp_MxA1_FWR_N,2),1);
SpFreq_set = zeros(size(Rsp_MxA1_FWR_N,2),1);
Phase_set = zeros(size(Rsp_MxA1_FWR_N,2),1);

spFreq2 = zeros(length(spFreq)+1,1);
spFreq2(2:end) = spFreq; 
%spFreq2(length(spFreq)+2) = 0

spSet2 = zeros(length(spSet)+1,1);
spSet2(2:end) = spSet; 
%spSet2(length(spFreq)+2) = 0;

phaseSet2 = zeros(length(phaseSet)+1,1);
phaseSet2(2:end) = phaseSet; 
%phaseSet2(length(spFreq)+2) = 0;

tic
Scttr_FWR_AVG(:,1) = Rsp_MxA1_FWR_N_Avg;
Scttr_FWR_AVG(:,2) = Rsp_MxA2_FWR_N_Avg;
Scttr_FWR_AVG(:,3) = spSet2;
Scttr_FWR_AVG(:,4) = spFreq2;
Scttr_FWR_AVG(:,5) = phaseSet2;

Scttr_SQRD_AVG(:,1) = Rsp_MxA1_SQRD_N_Avg;
Scttr_SQRD_AVG(:,2) = Rsp_MxA2_SQRD_N_Avg;
Scttr_SQRD_AVG(:,3) = spSet2;
Scttr_SQRD_AVG(:,4) = spFreq2;
Scttr_SQRD_AVG(:,5) = phaseSet2;
toc




%% Prob Distriubtion and Weibull Fitting FWR
%Mag_Factor_FWR = max(max(wblplotsA1_FWR(2:end,1:150)))/max(max(ProbDistA1_FWR(:,:)));
%Mag_Factor_SQRD = max(max(wblplotsA1_SQRD(2:end,1:150)))/max(max(ProbDistA1_SQRD(:,:)));
Mag_Factor_FWR = 300; %this is 300 bc there are 300 bins in the Prob distributions
ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_N,2)),1);
figure(2)
clf
for i = 1:9
    subplot(3,3,i)
    plot(xbinz,ProbDistA1_FWR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWR(:,ixf(i))/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_FWR(:,ixf(i))/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR','Scrambled Filter FWR','WBL fit Scarmbled Filter FWR')
    title(['FWR NvS w WBL Fits. Filter: ',num2str(ixf(i)),', SpFq: ', num2str(spFreq2(ixf(i)))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.3 0 (max(ProbDistA2_FWR(:,ixf(i)))+ 0.2)])
    
end
%% Prob Distriubtion and Weibull Fitting SQRD
Mag_Factor_SQRD = 300; %this is 300 bc there are 300 bins in the Prob distributions

figure(3)
clf
for i = 1:9
    subplot(3,3,i)
    plot(xbinz,ProbDistA1_SQRD(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_SQRD(:,ixf(i))/Mag_Factor_SQRD),':r')
    plot(xbinz,ProbDistA2_SQRD(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_SQRD(:,ixf(i))/Mag_Factor_SQRD),'--b')
    legend('Normal Filter SQRD','WBL fit Normal Filter SQRD','Scrambled Filter SQRD','WBL fit Scarmbled Filter SQRD')
    title(['SQRD NvS w WBL Fits. Filter: ',num2str(ixf(i)),', SpFq: ', num2str(spFreq2(ixf(i)))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.25 0 (max(ProbDistA2_FWR(:,ixf(i)))+ 0.2)])
end

%% DISPLAY FWR Prob Dist and WBL fitting and Filters 

ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_N,2)),1);
figure(4)
clf
for i = 1:5
    subplot(2,5,i)
    plot(xbinz,ProbDistA1_FWR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWR(:,ixf(i))/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_FWR(:,ixf(i))/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR',...
        'Scrambled Filter FWR','WBL fit Scarmbled Filter FWR')
    title({['FWR NvS w Fits ' num2str(ixf(i))],...
        ['SpFq: ' num2str(spFreq2(ixf(i)))]})
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.3 0 (max(ProbDistA2_FWR(:,ixf(i)))+ 0.2)])
    
end  
for i = 1:5
    filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(2,5,i+5)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title({'Filter | SpFq | Ori | Phase ',...
        [num2str(ixf(i)) ' | ' num2str(spFreq2(ixf(i))) ' | ' num2str(spSet2(ixf(i))) ' | ' num2str(phaseSet2(ixf(i)))]})
end

%% DISPLAY SQRD Prob Dist and WBL fitting and Filters 
ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_N,2)),1);
figure(5)
clf
for i = 1:4
    subplot(2,4,i)
    plot(xbinz,ProbDistA1_SQRD(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_SQRD(:,ixf(i))/Mag_Factor_SQRD),':r')
    plot(xbinz,ProbDistA2_SQRD(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_SQRD(:,ixf(i))/Mag_Factor_SQRD),'--b')
    legend('Normal Filter SQRD','WBL fit Normal Filter SQRD',...
        'Scrambled Filter SQRD','WBL fit Scarmbled Filter SQRD')
    title(['SQRD N v S w Fits ', num2str(ixf(i))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.1 0 (max(ProbDistA2_FWR(:,ixf(i)))+ 0.2)])
    
end  
for i = 1:4
    filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(2,4,i+4)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title(['Filter: ', num2str(ixf(i)),...
        ', SpFq: ', num2str(spFreq2(ixf(i))), ...
        ', Ori: ', num2str(spSet2(ixf(i))),...
        ', Phase: ', num2str(phaseSet2(ixf(i)))])
end

%% BOX  PLOTS
figure (6)
clf
subplot(2,2,1)
boxplot([WBLA1_FWR(2:end,2),WBLA2_FWR(2:end,2)],'Labels',{'Normal Filter FWR','Scrambled Filter FWR'})
axis([0 3 0 2.5])
title(['Box Plot FWR NvS WBL_B'])

subplot(2,2,2)
x2 = linspace(0,max(WBLA1_FWR(:,2)+10),size(WBLA1_FWR,1));

scatter(WBLA1_FWR(:,2),WBLA2_FWR(:,2),30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2,'-r')
axis([0.8 1.5 0.5 2.5])
xlabel('FWR Normal Filter')
ylabel('FWR Scrambled Filter')
title(['Scatter Plot FWR NvS WBL_B'])
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','3/6pi','5/6 pi'},'FontSize',18)
c.Label.String = 'SpFreq';

subplot(2,2,3)
boxplot([WBLA1_SQRD(2:end,2),WBLA2_SQRD(2:end,2)],'Labels',{'Normal Filter SQRD','Scrambled Filter SQRD'})
axis([0 3 0 1.5])
title(['Box Plot SQRD NvS WBL_B'])

subplot(2,2,4)
scatter(WBLA1_SQRD(:,2),WBLA2_SQRD(:,2),30,Scttr_SQRD_AVG(:,4),'filled')
hold on
plot(x2,x2,'-r')
axis([0.4 0.8  0.25 1.25])
xlabel('SQRD Normal Filter')
ylabel('SQRD Scrambled Filter')
title(['Scatter Plot SQRD NvS WBL_B'])
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','3/6pi','5/6 pi'},'FontSize',18)
c.Label.String = 'SpFreq';

%% Scatters of AVG Reponses w orientation 

figure(7)
clf

x2 = linspace(0,max(Rsp_MxA1_FWR_N_Avg(:,2)+10),size(Rsp_MxA1_FWR_N_Avg,1));
scatter(Rsp_MxA1_FWR_N_Avg,Rsp_MxA2_FWR_N_Avg,30,Scttr_FWR_AVG(:,3),'filled')
hold on
plot(x2,x2)
%axis([0 max(Rsp_MxA1_FWR_N_Avg) 0 max(Rsp_MxA1_FWR_N_Avg)])
axis([0 0.25 0 0.25])
xlabel('FWR Normal Filter')
ylabel('FWR Scrambled Filter')
title(['Scatter Plot FWR NvS Avg and Orientation'],'FontSize',14)
c = colorbar('Ticks',[0,(0.5*max(Scttr_FWR_AVG(:,3))),max(Scttr_FWR_AVG(:,3))],'TickLabels',{'0','pi','2 pi'},'FontSize',18)
c.Label.String = 'Orientation';

%% Scatter Plot FWR Avgs w Spatial freq
figure(8)
clf
scatter(Scttr_FWR_AVG(:,1),Scttr_FWR_AVG(:,2),30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
%axis([0 max(Rsp_MxA1_FWR_N_Avg) 0 max(Rsp_MxA1_FWR_N_Avg)])
axis([0 0.1 0 0.1])
axis square
xlabel('FWR Normal Filter Avg Response','FontSize',15)
ylabel('FWR Scrambled Filter Avg Response','FontSize',15)
title(['Scatter Plot FWR NvS Avg and Spatial Freq',],'FontSize',18)
c = colorbar('Ticks',[0,(0.5*max(Scttr_FWR_AVG(:,4))),max(Scttr_FWR_AVG(:,4))],'TickLabels',{'0','pi','5/6 pi'},'FontSize',18)
c.Label.String = 'Space Parameter';

%% Scatter Plot FWR Avgs w phase

figure(9)
clf
scatter(Rsp_MxA1_FWR_N_Avg,Rsp_MxA2_FWR_N_Avg,30,Scttr_FWR_AVG(:,5),'filled')
hold on
plot(x2,x2)
%axis([0 max(Rsp_MxA1_FWR_N_Avg) 0 max(Rsp_MxA1_FWR_N_Avg)])
axis([0 0.1 0 0.1])
axis square
xlabel('FWR Normal Filter')
ylabel('FWR Scrambled Filter')
title(['Scatter Plot FWR NvS Avg and Phase'],'FontSize',14)
c = colorbar('Ticks',[0,(0.5 * max(phaseSet)),max(phaseSet)],'TickLabels',{'0','pi','2 pi'},'FontSize',18)
c.Label.String = 'Phase';

%% Denser version of scatter plots
% Here we will pick 20 random responses from each filter
nf=2;
Rsp_MxA1_FWR_N_SS = zeros(nf,size(Rsp_MxA1_FWR_N,2));
Rsp_MxA2_FWR_N_SS = zeros(nf,size(Rsp_MxA1_FWR_N,2));
Rsp_MxA1_SQRD_N_SS = zeros(nf,size(Rsp_MxA1_FWR_N,2));
Rsp_MxA2_SQRD_N_SS = zeros(nf,size(Rsp_MxA1_FWR_N,2));


for i = 1:size(Rsp_MxA1_FWR_N,2)
    ixf2 = randi([2,size(Rsp_MxA1_FWR_N,1)],nf,1);
    Rsp_MxA1_FWR_N_SS(:,i) = Rsp_MxA1_FWR_N(ixf2,i);
    Rsp_MxA2_FWR_N_SS(:,i) = Rsp_MxA2_FWR_N(ixf2,i);
    Rsp_MxA1_SQRD_N_SS(:,i) = Rsp_MxA1_SQRD_N(ixf2,i);
    Rsp_MxA2_SQRD_N_SS(:,i) = Rsp_MxA2_SQRD_N(ixf2,i);
end

Rsp_MxA1_FWR_N_SSv = Rsp_MxA1_FWR_N_SS(:);
Rsp_MxA2_FWR_N_SSv = Rsp_MxA2_FWR_N_SS(:);

SpFreq_set_SSv1 = zeros(size(Rsp_MxA1_FWR_N_SSv,2),1);
for i = 1:size(Rsp_MxA1_FWR_N,2)
    SpFreq_set_SSv1(:,(nf*i)-(nf-1):nf*i) = SpFreq_set(i);
end
SpFreq_set_SSv = SpFreq_set_SSv1/256;

Rsp_Mx_FWR_SSv = zeros(size(Rsp_MxA1_FWR_N_SSv,1),3);
Rsp_Mx_SQRD_SSv = zeros(size(Rsp_MxA1_FWR_N_SSv,1),3);

Rsp_Mx_FWR_SSv(:,1) = Rsp_MxA1_FWR_N_SS(:);
Rsp_Mx_FWR_SSv(:,2) = Rsp_MxA2_FWR_N_SS(:);
Rsp_Mx_FWR_SSv(:,3) = SpFreq_set_SSv;

Rsp_Mx_SQRD_SSv(:,1) = Rsp_MxA1_SQRD_N_SS(:);
Rsp_Mx_SQRD_SSv(:,2) = Rsp_MxA1_SQRD_N_SS(:);
Rsp_Mx_SQRD_SSv(:,3) = SpFreq_set_SSv;

figure(11)
clf
scatter(Rsp_Mx_FWR_SSv(:,1),Rsp_Mx_FWR_SSv(:,2),15,Rsp_Mx_FWR_SSv(:,3))
hold on
plot(x2,x2,'--r','LineWidth',2)
scatter(Rsp_MxA1_FWR_N_Avg,Rsp_MxA2_FWR_N_Avg,50,Scttr_FWR_AVG(:,4),'filled')
%axis([0 max(Rsp_MxA1_FWR_N_Avg) 0 max(Rsp_MxA1_FWR_N_Avg)])
axis([0 0.25 0 0.25])
xlabel('FWR Normal Filter')
ylabel('FWR Scrambled Filter')
title(['Scatter Plot FWR NvS Avg and Spatial Freq',],'FontSize',14)
c = colorbar('Ticks',[0,0.0061,0.0102],'TickLabels',{'0','pi','5/6 pi'},'FontSize',18)
c.Label.String = 'Spatial Freq';

%% Averge Responses of FWR NvS
figure(12)
clf
plot(Rsp_MxA1_FWR_N_Avg)
hold on
plot(Rsp_MxA2_FWR_N_Avg)
legend({'FWR Normal Avg','FWR Shuffled Avg'},'FontSize',18)
title(['FWR Norm v Shuffled Mean Responses']);
xlabel('FWR Normal Filter')
ylabel('FWR Scrambled Filter')
axis([0 2500 0 0.25])

%% WBL B w Spatial freq
figure(13)
clf
%scatter(WBLA1_FWR(700:1400,2),WBLA2_FWR(700:1400,2),30,Scttr_FWR_AVG(700:1400,4),'filled')
scatter(WBLA1_FWR(:,2),WBLA2_FWR(:,2),30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2)
%axis([0 max(Rsp_MxA1_FWR_N_Avg) 0 max(Rsp_MxA1_FWR_N_Avg)])
axis([0 2.5 0 2.5])
xlabel('FWR Normal Filter')
ylabel('FWR Scrambled Filter')
title(['Scatter Plot FWR NvS Avg and Spatial Freq',],'FontSize',14)
c = colorbar('Ticks',[0,0.0061,0.0102],'TickLabels',{'0','pi','5/6 pi'},'FontSize',18)
c.Label.String = 'Spatial Freq';

%% Scatter plots of Standard Deviation Skew and Kurtosis
tic
Rsp_MxA1_FWR_STD1 = std(Rsp_MxA1_FWR_N);
Rsp_MxA2_FWR_STD1 = std(Rsp_MxA2_FWR_N);

Rsp_MxA1_SQRD_STD1 = std(Rsp_MxA1_SQRD_N);
Rsp_MxA2_SQRD_STD1 = std(Rsp_MxA2_SQRD_N);
toc
%%
Rsp_MxA1_FWR_SKW1 = skewness(Rsp_MxA1_FWR_N);
Rsp_MxA2_FWR_SKW1 = skewness(Rsp_MxA2_FWR_N);

Rsp_MxA1_SQRD_SKW1 = skewness(Rsp_MxA1_SQRD_N);
Rsp_MxA2_SQRD_SKW1 = skewness(Rsp_MxA2_SQRD_N);

Rsp_MxA1_FWR_KUR1 = kurtosis(Rsp_MxA1_FWR_N);
Rsp_MxA2_FWR_KUR1 = kurtosis(Rsp_MxA2_FWR_N);

Rsp_MxA1_SQRD_KUR1 = kurtosis(Rsp_MxA1_SQRD_N);
Rsp_MxA2_SQRD_KUR1 = kurtosis(Rsp_MxA2_SQRD_N);
%% dividing out STD

Rsp_MxA1_FWR_STD = Rsp_MxA1_FWR_N./Rsp_MxA1_FWR_STD1;
Rsp_MxA2_FWR_STD = Rsp_MxA2_FWR_N./Rsp_MxA2_FWR_STD1;

%% Probability ditribution of STD divided responses
ProbDistA1_FWR_STD = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));
ProbDistA2_FWR_STD = zeros(nbinz,size(Rsp_MxA1_FWR_N,2));

tic
for n = 2:size(Rsp_MxA1_FWR_N,2)
    hprdA1 = histogram(Rsp_MxA1_FWR_STD(:,n),'BinEdges',[0:10/nbinz:10],'Normalization','probability');
    ProbDistA1_FWR_STD(:,n) = hprdA1.Values;
    hprdA2 = histogram(Rsp_MxA2_FWR_STD(:,n),'BinEdges',[0:10/nbinz:10],'Normalization','probability');
    ProbDistA2_FWR_STD(:,n) = hprdA2.Values;
    
end
toc 

%% Displaying random filter prob dist 
x = randi([0,size(ProbDistA1_FWR,2)],36,1);
figure(14)
clf
for i = 1:4    
    subplot(2,2,i)
    plot(xbinz,ProbDistA1_FWR(:,x(i)),'-b')
    hold on
    plot(xbinz,ProbDistA2_FWR(:,x(i)),'-r')
     title(['FWR Filter: ', num2str(x(i)),...
        ', SpFq: ', num2str(spFreq2(x(i)))])
    legend('Normal','Scrambled')
    axis([0 0.3 0 (max(ProbDistA2_FWR(:,ixf(i)))+ 0.2)])
end
%% Displaying random filters with STD divided
figure(15)
clf
for i = 1:4
    subplot(2,2,i)
    plot(ProbDistA1_FWR_STD(:,x(i)),'-b')
    hold on
    plot(ProbDistA2_FWR_STD(:,x(i)),'-r')
     title(['FWR/STD Filter: ', num2str(x(i)),...
        ', SpFq: ', num2str(256*(Scttr_FWR_AVG(x(i),4)))])
    legend('Normal','Scrambled')
    xlabel('Proportion')
    ylabel('FWR Scrambled Filter')
    %axis([0 0.03 0 (max(ProbDistA1_FWR_STD(:,ixf(i)))+ 0.2)])
end

%% Scatter plot of STDs
figure(16)
clf

subplot(1,2,1)
scatter(Rsp_MxA1_FWR_N_Avg,Rsp_MxA2_FWR_N_Avg,30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([0 0.06 0 0.06])
set(gca,'FontSize',26);
%c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','pi','5/6 pi'},'FontSize',24)
%c.Label.String = 'Space Factor';
title(['FWR Shuffled v Normal Filter AVG Response'],'FontSize',26)
xlabel('Average Response Normal Filter','FontSize',22)
ylabel('Average Response Shuffled Filter','FontSize',22)
%axis equal 


subplot(1,2,2)
scatter(Rsp_MxA1_FWR_STD1,Rsp_MxA2_FWR_STD1,30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([0 0.06 0 0.06])
set(gca,'FontSize',22);
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'pi/6','pi/2','5/6 pi'},'FontSize',22)
c.Label.String = 'Space Factor';
title(['FWR Shuffled v Normal Filter STDs'],'FontSize',26)
xlabel('Standard Deviation Normal Filter','FontSize',22)
ylabel('Standard Deviation Shuffled Filter','FontSize',22)
%axis equal 



%%
subplot(3,1,2)

scatter(Rsp_MxA1_FWR_SKW1,Rsp_MxA2_FWR_SKW1,30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([0.5 3 0.5 3])
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','pi','5/6 pi'},'FontSize',18)
c.Label.String = 'Spatial Freq';
title(['FWR Skews'])

subplot(3,1,3)

scatter(Rsp_MxA1_FWR_KUR1,Rsp_MxA2_FWR_KUR1,30,Scttr_FWR_AVG(:,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([2 12 2 10])
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','pi','5/6 pi'},'FontSize',18)
c.Label.String = 'Spatial Freq';
title(['FWR Kurtosi'])




%% Displaying filters at spots
dots3 = zeros(length(dots2),2);
for i = 1:length(dots2)
    dots3(i,:) = dots2(i).Position;
end

filt_pt = zeros(length(dots3),1);
for i = 1:length(dots3)
    filt_pt(i) = find(dots3(i,1)==Scttr_FWR_AVG(:,1));
end
%%
figure(17)
clf
%subplot(2,2,1:2)
scatter(Scttr_FWR_AVG(2:end,1),Scttr_FWR_AVG(2:end,2),30,Scttr_FWR_AVG(2:end,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
%axis([0 max(Rsp_MxA1_FWR_N_Avg) 0 max(Rsp_MxA1_FWR_N_Avg)])
axis([0 0.06 0 0.06])
xlabel('Mean filter response','FontSize',22,'FontAngle','italic')
ylabel('Mean shuffled filter response','FontSize',22,'FontAngle','italic')
%title(['Scatter Plot FWR NvS Avg and Spatial Freq',],'FontSize',14)
%c = colorbar('Ticks',[-2.6,0],'TickLabels',{'5','2'},'FontSize',26)
%c.Label.String = 'Derivative Scale (Cyc/RF)';
axis square
set(gca,'FontSize',24,'Color','w')

for i = 1:length(dots3)
    centers = [dots3(i,:)];
    radii = 0.001;
 
    % Display the circles.
    viscircles(centers,radii,'LineWidth',4);
end
xtxt = [dots3(:,1)];
ytxt = [dots3(:,2)];


title(['Scatter Plot FWR NvS Avg and Spatial Freq',],'FontSize',18)
c = colorbar('Ticks',[0,(0.5*max(Scttr_FWR_AVG(:,4))),max(Scttr_FWR_AVG(:,4))],'TickLabels',{'0','pi','5/6 pi'},'FontSize',18)
c.Label.String = 'Space Parameter';
% for i = 1:length(dots3)
%     text(xtxt(i),ytxt(i),num2str(i),'Color','black','FontSize',24)
% end
%%
%{
subplot(2,2,2)
scatter(Rsp_MxA1_FWR_STD1(2:end),Rsp_MxA2_FWR_STD1(2:end),30,Scttr_FWR_AVG(2:end,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([0 0.06 0 0.06])
set(gca,'FontSize',22);
c = colorbar('Ticks',[min(Scttr_FWR_AVG(2:end,4)),max(Scttr_FWR_AVG(2:end,4))],'TickLabels',{'3','0.6'},'FontSize',26)
c.Label.String = 'Derivative scale (Cyc/RF)';
%title(['FWR Shuffled v Normal Filter STDs'],'FontSize',26)
xlabel('Filter STD','FontSize',22,'FontAngle','italic')
ylabel('Shuffled filter STD','FontSize',22,'FontAngle','italic')
set(gca,'FontSize',24,'Color','w')

axis square
for i = 1:length(dots3)
    centers = [dots3(i,:)];
    radii = 0.001;
 
    % Display the circles.
    viscircles(centers,radii,'LineWidth',4);
end
xtxt = [dots3(:,1)];
ytxt = [dots3(:,2)];
% for i = 1:length(dots3)
%     text(xtxt(i),ytxt(i),num2str(i),'Color','black','FontSize',24)
% end
%}


%Filter 1
i=1
filterSetMean = mean(filterSet(Der_circ_ind,filt_pt(i)));
filterSetSTD = std(filterSet(Der_circ_ind,filt_pt(i)));
filterSetnorm(Der_circ_ind,filt_pt(i)) = (filterSet(Der_circ_ind,filt_pt(i)) - filterSetMean)/filterSetSTD;

max_filterSet = max(filterSetnorm(Der_circ_ind,filt_pt(i)));
filterSetnorm1(Der_circ_ind,filt_pt(i)) = filterSetnorm(Der_circ_ind,filt_pt(i))/max_filterSet;

filterSetnorm1(Der_circ_ind0,filt_pt(i)) = 0;
subplot(2,2,i+2)
spyrDisp(filterSetnorm1(:,filt_pt(i)));
title({['Filter ' num2str(i)],...
    ['Cycles per receptive field: ' num2str(0.5*pi/spFreq2(filt_pt(2)))]})
set(gca,'FontSize',24,'Color','w')





%Filter 2
i=2
filterSetMean = mean(filterSet(Der_circ_ind,filt_pt(i)));
filterSetSTD = std(filterSet(Der_circ_ind,filt_pt(i)));
filterSetnorm(Der_circ_ind,filt_pt(i)) = (filterSet(Der_circ_ind,filt_pt(i)) - filterSetMean)/filterSetSTD;

max_filterSet = max(filterSetnorm(Der_circ_ind,filt_pt(i)));
filterSetnorm1(Der_circ_ind,filt_pt(i)) = filterSetnorm(Der_circ_ind,filt_pt(i))/max_filterSet;

filterSetnorm1(Der_circ_ind0,filt_pt(i)) = 0;
subplot(2,2,i+2)
spyrDisp(filterSetnorm1(:,filt_pt(i)));
title({['Filter ' num2str(i)],...
    ['Cycles per receptive field: ' num2str(0.5*pi/spFreq2(filt_pt(2)))]})
set(gca,'FontSize',24,'Color','w')
%%

%Prob distributions
subplot(3,2,3)
plot(xbinz,ProbDistA1_FWR(:,filt_pt(1)),'-r','LineWidth',2)
hold on
%plot(xbinz2,(wblplotsA1_FWR(:,filt_pt(i))/Mag_Factor_FWR),':r')
plot(xbinz,ProbDistA2_FWR(:,filt_pt(1)),'-b','LineWidth',2)
%plot(xbinz2,(wblplotsA2_FWR(:,filt_pt(i))/Mag_Factor_FWR),'--b')
legend('Normal Filter','Scrambled Filter')
%title({['Normal vs Shffled Filter ', num2str(filt_pt(1))],...
%    ['Cycles per receptive field: ' num2str(0.5*pi/spFreq2(filt_pt(1)))]})
xlabel('Reponse Magnitude Normalized','FontSize',22,'FontAngle','italic')
ylabel('Response Probability','FontSize',22,'FontAngle','italic')
axis([0 0.1 0 (max(ProbDistA2_FWR(:,filt_pt(2)))+ 0.05)])
axis square
set(gca,'FontSize',24,'Color','w')

subplot(3,2,4)
plot(xbinz,ProbDistA1_FWR(:,filt_pt(2)),'-r','LineWidth',2)
hold on
%plot(xbinz2,(wblplotsA1_FWR(:,filt_pt(i))/Mag_Factor_FWR),':r')
plot(xbinz,ProbDistA2_FWR(:,filt_pt(2)),'-b','LineWidth',2)
%plot(xbinz2,(wblplotsA2_FWR(:,filt_pt(i))/Mag_Factor_FWR),'--b')
legend('Normal Filter','Scrambled Filter')
%title({['Normal vs Shuffled Filter ' num2str(filt_pt(2))],...
%    ['Cycles per receptive field: ' num2str(0.5*pi/spFreq2(filt_pt(2)))]})
xlabel('Reponse Magnitude Normalized','FontSize',22,'FontAngle','italic')
%ylabel('Response Probability','FontSize',22,'FontAngle','italic')
axis([0 0.1 0 (max(ProbDistA2_FWR(:,filt_pt(2)))+ 0.05)])
axis square
set(gca,'FontSize',24,'Color','w')


%%  ABSTRACT Images
tic
%Accessing picture names
i=1;
%Preallocating space
picnames2 = cell(3000,1);

%loading names
for i=1:3000
    picnames2(i) = splits.train_files(i);
end

p = randi([0 3000],1);
picName = picnames2{p};
imOG = imread(strcat(ImgDir,picName));
imG = rgb2gray(im2double(imOG));
toc
%%
%extracting image size
%it is very important to set x/y_min = 1 because pictures start at 1,1
%not 0,0.
imSize = size(imG(:,:,1));
x_min = 1;
x_max = imSize(2)-70;
y_min = 1;
y_max = imSize(1)-70;

imSpyrs = zeros(40960,5);
smFSPyr1 = zeros(40960,5);
smFSPyr1_4 = zeros(40960,1);

figure(100), 
clf
subplot(5,5,1:10)

    imshow(imG)
    hold on 
    
        for i = 1:5
        ppx = randi([x_min,x_max],1);
        ppy = randi([y_min,y_max],1);
        
        imG1 = imcrop(imG,[ppx,ppy,63,63]);
        rectangle('Position',[(ppx),(ppy),63,63 ],'LineWidth',3,'EdgeColor','r')
        text((ppx+5),(ppy+10),num2str(i),'Color','red','FontSize',10)
        [ imSpyr, imPind, imSpyr2 ] = getSpyr(imG1);
        
        smFSPyr1(:,i) = imSpyr2;
    end
    
            
        mean_smFSPyr1 = mean(mean(smFSPyr1));
        std_smFSPyr1 = mean(var(smFSPyr1,0,1));
        
        smFSPyr1_1 = smFSPyr1-mean_smFSPyr1;
        smFSPyr1_2 = smFSPyr1_1/std_smFSPyr1;
        
        max_smFSPyr1_2 = max(smFSPyr1_2);
        smFSPyr1_3 = smFSPyr1_2./max_smFSPyr1_2;
        
        smFSPyr1_4 = abs(smFSPyr1_3);
        %smFSPyr1_4OG(:,i) = abs(smFSPyr1);
        
        imSpyrs = smFSPyr1_4;
    
for i = 1:5
    subplot(5,5,10+i)
    spyrDisp(imSpyrs(:,i));
    title(num2str(i),'Color','red','FontSize',10)
end    

ixf = randi(size(filterSet,2),(size(Rsp_MxA1_FWR_N,2)-1),1);

for i = 1:5
   filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(5,5,i+15)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title({'Filter | SpFq | Ori | Phase ',...
        [num2str(ixf(i)) ' | ' num2str(spFreq2(ixf(i))) ' | ' num2str(spSet2(ixf(i))) ' | ' num2str(phaseSet2(ixf(i)))]})
    
end

for i =1:5
    subplot(5,5,i+20)
    plot(xbinz,ProbDistA1_FWR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWR(:,ixf(i))/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_FWR(:,ixf(i))/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR',...
        'Scrambled Filter FWR','WBL fit Scarmbled Filter FWR')
    title({['FWR NvS w Fits ' num2str(ixf(i))],...
        ['SpFq: ' num2str(spFreq2(ixf(i)))]})
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.3 0 (max(ProbDistA2_FWR(:,ixf(i)))+ 0.2)])
end
    
%% new Image: Make head of bird fit in 70 by 70 box

figure(101), clf

subplot(3,1,1)
imshow(imG(1:300,1:300))
ppx = 33;
ppy = 100;
rectangle('Position',[ppx,ppy,64,64],'LineWidth',3,'EdgeColor','r')
hold on 


img2 = imG;
img3 = img2(ppy:ppy+64,ppx:ppx+64);
img4 = imresize(img3,[64,64]);
[ imSpyr, imPind, imSpyr2 ] = getSpyr(img4);
smFSPyr1 = imSpyr2;

mean_smFSPyr1 = mean(mean(smFSPyr1));
std_smFSPyr1 = mean(var(smFSPyr1,0,1));

smFSPyr1_1 = smFSPyr1-mean_smFSPyr1;
smFSPyr1_2 = smFSPyr1_1/std_smFSPyr1;

max_smFSPyr1_2 = max(smFSPyr1_2);
smFSPyr1_3 = smFSPyr1_2./max_smFSPyr1_2;

smFSPyr1_4 = abs(smFSPyr1_3);
smFSPyr1_4(Der_circ_ind0) = 0;
%smFSPyr1_4OG(:,i) = abs(smFSPyr1);

imSpyrs = smFSPyr1_4;

subplot(3,1,2)
spyrDisp(smFSPyr1);
xlabel('0°               45°                 90°              135°','FontAngle','italic')
ylabel('Fine            Coarse','FontAngle','italic');
set(gca,'FontSize',24,'Color','w')

subplot(3,1,3)
spyrDisp(smFSPyr1_4);
truesize([500 500]);
xlabel('0°               45°                 90°              135°','FontAngle','italic')
ylabel('Fine            Coarse','FontAngle','italic');
set(gca,'FontSize',24)


%% SET OF 4 FILTER 
filterSet_smll = filterSet(:,2:4);
figure(103)
for i = 1:3
    filterSetMean = mean(filterSet_smll(Der_circ_ind,i));
    filterSetSTD = std(filterSet(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(1,4,i)
    spyrDisp(filterSetnorm1(:,i+1));
    
end 
subplot(1,4,4)

%%
%%
%%
%%
%%
%% Distribution for Scale and Orientation selectivitity 


%Aperture must remain the same 
tic
ap = 6.5;

oriSet_r = 2*pi*rand(3000,1);

total = true;
spaceSet = [];
oriSet = [oriSet_r];
space2Set = [];
ori2Set = [];
space3Set = [];

tic
filterSetOR =  sPyrDer([1 0 100 pi], [1 0 0 100], [1 100 100], ap);

    for i = 1:length(oriSet);
        filterSetOR = [filterSetOR sPyrDer([1 0 100 pi], [0 1 oriSet(i) .1], [1 100 100], ap)];
    end
toc
% Filter Permutation
tic

  
 In_perm = randperm(size(Der_circ_ind,1));
 Der_circ_indP = Der_circ_ind(In_perm);
 
 filterSet_pOR = filterSetOR;
 
 for r = 1:size(filterSet,2)
     filterSet_pOR(Der_circ_ind,r) = filterSet(Der_circ_indP,r);
 end
toc    
    
%%
filterSetnormOR = zeros(40960,size(filterSetOR,2));
filterSetnorm1OR = zeros(40960,size(filterSetOR,2));
filterSet_smll = filterSetOR(:,:);
toc
figure(103)
clf
tic

for i = 1:50
    filterSetMean = mean(filterSet_smll(Der_circ_ind,i));
    filterSetSTD = std(filterSet_smll(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet_smll(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    %filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(5,10,i)
    spyrDisp(filterSetnorm1(:,i));
    
end 
toc

%% shuffled
figure(103)
clf
tic
filterSet_smll = filterSet_pOR;

for i = 1:50
    filterSetMean = mean(filterSet_smll(Der_circ_ind,i));
    filterSetSTD = std(filterSet_smll(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet_smll(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    %filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(5,10,i)
    spyrDisp(filterSetnorm1(:,i));
    
end 
toc




%% Passing Spyrs through Mult Filters - Multiplication 
%Matrix multiplication
tic
Rsp_MxA1OR = Master_MatrixT * filterSetOR;
toc
% 70 secs
tic
Rsp_MxA2OR = Master_MatrixT * filterSet_pOR;
toc
% 80 secs
%% Rectification

%Full Wave rectification  
tic
Rsp_MxA1_FWROR = abs(Rsp_MxA1OR);
Rsp_MxA2_FWROR = abs(Rsp_MxA2OR);
toc

 

%Squaring rectification
tic
Rsp_MxA1_SQRDOR = Rsp_MxA1OR.^2;
Rsp_MxA2_SQRDOR = Rsp_MxA2OR.^2;
toc
 %11 secs

%% 
 
MaxA1OR = round(max(max(Rsp_MxA1_FWROR(:,2:(size(Rsp_MxA1_FWROR,2)-1)))));
MaxA2OR = round(max(max(Rsp_MxA1_SQRDOR(:,2:(size(Rsp_MxA1_SQRDOR,2)-1)))));
%MaxA3 = round(max(max(Rsp_MxA1_HWR(:,2:(size(Rsp_MxA1_HWR,2)-1)))));

Rsp_MxA1_FWR_NOR = zeros(size(Rsp_MxA1_FWROR));
Rsp_MxA2_FWR_NOR = zeros(size(Rsp_MxA1_FWROR));
Rsp_MxA1_SQRD_NOR = zeros(size(Rsp_MxA1_FWROR));
Rsp_MxA2_SQRD_NOR = zeros(size(Rsp_MxA1_FWROR));


%These two numbers are important and must be kept the same for all
%normalizations
nbinz = 300;
SampleSz = size(Rsp_MxA1_FWROR,1);
tic
%FWR
Rsp_MxA1_FWR_NOR(:,:) = Rsp_MxA1_FWROR(:,:)./MaxA1OR;
Rsp_MxA2_FWR_NOR(:,:) = Rsp_MxA2_FWROR(:,:)./MaxA1OR;

%SQRD
Rsp_MxA1_SQRD_NOR(:,:) = Rsp_MxA1_SQRDOR(:,:)./MaxA2OR;
Rsp_MxA2_SQRD_NOR(:,:) = Rsp_MxA2_SQRDOR(:,:)./MaxA2OR;

%% Statistical transformations 

% Average response Per filter
tic
Rsp_MxA1_FWR_NOR_Avg = mean(Rsp_MxA1_FWR_NOR(:,:));
Rsp_MxA2_FWR_NOR_Avg = mean(Rsp_MxA2_FWR_NOR(:,:));

Rsp_MxA1_SQRD_NOR_Avg = mean(Rsp_MxA1_SQRD_NOR(:,:));
Rsp_MxA2_SQRD_NOR_Avg = mean(Rsp_MxA2_SQRD_NOR(:,:));

%This should give us a matrix of 20x3502
toc

% Standard Deviaiton of Orientation Derivatives 
tic
Rsp_MxA1_FWR_STD1OR = std(Rsp_MxA1_FWR_NOR);
Rsp_MxA2_FWR_STD1OR = std(Rsp_MxA2_FWR_NOR);

Rsp_MxA1_SQRD_STD1OR = std(Rsp_MxA1_SQRD_NOR);
Rsp_MxA2_SQRD_STD1OR = std(Rsp_MxA2_SQRD_NOR);
toc


%% Orientation Derivatives Reorganization
%{
ORI_setOR = zeros(size(Rsp_MxA1_FWR_N,2),1);
SpFreq_setOR = zeros(size(Rsp_MxA1_FWR_N,2),1);
Phase_setOR = zeros(size(Rsp_MxA1_FWR_N,2),1);

spFreq2 = zeros(length(spFreq)+1,1);
spFreq2(2:end) = spFreq; 
spFreq2(length(spFreq)+2) = 0;

spSet2 = zeros(length(spSet)+1,1);
spSet2(2:end) = spSet; 
spSet2(length(spFreq)+2) = 0;

phaseSet2 = zeros(length(phaseSet)+1,1);
phaseSet2(2:end) = phaseSet; 
phaseSet2(length(spFreq)+2) = 0;

tic
Scttr_FWR_AVG(:,1) = Rsp_MxA1_FWR_N_Avg;
Scttr_FWR_AVG(:,2) = Rsp_MxA2_FWR_N_Avg;
Scttr_FWR_AVG(:,3) = spSet2;
Scttr_FWR_AVG(:,4) = spFreq2;
Scttr_FWR_AVG(:,5) = phaseSet2;

Scttr_SQRD_AVG(:,1) = Rsp_MxA1_SQRD_N_Avg;
Scttr_SQRD_AVG(:,2) = Rsp_MxA2_SQRD_N_Avg;
Scttr_SQRD_AVG(:,3) = spSet2;
Scttr_SQRD_AVG(:,4) = spFreq2;
Scttr_SQRD_AVG(:,5) = phaseSet2;
toc
%}
%%

figure(26)
clf
subplot(1,2,1)
scatter(Rsp_MxA1_FWR_STD1OR,Rsp_MxA2_FWR_STD1OR,30,Scttr_FWR_AVG(1:3001,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([0.0 0.15 0.0 0.1])
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','pi','5/6 pi'},'FontSize',24)
c.Label.String = 'Space Factor';
title(['Full Wave Rectified Scattered v Normal Filter STDs'],'FontSize',18)
xlabel('Standard Deviation Normal Filter','FontSize',12)
ylabel('Standard Deviation Shuffled Filter','FontSize',12)

subplot(1,2,2)

scatter(Rsp_MxA1_FWR_NOR_Avg,Rsp_MxA2_FWR_NOR_Avg,30,Scttr_FWR_AVG(1:3001,4),'filled')
hold on
plot(x2,x2,'--r','LineWidth',2)
axis([0 0.2 0 0.18])
c = colorbar('Ticks',[0,1.507,2.61],'TickLabels',{'0','pi','5/6 pi'},'FontSize',24)
c.Label.String = 'Space Factor';
title(['Full Wave Rectified Scattered v Normal Filter AVG Response'],'FontSize',18)
xlabel('Average Response Normal Filter','FontSize',12)
ylabel('Average Response Shuffled Filter','FontSize',12)

%% Obtaining Prob Distributions

% from the histogram function, we normalize and extract Values 
%This method is faster
tic
ProbDistA1_FWROR = zeros(nbinz,size(Rsp_MxA1_FWR_NOR,2));
ProbDistA2_FWROR = zeros(nbinz,size(Rsp_MxA1_FWR_NOR,2));
ProbDistA1_HWROR = zeros(nbinz,size(Rsp_MxA1_FWR_NOR,2));
ProbDistA2_HWROR = zeros(nbinz,size(Rsp_MxA1_FWR_NOR,2));
ProbDistA1_SQRDOR = zeros(nbinz,size(Rsp_MxA1_FWR_NOR,2));
ProbDistA2_SQRDOR = zeros(nbinz,size(Rsp_MxA1_FWR_NOR,2));
toc

tic
for n = 2:size(Rsp_MxA1_FWR_NOR,2)
    hprdA1OR = histogram(Rsp_MxA1_FWR_NOR(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA1_FWROR(:,n) = hprdA1OR.Values;
    hprdA2OR = histogram(Rsp_MxA2_FWR_NOR(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA2_FWROR(:,n) = hprdA2OR.Values;
end
toc 
%~6.5 min

tic
for n = 2:size(Rsp_MxA1_SQRD_NOR,2)
    hprdA1OR = histogram(Rsp_MxA1_SQRD_NOR(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA1_SQRDOR(:,n) = hprdA1OR.Values;
    hprdA2OR = histogram(Rsp_MxA2_SQRD_NOR(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA2_SQRDOR(:,n) = hprdA2OR.Values;
end
toc
%~14.5 min

%% Check integration is correct 
qq = zeros(size(filterSet,2),6);
tic
for t = 1:size(ProbDistA1_FWROR,2)
    qq(t,1) = sum(ProbDistA1_FWROR(:,t));
end 

for t = 1:size(ProbDistA1_FWROR,2)
    qq(t,2) = sum(ProbDistA2_FWROR(:,t));
end 



for t = 1:size(ProbDistA1_FWROR,2)
    qq(t,5) = sum(ProbDistA1_SQRDOR(:,t));
end 

for t = 1:size(ProbDistA1_FWROR,2)
    qq(t,6) = sum(ProbDistA2_SQRDOR(:,t));
end 
toc

%% FITTING WBL to PROBABILITY DISTRIBUTIONS

%This won't accept any zeros so let's first add a constant 
ixf = randi(size(filterSet,2),(size(Rsp_MxA1_FWR_N,2)-1),1);
xbinz = linspace(0,1,nbinz);

% Get all A and B values from Histogram
WBLA1_FWROR = zeros(size(Rsp_MxA1_FWR_NOR,2),2);
WBLA2_FWROR = zeros(size(Rsp_MxA2_FWR_NOR,2),2);
WBLA1_SQRDOR = zeros(size(Rsp_MxA1_FWR_NOR,2),2);
WBLA2_SQRDOR = zeros(size(Rsp_MxA2_FWR_NOR,2),2);
tic
for p = 2:size(Rsp_MxA1_FWR_NOR,2)
    try 
        hfitT1OR = fitdist(Rsp_MxA1_FWR_NOR(:,p),'weibull');
    catch me
        continue
    end
    WBLA1_FWROR(p,1) = hfitT1OR.A;
    WBLA1_FWROR(p,2) = hfitT1OR.B;
end

for p = 2:size(Rsp_MxA2_FWR_NOR,2)
    try 
        hfitT1OR = fitdist(Rsp_MxA2_FWR_NOR(:,p),'weibull');
    catch me
        continue
    end
    WBLA2_FWROR(p,1) = hfitT1OR.A;
    WBLA2_FWROR(p,2) = hfitT1OR.B;
end

for p = 2:size(Rsp_MxA1_FWR_NOR,2)
    try 
        hfitT1OR = fitdist(Rsp_MxA1_SQRD_NOR(:,p),'weibull');
    catch me
        continue
    end
    WBLA1_SQRD(p,1) = hfitT1OR.A;
    WBLA1_SQRD(p,2) = hfitT1OR.B;
end

for p = 2:size(Rsp_MxA1_FWR_N,2)
    try 
        hfitT1OR = fitdist(Rsp_MxA2_SQRD_NOR(:,p),'weibull');
    catch me
        continue
    end
    WBLA2_SQRDOR(p,1) = hfitT1OR.A;
    WBLA2_SQRDOR(p,2) = hfitT1OR.B;
end
toc

%~1min
%% CREATE weibull plots
nbinz2 = 300;
xbinz2 = linspace(0,1,nbinz2);

tic
wblplotsA1_FWROR = zeros(nbinz2,size(Rsp_MxA1_FWR_NOR,2));
wblplotsA2_FWROR = zeros(nbinz2,size(Rsp_MxA1_FWR_NOR,2));
wblplotsA1_SQRDOR = zeros(nbinz2,size(Rsp_MxA1_SQRD_NOR,2));
wblplotsA2_SQRDOR = zeros(nbinz2,size(Rsp_MxA1_SQRD_NOR,2));

for p = 1:size(Rsp_MxA1_FWR_NOR,2)
    wblplotsA1_FWROR(:,p)  = wblpdf(xbinz2,WBLA1_FWROR(p,1),WBLA1_FWROR(p,2));
    wblplotsA2_FWROR(:,p)  = wblpdf(xbinz2,WBLA2_FWROR(p,1),WBLA2_FWROR(p,2));
    wblplotsA1_SQRDOR(:,p)  = wblpdf(xbinz2,WBLA1_SQRDOR(p,1),WBLA1_SQRDOR(p,2));
    wblplotsA2_SQRDOR(:,p)  = wblpdf(xbinz2,WBLA2_SQRDOR(p,1),WBLA2_SQRDOR(p,2));
end
toc


%% Prob Distriubtion and Weibull Fitting FWR

Mag_Factor_FWR = 300; %this is 300 bc there are 300 bins in the Prob distributions
ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_NOR,2)),1);
figure(2)
clf
for i = 1:9
    subplot(3,3,i)
    plot(xbinz,ProbDistA1_FWROR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWROR(:,ixf(i))/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWROR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_FWROR(:,ixf(i))/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR','Scrambled Filter FWR','WBL fit Scarmbled Filter FWR')
    title(['FWR NvS w WBL Fits. Filter: ',num2str(ixf(i))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 1 0 (max(ProbDistA2_FWROR(:,ixf(i)))+ 0.005)])
    
end
%% Prob Distriubtion and Weibull Fitting SQRD
Mag_Factor_SQRD = 300; %this is 300 bc there are 300 bins in the Prob distributions

figure(3)
clf
for i = 1:9
    subplot(3,3,i)
    plot(xbinz,ProbDistA1_SQRDOR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),':r')
    plot(xbinz,ProbDistA2_SQRDOR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),'--b')
    legend('Normal Filter SQRD','WBL fit Normal Filter SQRD','Scrambled Filter SQRD','WBL fit Scarmbled Filter SQRD')
    title(['SQRD NvS w WBL Fits. Filter: ',num2str(ixf(i))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.5 0 (max(ProbDistA2_SQRDOR(:,ixf(i)))+ 0.005)])
end




%% DISPLAY FWR Prob Dist and WBL fitting and Filters 

ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_NOR,2)),1);
figure(4)
clf
for i = 1:5
    subplot(2,5,i)
    plot(xbinz,ProbDistA1_FWROR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWROR(:,ixf(i))/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWROR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_FWROR(:,ixf(i))/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR',...
        'Scrambled Filter FWR','WBL fit Scarmbled Filter FWR')
    title({['FWR NvS w Fits ' num2str(ixf(i))],...
        ['SpFq: ' num2str(spFreq2(ixf(i)))]})
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 1 0 (max(ProbDistA2_FWROR(:,ixf(i)))+ 0.005)])
    
end  
for i = 1:5
    filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(2,5,i+5)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title({'Filter ',...
        [num2str(ixf(i))]})
end





%% Display Prob Dist and Wbl Plots

ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_NOR,2)),1);
Mag_Factor_SQRD = 300;
figure(5)
clf
for i = 1:4
    subplot(2,4,i)
    plot(xbinz,ProbDistA1_SQRDOR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),':r')
    plot(xbinz,ProbDistA2_SQRDOR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),'--b')
    legend('Normal Filter SQRD','WBL fit Normal Filter SQRD',...
        'Scrambled Filter SQRD','WBL fit Scarmbled Filter SQRD')
    title(['SQRD N v S w Fits ', num2str(ixf(i))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.5 0 (max(ProbDistA2_SQRDOR(:,ixf(i)))+ 0.01)])
    
end  

for i = 1:4
    filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(2,4,i+4)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title(['Filter: ', num2str(ixf(i))])
end




%%
%%
%%
%%
%% Supplemtary material 3 filters and their response distribution 

%Aperture must remain the same 
tic
ap = 6.5;

oriSet_r = 2*pi*rand(3000,1);

total = true;
spaceSet = [pi*1/2];
oriSet = [0];
space2Set = [];
ori2Set = [];
space3Set = [];
    
filterSetSH =  sPyrDer([1 0 100 pi], [1 0 0 100], [1 100 100], ap);

for i = 1:length(oriSet);
    filterSetSH = [filterSetSH sPyrDer([1 spaceSet(i) 2.4 pi/2], [1 0 0 100], [1 100 100], ap)];
end

for i = 1:length(oriSet);
    filterSetSH = [filterSetSH sPyrDer([1 0 100 pi], [0 1 oriSet(i) .1], [1 100 100], ap)];
end

filterSetSH = [filterSetSH sPyrDer([1 0 100 pi], [1 0 0 100], [1 3 1], ap)];
    toc
    
filterSetnormSH = zeros(40960,size(filterSetSH,2));
filterSetnorm1SH = zeros(40960,size(filterSetSH,2));
filterSet_smllSH = filterSetSH(:,:);
toc

tic


In_perm = randperm(size(Der_circ_ind,1));
Der_circ_indP = Der_circ_ind(In_perm);

filterSet_pSH =  zeros(size(filterSetSH));

for r = 1:size(filterSetSH,2)
    filterSet_pSH(Der_circ_ind,r) = filterSetSH(Der_circ_indP,r);
end
toc






%%
figure(104)
clf
tic

for i = 2:4
    filterSetMean = mean(filterSet_smllSH(Der_circ_ind,i));
    filterSetSTD = std(filterSet_smllSH(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet_smllSH(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1SH(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    %filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(3,3,i-1)
    spyrDisp(filterSetnorm1SH(:,i));
    
end

hold on 

for i = 2:4
    filterSetMean = mean(filterSet_pSH(Der_circ_ind,i));
    filterSetSTD = std(filterSet_pSH(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet_pSH(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1SH(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    %filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(2,3,2+i)
    spyrDisp(filterSetnorm1SH(:,i));
    
end 

%% Passing Spyrs through Mult Filters - Multiplication 
%Matrix multiplication
tic
Rsp_MxA1SH = Master_MatrixT * filterSetSH;
toc
% 70 secs
tic
Rsp_MxA2SH = Master_MatrixT * filterSet_pSH;
toc
% 80 secs
%% Rectification

%Full Wave rectification  
tic
Rsp_MxA1_FWRSH = abs(Rsp_MxA1SH);
Rsp_MxA2_FWRSH = abs(Rsp_MxA2SH);
toc

 

%Squaring rectification
tic
Rsp_MxA1_SQRDSH = Rsp_MxA1SH.^2;
Rsp_MxA2_SQRDSH = Rsp_MxA2SH.^2;
toc
 %11 secs

%% Normalization
tic
MaxA1SH = round(max(max(Rsp_MxA1_FWRSH(:,2:(size(Rsp_MxA1_FWRSH,2)-1)))));
MaxA2SH = round(max(max(Rsp_MxA1_SQRDSH(:,2:(size(Rsp_MxA1_SQRDSH,2)-1)))));
%MaxA3 = round(max(max(Rsp_MxA1_HWR(:,2:(size(Rsp_MxA1_HWR,2)-1)))));

Rsp_MxA1_FWR_NSH = zeros(size(Rsp_MxA1_FWRSH));
Rsp_MxA2_FWR_NSH = zeros(size(Rsp_MxA1_FWRSH));
Rsp_MxA1_SQRD_NSH = zeros(size(Rsp_MxA1_FWRSH));
Rsp_MxA2_SQRD_NSH = zeros(size(Rsp_MxA1_FWRSH));


%These two numbers are important and must be kept the same for all
%normalizations
nbinz = 300;
SampleSz = size(Rsp_MxA1_FWRSH,1);
tic
%FWR
Rsp_MxA1_FWR_NSH(:,:) = Rsp_MxA1_FWRSH(:,:)./MaxA1SH;
Rsp_MxA2_FWR_NSH(:,:) = Rsp_MxA2_FWRSH(:,:)./MaxA1SH;

%SQRD
Rsp_MxA1_SQRD_NSH(:,:) = Rsp_MxA1_SQRDSH(:,:)./MaxA2SH;
Rsp_MxA2_SQRD_NSH(:,:) = Rsp_MxA2_SQRDSH(:,:)./MaxA2SH;
toc
%% Statistical transformations 

% Average response Per filter
tic
Rsp_MxA1_FWR_NSH_Avg = mean(Rsp_MxA1_FWR_NSH(:,:));
Rsp_MxA2_FWR_NSH_Avg = mean(Rsp_MxA2_FWR_NSH(:,:));

Rsp_MxA1_SQRD_NSH_Avg = mean(Rsp_MxA1_SQRD_NSH(:,:));
Rsp_MxA2_SQRD_NSH_Avg = mean(Rsp_MxA2_SQRD_NSH(:,:));

%This should give us a matrix of 20x3502
toc

% Standard Deviaiton of Orientation Derivatives 
tic
Rsp_MxA1_FWR_STD1SH = std(Rsp_MxA1_FWR_NSH);
Rsp_MxA2_FWR_STD1SH = std(Rsp_MxA2_FWR_NSH);

Rsp_MxA1_SQRD_STD1SH = std(Rsp_MxA1_SQRD_NSH);
Rsp_MxA2_SQRD_STD1SH = std(Rsp_MxA2_SQRD_NSH);
toc

%% Obtaining Prob Distributions

% from the histogram function, we normalize and extract Values 
%This method is faster
tic
ProbDistA1_FWRSH = zeros(nbinz,size(Rsp_MxA1_FWR_NSH,2));
ProbDistA2_FWRSH = zeros(nbinz,size(Rsp_MxA1_FWR_NSH,2));
ProbDistA1_SQRDSH = zeros(nbinz,size(Rsp_MxA1_FWR_NSH,2));
ProbDistA2_SQRDSH = zeros(nbinz,size(Rsp_MxA1_FWR_NSH,2));
toc

tic
for n = 2:size(Rsp_MxA1_FWR_NSH,2)
    hprdA1SH = histogram(Rsp_MxA1_FWR_NSH(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA1_FWRSH(:,n) = hprdA1SH.Values;
    hprdA2SH = histogram(Rsp_MxA2_FWR_NSH(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA2_FWRSH(:,n) = hprdA2SH.Values;
end
toc 
%~6.5 min

tic
for n = 2:size(Rsp_MxA1_SQRD_NSH,2)
    hprdA1SH = histogram(Rsp_MxA1_SQRD_NSH(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA1_SQRDSH(:,n) = hprdA1SH.Values;
    hprdA2SH = histogram(Rsp_MxA2_SQRD_NSH(:,n),'BinEdges',[0:1/nbinz:1],'Normalization','probability');
    ProbDistA2_SQRDSH(:,n) = hprdA2SH.Values;
end
toc
%~14.5 min


%% FITTING WBL to PROBABILITY DISTRIBUTIONS

%This won't accept any zeros so let's first add a constant 
ixf = randi(size(filterSet,2),(size(Rsp_MxA1_FWR_NSH,2)-1),1);
xbinz = linspace(0,1,nbinz);

% Get all A and B values from Histogram
WBLA1_FWRSH = zeros(size(Rsp_MxA1_FWR_NSH,2),2);
WBLA2_FWRSH = zeros(size(Rsp_MxA2_FWR_NSH,2),2);
WBLA1_SQRDSH = zeros(size(Rsp_MxA1_FWR_NSH,2),2);
WBLA2_SQRDSH = zeros(size(Rsp_MxA2_FWR_NSH,2),2);
tic
for p = 2:size(Rsp_MxA1_FWR_NSH,2)
    try 
        hfitT1SH = fitdist(Rsp_MxA1_FWR_NSH(:,p),'weibull');
    catch me
        continue
    end
    WBLA1_FWRSH(p,1) = hfitT1SH.A;
    WBLA1_FWRSH(p,2) = hfitT1SH.B;
end

for p = 2:size(Rsp_MxA2_FWR_NSH,2)
    try 
        hfitT1SH = fitdist(Rsp_MxA2_FWR_NSH(:,p),'weibull');
    catch me
        continue
    end
    WBLA2_FWRSH(p,1) = hfitT1SH.A;
    WBLA2_FWRSH(p,2) = hfitT1SH.B;
end

for p = 2:size(Rsp_MxA1_FWR_NSH,2)
    try 
        hfitT1SH = fitdist(Rsp_MxA1_SQRD_NSH(:,p),'weibull');
    catch me
        continue
    end
    WBLA1_SQRDSH(p,1) = hfitT1SH.A;
    WBLA1_SQRDSH(p,2) = hfitT1SH.B;
end

for p = 2:size(Rsp_MxA1_FWR_NSH,2)
    try 
        hfitT1SH = fitdist(Rsp_MxA2_SQRD_NSH(:,p),'weibull');
    catch me
        continue
    end
    WBLA2_SQRDSH(p,1) = hfitT1SH.A;
    WBLA2_SQRDSH(p,2) = hfitT1SH.B;
end
toc

%~1min
%% CREATE weibull plots
nbinz2 = 300;
xbinz2 = linspace(0,1,nbinz2);

tic
wblplotsA1_FWRSH = zeros(nbinz2,size(Rsp_MxA1_FWR_NSH,2));
wblplotsA2_FWRSH = zeros(nbinz2,size(Rsp_MxA1_FWR_NSH,2));
wblplotsA1_SQRDSH = zeros(nbinz2,size(Rsp_MxA1_SQRD_NSH,2));
wblplotsA2_SQRDSH = zeros(nbinz2,size(Rsp_MxA1_SQRD_NSH,2));

for p = 1:size(Rsp_MxA1_FWR_NSH,2)
    wblplotsA1_FWRSH(:,p)  = wblpdf(xbinz2,WBLA1_FWRSH(p,1),WBLA1_FWRSH(p,2));
    wblplotsA2_FWRSH(:,p)  = wblpdf(xbinz2,WBLA2_FWRSH(p,1),WBLA2_FWRSH(p,2));
    wblplotsA1_SQRDSH(:,p)  = wblpdf(xbinz2,WBLA1_SQRDSH(p,1),WBLA1_SQRDSH(p,2));
    wblplotsA2_SQRDSH(:,p)  = wblpdf(xbinz2,WBLA2_SQRDSH(p,1),WBLA2_SQRDSH(p,2));
end
toc


%% Prob Distriubtion and Weibull Fitting FWR

Mag_Factor_FWR = 300; %this is 300 bc there are 300 bins in the Prob distributions
%ixf = randi([2,size(filterSetSH,2)],(size(Rsp_MxA1_FWR_NSH,2)),1);
figure(32)
clf

for i = 2:4
    filterSetMean = mean(filterSet_smllSH(Der_circ_ind,i));
    filterSetSTD = std(filterSet_smllSH(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet_smllSH(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1SH(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    %filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(3,3,i-1)
    spyrDisp(filterSetnorm1SH(:,i));
    title(['FWR Norm v Shuffled with Fits. Filter: ',num2str(i-1)],'FontSize',18)
    
end


for i = 2:4
    filterSetMean = mean(filterSet_pSH(Der_circ_ind,i));
    filterSetSTD = std(filterSet_pSH(Der_circ_ind,i));
    filterSetnorm(Der_circ_ind,i) = (filterSet_pSH(Der_circ_ind,i) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,i));
    filterSetnorm1SH(Der_circ_ind,i) = filterSetnorm(Der_circ_ind,i)/max_filterSet;
    
    %filterSetnorm1(Der_circ_ind0,i) = 0;
    subplot(3,3,i+2)
    spyrDisp(filterSetnorm1SH(:,i));
    
end 

for i = 2:4
    subplot(3,3,i+5)
    plot(xbinz,ProbDistA1_FWRSH(:,i),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWRSH(:,i)/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWRSH(:,i),'-b')
    plot(xbinz2,(wblplotsA2_FWRSH(:,i)/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR','Scrambled Filter FWR','WBL fit Scarmbled Filter FWR','FontSize',18)
    xlabel('Reponse Magnitude Normalized','FontSize',18)
    ylabel('Response Probability','FontSize',18)
    axis([0 0.25 0 (max(ProbDistA2_FWRSH(:,i))+ 0.005)])
    set(gca,'FontSize',18)
end

%% Prob Distriubtion and Weibull Fitting SQRD
Mag_Factor_SQRD = 300; %this is 300 bc there are 300 bins in the Prob distributions

figure(3)
clf

for i = 1:9
    subplot(3,3,i)
    plot(xbinz,ProbDistA1_SQRDOR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),':r')
    plot(xbinz,ProbDistA2_SQRDOR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),'--b')
    legend('Normal Filter SQRD','WBL fit Normal Filter SQRD','Scrambled Filter SQRD','WBL fit Scarmbled Filter SQRD')
    title(['SQRD NvS w WBL Fits. Filter: ',num2str(ixf(i))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.5 0 (max(ProbDistA2_SQRDOR(:,ixf(i)))+ 0.005)])
end




%% DISPLAY FWR Prob Dist and WBL fitting and Filters 

ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_NOR,2)),1);
figure(4)
clf
for i = 1:5
    subplot(2,5,i)
    plot(xbinz,ProbDistA1_FWROR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_FWROR(:,ixf(i))/Mag_Factor_FWR),':r')
    plot(xbinz,ProbDistA2_FWROR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_FWROR(:,ixf(i))/Mag_Factor_FWR),'--b')
    legend('Normal Filter FWR','WBL fit Normal Filter FWR',...
        'Scrambled Filter FWR','WBL fit Scarmbled Filter FWR')
    title({['FWR NvS w Fits ' num2str(ixf(i))],...
        ['SpFq: ' num2str(spFreq2(ixf(i)))]})
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 1 0 (max(ProbDistA2_FWROR(:,ixf(i)))+ 0.005)])
    
end  
for i = 1:5
    filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(2,5,i+5)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title({'Filter ',...
        [num2str(ixf(i))]})
end





%% Display Prob Dist and Wbl Plots

ixf = randi([2,size(filterSet,2)],(size(Rsp_MxA1_FWR_NOR,2)),1);
Mag_Factor_SQRD = 300;
figure(5)
clf
for i = 1:4
    subplot(2,4,i)
    plot(xbinz,ProbDistA1_SQRDOR(:,ixf(i)),'-r')
    hold on 
    plot(xbinz2,(wblplotsA1_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),':r')
    plot(xbinz,ProbDistA2_SQRDOR(:,ixf(i)),'-b')
    plot(xbinz2,(wblplotsA2_SQRDOR(:,ixf(i))/Mag_Factor_SQRD),'--b')
    legend('Normal Filter SQRD','WBL fit Normal Filter SQRD',...
        'Scrambled Filter SQRD','WBL fit Scarmbled Filter SQRD')
    title(['SQRD N v S w Fits ', num2str(ixf(i))])
    xlabel('Reponse Magnitude Normalized')
    ylabel('Response Probability')
    axis([0 0.5 0 (max(ProbDistA2_SQRDOR(:,ixf(i)))+ 0.01)])
    
end  

for i = 1:4
    filterSetMean = mean(filterSet(Der_circ_ind,ixf(i)));
    filterSetSTD = std(filterSet(Der_circ_ind,ixf(i)));
    filterSetnorm(Der_circ_ind,ixf(i)) = (filterSet(Der_circ_ind,ixf(i)) - filterSetMean)/filterSetSTD;
    
    max_filterSet = max(filterSetnorm(Der_circ_ind,ixf(i)));
    filterSetnorm1(Der_circ_ind,ixf(i)) = filterSetnorm(Der_circ_ind,ixf(i))/max_filterSet;
    
    filterSetnorm1(Der_circ_ind0,ixf(i)) = 0;
    subplot(2,4,i+4)
    spyrDisp(filterSetnorm1(:,ixf(i)));
    title(['Filter: ', num2str(ixf(i))])
end






    
