 %% New SPYR 2020 
% This uses the new filter of 5 levels

%% Alocate space and Load File

%% File Location & F
ImgDir3 = '/arc/1.2/p3/ruidiazp/Documents/MATLAB/matlabPyrTools2/PyrTools/Dwdd/images/Textures/';
ImgDir2 = dir(ImgDir3);

ImgDir = cell(((length(ImgDir2)/2)-2),1);

for n = ((length(ImgDir2)/2)+3):length(ImgDir2)
    ImgDir{n-552} = ImgDir2(n).name;
end

% check images are being read properly 
xIm = randi((length(ImgDir2)/2)-2,1);
figure(2),clf
imshow([ImgDir3,'/',ImgDir{xIm}]);

%% Creating new Spyrs

filts = 'sp3Filters';

% Number of pictures
N1 = input('How many pictures? ');

% Number of frames per picture
N1s = input('How many sections per picture? ');

% Preallocating Space
% Master_Matrix = zeros(40960,N*NN);
Master_MatrixOG = zeros(40960,N*NN);

% Prealocating space for pyramids
smFSPyr1_4 = zeros(40960,NN);
smFSPyr1_4OG = zeros(40960,NN); 

%% Extracting info
tic
S = 4; %scale
O = 4; %orientation
N = 128; %size of field - we will stick with 128
N1s = 10;

%Reading all images

for i = 1:length(ImgDir)
    img1 = imread(join([ImgDir3,ImgDir{i}]));
    img2 = hsv2rgb(rgb2hsv(img1).*repmat(reshape([1,0,1],[1 1 3]), [size(img1,1) size(img1,2), 1]));
    img3 = squeeze(img2(:,:,1));
    
    xmax = size(img3,1) - N;
    ymax = size(img3,2) - N;
    
    rectSectx = randi(xmax,N1s,1);
    rectSecty = randi(ymax,N1s,1);
    
    for ii = 1:length(rectSectx)
        img4 = img3(rectSectx(ii):rectSectx(ii)+128,rectSecty(ii):rectSecty(ii)+128);
        [spyrTest, pind, pyr] = getSpyr4(1*img4, S, O);
        
        SpyrMatrix1(:,ii) = spyrTest;
        
    end
    
    SpyrMatrix2(:,(((i*N1s)-9):(i*N1s))) = SpyrMatrix1;
    
end
toc

%Time taken: 7min

%% Saving Function
save('SPYR_580pic_10sects','SpyrMatrix2','-v7.3');



%%
%Load Mat File
S = 4; %scale
O = 4; %orientation
N = 128; %size of field - we will stick with 128
N1s = 10;

tic
load('SPYR_580pic_10sects.mat');
toc
%~1.5 min

%% Display random spyr
x = randi(100,1);

figure(3)
subplot(2,2,1)
imshow(join([ImgDir3,ImgDir{x}]))
subplot(2,2,3)
spyrDisp4(SpyrMatrix2(:,x)/max(abs(SpyrMatrix2(:,x))), pind, pyr);
subplot(2,2,2)
spyrDisp4(SpyrMatrix2(:,x), pind, pyr);

SpM_m = mean(SpyrMatrix2(:,x));
SpM_std = var(SpyrMatrix2(:,x),0,1);
SpM1 = (SpyrMatrix2(:,x)-SpM_m)/SpM_std;

subplot(2,2,4)
spyrDisp4(SpM1/max(SpM1), pind, pyr);



%% Normalize and such 
tic
SpM_m = mean(SpyrMatrix2(:,:));
SpM_std = var(SpyrMatrix2(:,:),0,1);
SpM_max = max(abs(SpyrMatrix2(:,:)));

for i = 1:size(SpyrMatrix2,2)
    SpM2(:,i) = ((SpyrMatrix2(:,i)-SpM_m(i))/SpM_std(i)); %Centered 
    SpM2_abs(:,i) = abs(SpM2(:,i));
    
    SpM1(:,i) = SpM2(:,i)/(max(SpM2(:,i))); %normalized
    SpM1_abs(:,i) = abs(SpM1(:,i));
end
toc
%%
figure(4)
x = randi(100,1);
subplot(1,2,1)
spyrDisp4(SpM1(:,x)/max(SpM1(:,x)), pind, pyr);
subplot(1,2,2)
spyrDisp4(SpM2(:,x)/max(SpM2(:,x)), pind, pyr);


%% Basis Set: Create a set of randomly generated Filters

tic
freqN = 1:6;

%Aperture must remain the same 
BasisSet = input('Pick number of filters: ');
ap = 6.5;

%Create Random orientation
spSet = 2*pi*rand(BasisSet,1);

%Create Random Spatial Frequence
spFreq = (pi*1/6) + ((pi*4/6)*rand(BasisSet,1));

%Create Random phase set
phaseSet = 2*pi*(rand(BasisSet,1));

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