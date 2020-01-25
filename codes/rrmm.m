%% Written by Muhammet Balcilar, France
% all rights reverved


clc
addpath('docde');


filename = 'canny.jpg';
image = imread(filename);

%% Region Growing method
[Region1, Region2, Region3, Region4]=RegionGrowing(image);

figure; imshow(uint8(image));
hold on;
DrawLine(Region1);
hold off;
title('Scale 1');

figure; imshow(uint8(image));
hold on;
DrawLine(Region2);
hold off;
title('Scale 2');

figure; imshow(uint8(image));
hold on;
DrawLine(Region3);
hold off;
title('Scale 3');

figure; imshow(uint8(image));
hold on;
DrawLine(Region4);
hold off;
title('Scale 4');


%% Region merging method 
% method has two parameters minimum adjacent pixel which means connected
% regions has at least this number of pixel connected
mnadj=10;
% threshold value to make decision to merge two region or nor
RMThresh=3.5;

RegionResult=RegionMerging(image,Region1,mnadj,RMThresh);

% figure results
figure; imshow(uint8(image));
hold on;
DrawLine(RegionResult);
hold off;
title('Final segmentation');






function DrawLine(SegI)
% input: SegI   --- Segments with segment labels
% output: ImgS  --- Segments with line boundaries

[m, n] = size(SegI);
ImgS = zeros(m, n);
% RGB = label2rgb(SegI);
% figure; imshow(RGB);
% hold on;
BRegion = imdilate(SegI, [0 1 0; 1 1 1; 0 1 0]);
Boundary = BRegion & ~SegI;


for i = 1:max(SegI(:))
    S = zeros(m, n);
    [x, y] = find(SegI == i);
    for j = 1:length(x)
        S(x(j), y(j)) = 1;
    end
    [B,L] = bwboundaries(S,'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1);
    end
    ImgS = ImgS + S;
end
end
function [Region1, Region2, Region3, Region4]=RegionGrowing(image_org)

%% Written by Muhammet Balcilar, France
% all rights reverved


[m,n,d] = size(image_org);


X = reshape(double(image_org), m*n,d);
[tmp,M,tmp2,P] = kmeansO(X,[],16,0,0,0,0);
map = reshape(P, m, n);

for w = 1:4
    W = GenerateWindow(w);
    JI{w} = JImage(map, W);    
end


ImgQ = class2Img(map, image_org);
Region = zeros(m, n);

u = mean(JI{4}(:));
s = std(JI{4}(:));
Region = ValleyD(JI{4},  4, u, s); 
Region = ValleyG1(JI{4}, Region);  
Region = ValleyG1(JI{3}, Region);  
Region = ValleyG2(JI{1}, Region);  
Region4 = Region;


w = 3;
Region = SpatialSeg(JI{3}, Region, w);
Region = ValleyG1(JI{2}, Region);
Region = ValleyG2(JI{1}, Region);
Region3 = Region;

w = 2;
Region = SpatialSeg(JI{2}, Region, w);
Region = ValleyG1(JI{1}, Region);
Region = ValleyG2(JI{1}, Region);
Region2 = Region;

w = 1;
Region = SpatialSeg(JI{1}, Region, w);
Region = ValleyG2(JI{1}, Region);
Region1 = Region;
 
end
function RegionResult=RegionMerging(image,Region,mnadj,RMThresh)
%% Written by Muhammet Balcilar, France
% all rights reverved



%% main bloc of algorithm

% find which rregions are connected to each other
[N, Con]=findNeighbour(Region,mnadj);

% find every regions sttistci which means means and covariances
Stat=findStatistic(Region,image);

% calculate initial regions S similarity values
Sval=calcSval(Stat,Con);

% separate R,G, and B channel of image
R=image(:,:,1);
G=image(:,:,2);
B=image(:,:,3);
R=double(R(:));
G=double(G(:));
B=double(B(:));

% while minimum similarity of connected regions are still less than given
% certain threshold continue the process
while min(Sval)<RMThresh
    % find two region which are connected and has minimum S value
    [a, b]=min(Sval);
    % take the region which will be alive
    i=Con(b,1);
    % take the region number which will be dead
    j=Con(b,2);
    % remove this connection and s value since this two region are no
    % longer separate regions they will merge to each other
    Con(b,:)=[];
    Sval(b)=[];
    % make all j th region pixel as ith region
    Region(Region==j)=i;
    % find all pixel which assigned new merged ith region
    I=find(Region==i);
    
    % calculate new merged ith regions statistics
    tmp=[R(I) G(I) B(I)];
    Stat{i}.mean=mean(tmp);
    Stat{i}.cov=cov(tmp);
    
    % make all jth class name as ith since jth class no longer alive
    Con(Con==j)=i;
    % find new merged ith region connections
    I=find(Con(:,1)==i | Con(:,2)==i);
    % calculate new similartiy between the new merged ith class and its
    % neighbourhood
    for itr=1:length(I)
        muA =Stat{Con(I(itr),1)}.mean';
        muB =Stat{Con(I(itr),2)}.mean';
        covA=Stat{Con(I(itr),1)}.cov;
        covB=Stat{Con(I(itr),2)}.cov;        
        Sval(I(itr))=(muA-muB)'*inv(covA+covB)*(muA-muB);
    end
    
end
% rename all region name from 1 to the end ascendend
I=unique(Region(:));
tmp=Region;
for i=1:length(I)
    tmp(Region==I(i))=i;
end
RegionResult=tmp;



end
function Sval=calcSval(Stat,Con)

%% Written by Muhammet Balcilar, France
% all rights reverved


Sval=zeros(size(Con,1),1);

for itr=1:size(Con,1)
    muA =Stat{Con(itr,1)}.mean';    
    muB =Stat{Con(itr,2)}.mean';
    covA=Stat{Con(itr,1)}.cov;    
    covB=Stat{Con(itr,2)}.cov;   
    
    Sval(itr)=(muA-muB)'*inv(covA+covB)*(muA-muB);
end
end
function [N, Con]=findNeighbour(Region,mnadj)

%% Written by Muhammet Balcilar, France
% all rights reverved

n=length(unique(Region(:)));
N=zeros(n,n);
Con=[];

for i=1:size(Region,1)-1
    for j=1:size(Region,2)-1
        tmp=Region(i:i+1,j:j+1);
        I=unique(tmp(:));
        if length(I)>1
            
            if N(I(1),I(2))==mnadj-1
                Con=[Con;min(I(1),I(2)) max(I(1),I(2)) ];
            end
            N(I(1),I(2))=N(I(1),I(2))+1;
            N(I(2),I(1))=N(I(2),I(1))+1;
        end
    end
end
end
function Stat=findStatistic(Region,image)
%% Written by Muhammet Balcilar, France
% all rights reverved
n=length(unique(Region(:)));

R=image(:,:,1);
G=image(:,:,2);
B=image(:,:,3);

R=double(R(:));
G=double(G(:));
B=double(B(:));

region=Region(:);

for i=1:n
    I=find(region==i);
    tmp=[R(I) G(I) B(I)];
    Stat{i}.mean=mean(tmp);
    Stat{i}.cov=cov(tmp);
end
end

