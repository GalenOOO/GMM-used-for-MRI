clc;
clear;
close all;
%% Loading ski_image.jpg
likelihood=0;
j1=(imread('t1_pn0_img.png'));
j2=(imread('t2_pn0_img.png'));
r1 = im2double(j1(:));
r2 = im2double(j2(:));
X=[r1,r2];

%% Initoalizing the parameters for EM algorithm
[muj1,maskj1]=kmeanspic(j1,4);
muj1=muj1/255;
pointNum=size(X,1);
mix=zeros(1,4);
for i =1:4
    mix(i)=sum(maskj1(:)==i)/pointNum;
end

[muj2,maskj2]=kmeanspic(j2,4);
muj2=muj2/255;
mu=[muj1;muj2];

% mu1=[0;0];
% mu2=[0.8220;0.3017];
% mu3=[0.2744;0.8820];
% mu4=[0.5181;0.5733];
% mu=[mu1,mu2,mu3,mu4];
% clear mu1 mu2 mu3 mu4;
% mix = [0.5179;0.3643;0.0509;0.0670];

cov1=eye(2,2);
cov2=eye(2,2);
cov3=eye(2,2);
cov4=eye(2,2);
itr =1;
gamma = zeros(pointNum,4);
likelihood_vect = [];
%% Updating means, weights,covariance and gamma matrix until convergence
while itr<70
    Nk_1=0;newmu1=[0;0];
    Nk_2=0;newmu2=[0;0];
    Nk_3=0;newmu3=[0;0];
    Nk_4=0;newmu4=[0;0];

    prob = normal(X,mu,cov1,cov2,cov3,cov4);
    for i=1:pointNum
        gamma(i,1)=(mix(1)*(prob(i,1)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        gamma(i,2)=(mix(2)*(prob(i,2)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        gamma(i,3)=(mix(3)*(prob(i,3)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        gamma(i,4)=(mix(4)*(prob(i,4)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        Nk_1=Nk_1+gamma(i,1);
        Nk_2=Nk_2+gamma(i,2);
        Nk_3=Nk_3+gamma(i,3);
        Nk_4=Nk_4+gamma(i,4);
        newmu1=newmu1+gamma(i,1)*X(i,:)';
        newmu2=newmu2+gamma(i,2)*X(i,:)';
        newmu3=newmu3+gamma(i,3)*X(i,:)';
        newmu4=newmu4+gamma(i,4)*X(i,:)';
    end
    newmu1=newmu1/Nk_1;
    newmu2=newmu2/Nk_2;
    newmu3=newmu3/Nk_3;
    newmu4=newmu4/Nk_4;   
    cov1=zeros(2,2)+0.00001*eye(2);
    cov2=zeros(2,2)+0.00001*eye(2);
    cov3=zeros(2,2)+0.00001*eye(2);
    cov4=zeros(2,2)+0.00001*eye(2);
    for i=1:pointNum
        cov1=cov1+gamma(i,1)*(X(i,:)'-newmu1)*(X(i,:)-newmu1');
        cov2=cov2+gamma(i,2)*(X(i,:)'-newmu2)*(X(i,:)-newmu2');
        cov3=cov3+gamma(i,3)*(X(i,:)'-newmu3)*(X(i,:)-newmu3');
        cov4=cov4+gamma(i,4)*(X(i,:)'-newmu4)*(X(i,:)-newmu4');
    end
    cov1=cov1/Nk_1;
    cov2=cov2/Nk_2;
    cov3=cov3/Nk_3;
    cov4=cov4/Nk_4;
    itr

    mu=[newmu1,newmu2,newmu3,newmu4];    
    mix(1)=Nk_1/pointNum;
    mix(2)=Nk_2/pointNum;
    mix(3)=Nk_3/pointNum;
    mix(4)=Nk_4/pointNum;
    
    prob=normal(X,mu,cov1,cov2,cov3,cov4);
    error=0;
    old_likelihood = likelihood;
    likelihood=0;
    % Calculating log-likelihood to find error
    for i=1:pointNum
        likelihood = likelihood+log((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
    end
    error=abs(old_likelihood-likelihood)
%     likelihood_vect = [likelihood_vect old_likelihood];
    itr=itr+1;
end
%% Displaying the segmented image
gamma2 = [gamma(:,1),gamma(:,3),gamma(:,4),gamma(:,2)];
[x,index2]=max(gamma2,[],2);
bb = reshape(index2,[217,181]);
bb = bb-1;
imshow(bb,[])

% [a,index]=max(gamma,[],2);
% b=reshape(index,[217,181]);
% imshow(b,[]);

c=imread('gt_img.png');
figure
imshow(c,[]);
cc = double(c);

mcr = sum(sum((bb-cc)~=0,2))/pointNum

% im1 = imread('t1_pn0_img.png');
% im2 = imread('t2_pn0_img.png');
% figure
% subplot(121)
% imshow(im1,[])
% subplot(122)
% imshow(im2,[])


% gamma=gamma*255;
% im = uint8(reshape(gamma,[321,481,3]));
% figure(2);
% imshow(im);
% hold on;
% iter_ax = 1:69;
% figure(3);
% plot(iter_ax, likelihood_vect,'-','LineWidth',3);
% title('Log Likelihood vs Number of Iterations');
% ylabel('Updated Log Likelihood');
% xlabel('Number of Iterations');


function [mu,mask]=kmeanspic(ima,k)
%   ¹ŠÄÜ£ºÔËÓÃk-meansËã·š¶ÔÍŒÏñœøÐÐ·Öžî
%   ÊäÈë: ima-ÊäÈëµÄ»Ò¶ÈÍŒÏñ           k-·ÖÀàÊý
%   Êä³ö:  mu-ŸùÖµÀàÏòÁ¿                mask-·ÖÀàºóµÄÍŒÏñ
ima=double(ima);
copy=ima;
ima=ima(:);
mi=min(ima);
ima=ima-mi+1;
s=length(ima);
% ŒÆËãÍŒÏñ»Ò¶ÈÖ±·œÍŒ
m=max(ima)+1;
h=zeros(1,m);
hc=zeros(1,m);
for i=1:s
    if(ima(i)>0) h(ima(i))=h(ima(i))+1;end;
end
ind=find(h);
hl=length(ind);
% ³õÊŒ»¯ÖÊÐÄ
mu=(1:k)*m/(k+1);
% start process
while(true)
    oldmu=mu;
    % ÏÖÓÐµÄ·ÖÀà
    for i=1:hl
        c=abs(ind(i)-mu);
        cc=find(c==min(c));
        hc(ind(i))=cc(1);
    end
    %ÖØÐÂŒÆËãŸùÖµ
    for i=1:k,
        a=find(hc==i);
        mu(i)=sum(a.*h(a))/sum(h(a));
    end
    if(mu==oldmu) break;end;
end
% calculate mask
s=size(copy);
mask=zeros(s);
for i=1:s(1),
    for j=1:s(2),
        c=abs(copy(i,j)-mu);
        a=find(c==min(c));
        mask(i,j)=a(1);
    end
end
mu=mu+mi-1;
end

function prob=normal(X,mu,cov1,cov2,cov3,cov4)
prob(:,1)=GaussPDF(X',mu(:,1),cov1);
prob(:,2)=GaussPDF(X',mu(:,2),cov2);
prob(:,3)=GaussPDF(X',mu(:,3),cov3);
prob(:,4)=GaussPDF(X',mu(:,4),cov4);
end

function prob = GaussPDF(Data, Mu, Sigma)
%
% ��ݸ�˹�ֲ��������ÿ����ݵĸ����ܶ� Probability Density Function (PDF)
% ���� -----------------------------------------------------------------
%   o Data:  D x N ��N��Dά���
%   o Mu:    D x 1 ��M��Gaussģ�͵����ĳ�ʼֵ
%   o Sigma: M x M ��ÿ��Gaussģ�͵ķ������ÿ����������ǶԽ���
%                                   ��һ����͵�λ����ĳ˻�
% Outputs ----------------------------------------------------------------
%   o prob:  1 x N array representing the probabilities for the
%            N datapoints.    
[dim,N] = size(Data);
Data = Data' - repmat(Mu',N,1);
prob = sum((Data*inv(Sigma)).*Data, 2);
prob = exp(-0.5*prob) / sqrt((2*pi)^dim * (abs(det(Sigma)^2)+realmin));
end