clc;
clear;
%% Loading ski_image.jpg
likelihood=0;
j=im2double(imread('t1_pn0_img.png'));
r=j(:);
X=r;
%% Initoalizing the parameters for EM algorithm
mu=[0,80/255,160/255,1];
cov1=eye(1);
cov2=eye(1);
cov3=eye(1);
cov4=eye(1);
itr =1;
mix = [1/4;1/4;1/4;1/4];
pointSum=size(X,1);
gamma = zeros(pointSum,4);
likelihood_vect = [];
%% Updating means, weights,covariance and gamma matrix until convergence
while itr<70
    Nk_1=0;newmu1=0;
    Nk_2=0;newmu2=0;
    Nk_3=0;newmu3=0;
    Nk_4=0;newmu4=0;

    prob = normal(X,mu,cov1,cov2,cov3,cov4);
    for i=1:pointSum
        gamma(i,1)=(mix(1)*(prob(i,1)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        gamma(i,2)=(mix(2)*(prob(i,2)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        gamma(i,3)=(mix(3)*(prob(i,3)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        gamma(i,4)=(mix(4)*(prob(i,4)))/((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
        Nk_1=Nk_1+gamma(i,1);
        Nk_2=Nk_2+gamma(i,2);
        Nk_3=Nk_3+gamma(i,3);
        Nk_4=Nk_4+gamma(i,4);
        newmu1=newmu1+gamma(i,1)*X(i,1);
        newmu2=newmu2+gamma(i,2)*X(i,1);
        newmu3=newmu3+gamma(i,3)*X(i,1);
        newmu4=newmu4+gamma(i,4)*X(i,1);
    end
    newmu1=newmu1/Nk_1;
    newmu2=newmu2/Nk_2;
    newmu3=newmu3/Nk_3;
    newmu4=newmu4/Nk_4;
    cov1=zeros(1,1)+0.00001*eye(1);
    cov2=zeros(1,1)+0.00001*eye(1);
    cov3=zeros(1,1)+0.00001*eye(1);
    cov4=zeros(1,1)+0.00001*eye(1);
    for i=1:pointSum
        cov1=cov1+gamma(i,1)*(X(i,1)-newmu1)*(X(i,1)-newmu1)';
        cov2=cov2+gamma(i,2)*(X(i,1)-newmu2)*(X(i,1)-newmu2)';
        cov3=cov3+gamma(i,3)*(X(i,1)-newmu3)*(X(i,1)-newmu3)';
        cov4=cov4+gamma(i,4)*(X(i,1)-newmu4)*(X(i,1)-newmu4)';
    end
    cov1=cov1/Nk_1;
    cov2=cov2/Nk_2;
    cov3=cov3/Nk_3;
    cov4=cov4/Nk_4;
    itr
    mu=[newmu1,newmu2,newmu3,newmu4];    
    mix(1)=Nk_1/pointSum;
    mix(2)=Nk_2/pointSum;
    mix(3)=Nk_3/pointSum;
    mix(4)=Nk_4/pointSum;
    
    prob=normal(X,mu,cov1,cov2,cov3,cov4);
    error=0;
    old_likelihood = likelihood;
    likelihood=0;
    % Calculating log-likelihood to find error
    for i=1:pointSum
        likelihood = likelihood+log((mix(1)*prob(i,1))+(mix(2)*prob(i,2))+(mix(3)*prob(i,3))+(mix(4)*prob(i,4)));
    end
    error=abs(old_likelihood-likelihood);
    likelihood_vect = [likelihood_vect old_likelihood];
    itr=itr+1;
end

