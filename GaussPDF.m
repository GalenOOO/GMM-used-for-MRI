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