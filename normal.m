function prob=normal(X,mu,cov1,cov2,cov3,cov4)
prob(:,1)=GaussPDF(X',mu(:,1),cov1);
prob(:,2)=GaussPDF(X',mu(:,2),cov2);
prob(:,3)=GaussPDF(X',mu(:,3),cov3);
prob(:,4)=GaussPDF(X',mu(:,4),cov4);