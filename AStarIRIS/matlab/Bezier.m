function P=Bezier(B,C,t)
N=size(B,2); %Order of the Bezier curve
%Bp=fliplr(B);
P=zeros(size(C,1),length(t));
for i=1:N
    P=P+polyval(B(i,:),t).*C(:,i);
end

