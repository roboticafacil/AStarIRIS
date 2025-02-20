function B=BernsteinPolynomials(N,a,b)
    syms tau real;
    B=zeros(N,N);
    for n=0:(N-1)
        B(n+1,:)=coeffs(expand(nchoosek(N-1,n)*(((tau-a)/(b-a))^n)*(((b-tau)/(b-a))^(N-1-n))),'All');
    end
end