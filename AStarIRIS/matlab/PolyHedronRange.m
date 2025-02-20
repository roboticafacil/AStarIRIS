classdef PolyHedronRange < PolyHedron
    properties
        range;
    end
    methods
        function obj=PolyHedronRange(range)
            obj.range=range;
            n=size(range,2);
            m=numel(range);
            obj.A=zeros(m,n);
            obj.b=zeros(n,1);
            for i=1:2:m
                ii=round((i-1)/2);
                obj.A(i,:)=[zeros(1,ii) -1 zeros(1,n-ii-1)];
                obj.A(i+1,:)=[zeros(1,ii) 1 zeros(1,n-ii-1)];
                obj.b(i)=-range(i);
                obj.b(i+1)=range(i+1);
            end
        end
        function Range=duplicate(obj)
            Range=PolyHedronRange(obj.range);
        end
    end
end