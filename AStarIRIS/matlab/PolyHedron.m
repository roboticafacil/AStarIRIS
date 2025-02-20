classdef PolyHedron < duplicableObject
    properties
        vertices;
        A;
        b;
        centroid;
    end
    methods
        function obj=PolyHedron(vertices,A,b)
            if nargin<1
                obj.vertices=[];
            else
                obj.vertices=vertices;
            end
            if nargin<3
                if (~isempty(obj.vertices))
                    [obj.A,obj.b]=vert2con(obj.vertices);
                else
                    obj.A=[];
                    obj.b=[];
                end
            else
                obj.A=A;
                obj.b=b;
            end
            obj.centroid=[];
        end
        function P=duplicate(obj)
            P=PolyHedron(obj.vertices,obj.A,obj.b);
        end
        function inside=insidePolyhedron(obj,q)
            inside=all(obj.A*q-obj.b<=0);    
        end
        function [ai,bi]=findActiveConstraints(obj,q,tol)
            if (~isempty(obj.b))
                idx=find(abs(obj.A*q-obj.b)<tol);
                ai=obj.A(idx,:);
                bi=obj.b(idx);
            else
                ai=[];
                bi=[];
            end
        end
        function addDimensionToPolyHedron(obj,range,order)
            %range is a list the with limit values of the new dimension and order is a list with the
            %new state ordering
            if isempty(range)
                n=size(range,2);
                obj.A=[obj.A zeros(size(obj.A,1),n)];
                %If there's no range on the new dimensions, then we can not
                %compute vertices
                obj.vertices=[];
                if nargin==3
                    if (~isempty(order))
                        obj.A=obj.A(:,order);
                    end
                end
            else
                Range=PolyHedronRange(range);
                n=size(Range.A,2);
                obj.A=[obj.A zeros(size(obj.A,1),n);zeros(obj.A,2) Range.A];
                obj.b=[obj.b;Range.b];
                if nargin==3
                    obj.A=obj.A(:,order);
                end
                if (~isempty(obj.vertices))
                    obj.vertices=con2vert(obj.A,obj.b);
                end
            end
        end
        function computeCentroid(obj)
            if (isempty(obj.vertices))
                obj.vertices=con2vert(obj.A,obj.b);
            end
            obj.centroid=mean(obj.vertices,2);
        end

        function printPolyHedronEigenFormat(obj,num)
            s=sprintf("Eigen::Matrix<double, %d, %d> A%d({",size(obj.A,1),size(obj.A,2),num);
            if (size(obj.A,1)>0)
                s=sprintf("%s{%f",s,obj.A(1,1));
                for j=2:size(obj.A,2)
                    s=sprintf("%s,%f",s,obj.A(1,j));
                end
                s=sprintf("%s}",s);
                for i=2:size(obj.A,1)
                    s=sprintf("%s,{%f",s,obj.A(i,1));
                    for j=2:size(obj.A,2)
                        s=sprintf("%s,%f",s,obj.A(i,j));
                    end
                    s=sprintf("%s}",s);
                end
            end
            s=sprintf("%s});\n",s);
            s=sprintf("%sEigen::Vector<double, %d> b%d(",s,size(obj.b,1),num);
            if (size(obj.b,1)>0)
                if (size(obj.b,1)==1)
                    s=sprintf("%s%f",s,obj.b(1));
                else
                    s=sprintf("%s{%f",s,obj.b(1));
                    for i=2:size(obj.b,1)
                        s=sprintf("%s,%f",s,obj.b(i));
                    end
                    s=sprintf("%s}",s);
                end
            end
            s=sprintf("%s);",s);
            disp(s);
        end
    end
    methods(Static)
        function p=join(p1,p2)
            A=[p1.A;p2.A];
            b=[p1.b;p2.b];
            p=PolyHedron([],A,b);
            if (~isempty(p1.vertices)&&(~isempty(p2.vertices)))
                p.vertices=con2vert(A,b);
            end
        end

    end
end