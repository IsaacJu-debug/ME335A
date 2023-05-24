function [pel]=quadtreeRec(xp, LV, X, i)
    x1 = LV(1,i);x2 = LV(2,i);x3 = LV(3,i);
%     X(2,x2), i
    A2 = ((X(2,x2)-X(2,x3))*(X(1,x1)-X(1,x2))+(X(1,x3)-X(1,x2))*(X(2,x1)-X(2,x2)));

    th1 = ((X(2,x2)-X(2,x3))*(xp(1)-X(1,x2))+(X(1,x3)-X(1,x2))*(xp(2)-X(2,x2))...
        )/A2;

    th2 = ((X(2,x3)-X(2,x1))*(xp(1)-X(1,x3))+(X(1,x1)-X(1,x3))*(xp(2)-X(2,x3))...
        )/A2;

    th3 = 1-th1-th2;
%     [th1, th2, th3]
    if min([th1, th2, th3])<0
        if th1<0
            [row1,col1] = find(LV==x2) ;[row2,col2] = find(LV==x3); 
            [row3,col3] = find(LV==x1) ;
        elseif th2<0
            [row1,col1] = find(LV==x3) ;[row2,col2] = find(LV==x1); 
            [row3,col3] = find(LV==x2) ;
        else
            [row1,col1] = find(LV==x1) ;[row2,col2] = find(LV==x2); 
            [row3,col3] = find(LV==x3) ;
        end
        ic = setdiff(intersect(col1, col2), i);
        if isempty(ic)==1
            pel=0;
        else
            pel = quadtreeRec(xp, LV, X, ic);
        end
%         if isempty(ic)==1
%             ic1 = setdiff(intersect(col3, col1), i);
%             ic2 = setdiff(intersect(col2, col3), i);
%             if ic1==track
%                 pel = quadtreeRec(xp, LV, X, ic2, i);
%             else
%                 pel = quadtreeRec(xp, LV, X, ic1, i);
%             end
%         else
%             pel = quadtreeRec(xp, LV, X, ic, i);
%         end
    else
        pel = i;
    end
end