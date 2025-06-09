function h=ATS_h(A,j,H)
%GTSP的逻辑时钟相位调整
    
        deltah=A(j,4)-A(j,5);
    
h=H+deltah;
end