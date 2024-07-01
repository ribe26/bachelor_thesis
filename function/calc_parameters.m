function [PK] = calc_parameters(a,P1,N1,N2,N3)
    %P1,N1,N2,N3から仰角a°のパラメータを計算する
    N=10;
    Fs=48000;
    
    P1_W=P1(N+1,1);
    P1_B=P1(N+1,2);
    P1_G=P1(N+1,3);
    
    N1_W=N1(N+1,1);
    N1_B=N1(N+1,2);
    N1_G=N1(N+1,3);
    
    N2_W=N2(N+1,1);
    N2_B=N2(N+1,2);
    N2_G=N2(N+1,3);
    
    N3_W=N3(N+1,1);
    N3_B=N3(N+1,2);
    N3_G=N3(N+1,3);
    
    
    for i=1:N
        P1_W=P1_W+a^(N+1-i)*P1(i,1);
        P1_B=P1_B+a^(N+1-i)*P1(i,2);
        P1_G=P1_G+a^(N+1-i)*P1(i,3);
    
        N1_W=N1_W+a^(N+1-i)*N1(i,1);
        N1_B=N1_B+a^(N+1-i)*N1(i,2);
        N1_G=N1_G+a^(N+1-i)*N1(i,3);
    
        N2_W=N2_W+a^(N+1-i)*N2(i,1);
        N2_B=N2_B+a^(N+1-i)*N2(i,2);
        N2_G=N2_G+a^(N+1-i)*N2(i,3);
    
    
        N3_W=N3_W+a^(N+1-i)*N3(i,1);
        N3_B=N3_B+a^(N+1-i)*N3(i,2);
        N3_G=N3_G+a^(N+1-i)*N3(i,3);
    
    end
    
    P1_W=P1_W*2*pi/Fs;
    N1_W=N1_W*2*pi/Fs;
    N2_W=N2_W*2*pi/Fs;
    N3_W=N3_W*2*pi/Fs;
    
    PK=[P1_W,P1_B,P1_G;N1_W,N1_B,N1_G;N2_W,N2_B,N2_G;N3_W,N3_B,N3_G];
end