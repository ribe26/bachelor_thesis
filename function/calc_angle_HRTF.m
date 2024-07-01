function [generated_HRTF] = calc_angle_HRTF(HRTF_list,ang,phi1,phi2)
    %UNTITLED3 この関数の概要をここに記述
    %   ph1のHRTFとphi2のHRTFから仰角ang°のHRTFを推定
    
    %補間に用いるHRTFのインデックスを求める
    HRTF1_index=phi1/10+1;
    HRTF2_index=phi2/10+1;
    
    %補間に用いるHRTFへの重みづけを計算する
    omomi1=(phi2-ang)/(phi2-phi1);
    omomi2=(ang-phi1)/(phi2-phi1);
    
    
    
    %重みづけに基づいて補間する
    generated_HRTF=[];
    for i=1:length(HRTF_list(:,HRTF1_index))
        angle1=angle(HRTF_list(i,HRTF1_index));
        angle2=angle(HRTF_list(i,HRTF2_index));
        abs1=abs(HRTF_list(i,HRTF1_index));
        abs2=abs(HRTF_list(i,HRTF2_index));
    
        mixed_abs=abs1*omomi1+abs2*omomi2;
        mixed_angle=angle1*omomi1+angle2*omomi2;
        z=mixed_abs*(cos(mixed_angle)+j*sin(mixed_angle));
        generated_HRTF=[generated_HRTF;z];
    end
end