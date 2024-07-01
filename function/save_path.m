function [s] = save_path(people,elev,azim,LR,status)
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述
    s=sprintf("grapth/%s/elev%d/azim%d/%s_%s.png",people,elev,azim,LR,status);
end