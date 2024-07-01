function [s] = save_path_elevs(people,azim,LR,elev)
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述
    s=sprintf("graph/elevs/%s/azim%d/%s/%d.png",people,azim,LR,elev);
end