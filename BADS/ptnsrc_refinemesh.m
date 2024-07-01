% This function refines the mesh 

function current =  ptnsrc_refinemesh(current,conf)
    current.length = current.length * conf.refinefactor;
end