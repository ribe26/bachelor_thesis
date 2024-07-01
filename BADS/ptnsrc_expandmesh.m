% This function refines the mesh 

function current =  ptnsrc_expandmesh(current,conf)
    current.length = current.length * conf.expandfactor;
end