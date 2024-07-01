
function current =  ptnsrc_updatecurrentpoint(polls,current)
    f_polls = polls.f;
    [~,I] = min(f_polls);
    current.x = polls.xcatch(I,:);
    current.f = polls.f(I);
end