function [current,polls] = ptnsrc_poll(current, polls, conf)
    x_current = current.x;               % Extract the current design vector X
    objfn = conf.objfn;                  % Extract the objective function
    for idx = 1:polls.n
        current.x_eval(idx,:) = x_current + current.length*polls.patternshape(idx,:);
        polls.xcatch(idx,:) = current.x_eval(idx,:); 
        f = feval(@(x)objfn(x), current.x_eval(idx,:));
        polls.f(idx,1) = reshape(f,1,numel(f));
   
%   current.f = min(polls.f);
    end
    current_f = current.f;
    current_f = repmat(current_f, polls.n, 1);
    polls.success_vector = (current_f > polls.f);
    polls.success = any(polls.success_vector);
   % current.xeval
end