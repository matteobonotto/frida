function OUT = vec_linspace(vec_start,vec_end,nstep)

OUT = zeros(size(vec_end,1),nstep);
for ii=1:size(vec_end,1)
    OUT(ii,:) = linspace(vec_start(ii), vec_end(ii), nstep);
end
