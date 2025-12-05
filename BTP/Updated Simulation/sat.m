function sat=sat(s)
epsilon=0.01;
if(s<0)
    s=-s;
end
if(s<epsilon)
    sat=s/epsilon;
else
    sat=sign(s);
end

end