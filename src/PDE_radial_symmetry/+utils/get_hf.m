function hf = get_hf(x,y)
hv = min(y)+(max(y)-min(y))./2;
idx = -1;
for i=1:length(y)
   if y(i) <= hv
       idx = i;
       break;
   end
end
if idx == 1
   hf = x(1); 
end

if idx == -1
   hf = x(end); 
end


if idx > 1
   idx_1 = idx-1;
   dy = y(idx)-y(idx_1);
   dx = x(idx)-x(idx_1);
   hf = x(idx_1) + dx*(-y(idx_1)+hv)/dy;
end
end

