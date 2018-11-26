function x = iff(is,b,c)

if isscalar(is)
   if is, x = b; else x = c; end
else
   x = zeros(size(is));
   x(is) = b(is);
   x(~is) = c(~is);
end
