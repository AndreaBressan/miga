function C = multinomialMulti(alpha,split)
%LEIBNIZCOEFF  Multi-index multinomial coefficient using sequential nchoosek
%
% alpha : m×1 multi-index (total derivative)
% split : m×n matrix, each column is alpha_i (split among factors)
%
% Returns:
%   C = alpha! / (alpha_1! alpha_2! ... alpha_n!)
%   computed using sequential nchoosek (integer-based)
remaining = alpha;
n = size(split,2);
C = 1;
for k = 1:n
    g = split(:,k);
    C = C * multiNchooseK(remaining,g);
    remaining = remaining - g;
end
end