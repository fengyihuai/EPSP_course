function transV = dwt_Ltrans(vector, l_dwt)
% waiting to perfection!!!
% dwt长度的反演化拟合公式
% 仅适用于level=4
l_filter = 10/2; % sym5

% transV = 2^l_dwt.*vector - 30*l_filter+2^l_dwt - 1 - 16;
transV = vector;
for i = 1:l_dwt
    transV = 2.*ceil((2.*(transV - l_filter)+1)./2);
end