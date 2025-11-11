% 对字典中的原子按照离散差分的绝对值进行排序
%
% 2014-03-21

function Dsort = dictsort(D)
[atompix, atomnum] = size(D);
gradatom = zeros(atomnum, 1);
for i = 1:atomnum
    atomtemp = D(:, i);
    atom2 = reshape(atomtemp, sqrt(atompix), sqrt(atompix));
    gradatom(i) = sum(sum(abs(atom2(1:end-1,:) - atom2(2:end,:)))) + sum(sum(abs(atom2(:,1:end-1) - atom2(:,2:end))));
end

[B, IX] = sort(gradatom);
Dsort = D(:, IX);