function [ntable,map]=tNormalizeTable(table)
% normalize changes a table dim x term x multiid to a
% dim x term table in which each row contains indices from
% 1 to  the number of unique identification of the operators
% (i.e., each index corresponds to a multiindex of the form
% (r,b,:) of the original table) arranged by decreasing
% frequency.
% By normalizing in this way the same code can then be applied
% to optimizing evaluation of partial derivatives and
%
sz=size(table,1:3);
dim=size(table,1);
ops=cell(dim,1);
id=cell(dim,1);
count=cell(dim,1);
for i=1:dim
    % find unique univariate operators and assign temporary id
    % by position
    [ops{i},~,id{i}]=unique(reshape(table(i,:,:),sz(2),sz(3)),'row');
    % find new ids by counting
    [count{i},perm] = groupcounts(id{i}); assert(all(perm'==1:numel(perm)))
    [count{i},ord]=sort(count{i},'descend');
    ops{i}=reshape(ops{i}(ord,:),1,numel(ord),sz(3));
    iord(ord)=1:numel(ord);%#ok
    id{i} =reshape(iord(id{i}),1,sz(2),1);
end
map=struct("operators",ops,"count",count);
ntable=cat(1,id{:});

otable=arrayfun(@(i)ops{i}(1,ntable(i,:),:),1:dim,'UniformOutput',false);
assert(all(cat(1,otable{:})-table==0,'all'))
end