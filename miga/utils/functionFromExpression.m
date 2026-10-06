function fh=functionFromExpression(exprs,vars)
%FUNCTIONFROMEXPRESSION wraps matlabFunction taking into account constant
%values and avoiding errors due to reshape or strange behavior on the
%number of outputs
    [vs,isC]=arrayfun(@isConstExpression,exprs);
    fh=cell(size(exprs));
    fhc=arrayfun(@(v)constantFunction(v,vars),vs(isC),'UniformOutput',false);
    fhn=arrayfun(@(expr)matlabFunction(expr,'Vars',vars),exprs(~isC),'UniformOutput',false);
    [fh{isC}]=fhc{:};
    [fh{~isC}]=fhn{:};
end

function [v,isConst]=isConstExpression(expr)
    try v=double(expr); isConst=true;
    catch, v=nan; isConst=false;
    end
end

function fh=constantFunction(c,vars) %#ok
    p1="fh= @(";
    p2=") c*ones(size(x1));";
    if numel(vars)>1; varString=sprintf('%s, ',vars(1:end-1)); else,varString=""; end
    eval( strcat(p1,varString,vars(end),p2));
end