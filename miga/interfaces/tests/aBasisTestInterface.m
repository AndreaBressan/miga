classdef aBasisTestInterface < aMigaTestCase
    properties
        basis
    end

    methods(Test,TestTags=["basis"])

        function verifyDomainValidity(testCase)
            domain = testCase.basis.domain;
            dim=size(domain,1);
            testCase.verifySize(domain,[dim,2]);
            testCase.verifyLessThanOrEqual(domain(:,1),domain(:,2),'domain is empty' );
        end

        function verifyEvaluateInterface(testCase)
            bas=testCase.basis;
            pts=sum(bas.domain,2)/2;
            [v,a]=bas.evaluate(pts);
            testCase.verifyClass(v,"double");
            testCase.verifyInstanceOf(a,"aSparsity");
            testCase.verifyEqual(a.shape, [1, testCase.basis.dimension]);
            [v2,a2]=bas.evaluateCombination(eye(bas.dimension),pts);
            [v3,a3]=bas.evaluateCombination(sparse(eye(bas.dimension)),pts);
            testCase.verifyEqual(a.toFull(v),a2.toFull(v2));
            testCase.verifyEqual(a.toFull(v),a3.toFull(v3));
        end

        function verifyBasisSupportValidity(testCase)
            % Verifies:
            % - Support is valid (left <= right)
            % - Support is within domain
            % - Number of supports equals number of basis functions
            sups = testCase.basis.supports();
            domain = testCase.basis.domain;
            dim=size(domain,1);
            nFun=testCase.basis.dimension;
            testCase.verifySize(sups,[dim,2,nFun],'supports have the wrong shape');

            testCase.verifyTrue(all(sups(:,2) >= sups(:,1)), ...
                'Invalid support intervals (left > right)');
            
            testCase.verifyTrue(...
                all( testCase.basis.periodic | sups(:,1,:) >= domain(:,1) ,'all') && ...
                all( testCase.basis.periodic | sups(:,2,:) <= domain(:,2) ,'all'), ...
                'Support exceeds domain bounds');
        end

        function testOutOfDomainHandling(testCase)
            domain = testCase.basis.domain;
            width=domain(:,2)-domain(:,1);
            per = testCase.basis.periodic;
            pts_outside=domain+[-width,width]/3;
            pts_inside =pts_outside+[width,-width];
            dim = size(domain,1);
            vin = testCase.basis.evaluate(pts_inside);
            for d=1:dim
                pts=pts_inside;
                pts(d,:)=pts_outside(d,:);
                v = testCase.basis.evaluate(pts);
                if ~per(d)
                    testCase.verifyEqual(v, zeros([ testCase.basis.rangeSize,2, testCase.basis.dimension]), ...
                    'Values outside domain should be zero');
                % With outOfDomainError=true, should raise error
                    testCase.basis.outOfDomainError = true;
                    testCase.verifyError(...
                        @() testCase.basis.evaluate(pts_outside), ...
                        'MIGA:POINTS_OUTSIDE_DOMAIN');
                    testCase.basis.outOfDomainError = false;
                else
                    testCase.verifyEqual(v, vin, ...
                    'Values of matching points should be equal');
                end
            end
   
        end

        function verifyBasisSupportValues(testCase)
            % Verifies:
            % - functions are evaluated to zero outside their supports
            dim = size(testCase.basis.domain,1);
            num = testCase.basis.dimension;
            sup = testCase.basis.supports;
            dom = testCase.basis.domain;
            per = testCase.basis.periodic;
            pt  = cat(2,...
                sup(:,1,:)-2*eps*abs(sup(:,1,:)),...
                (sup(:,1,:)+sup(:,2,:))/2,...
                sup(:,2,:)+2*eps*abs(sup(:,2,:)) );
            pt(~per,:,:) = max(dom(~per,1),min(dom(~per,2),pt(~per,:,:)));
            possible=cell(dim,1);
            for d=1:dim
                possible{d}=d+(0:2)*dim;
            end
            ind=cell(dim,1);
            [ind{:}]=ndgrid(possible{:});
            ind=cellfun(@(m) reshape(m,1,3^dim),ind,"UniformOutput",false);
            indices=reshape(cat(1,ind{:}),dim,3^dim);
            shift=dim*3;
            mid=ceil(3^dim/2);
            reference=zeros(prod(testCase.basis.rangeSize),3^dim);
            for b=1:num
                pts=pt(indices+(b-1)*shift);
                val = reshape(testCase.basis.evaluate(pts),prod(testCase.basis.rangeSize),3^dim,num);
                val=val(:,:,b);
                exclude=any(pts<=dom(:,1) | pts>=dom(:,2),1); % TODO handle periodic better
                val(:,mid)=0;
                val(:,exclude)=0; % TODO become smarter
                testCase.assertEqual(val,reference);
            end
        end

        function testPlotFunction(testCase)
            import matlab.unittest.diagnostics.Diagnostic;
            import matlab.unittest.diagnostics.FigureDiagnostic;
            f=testCase.invisibleFigure();
            if prod(testCase.basis.rangeSize)>1 || size(testCase.basis.domain,1)>2
                return;
            end
            testCase.verifyWarningFree(@()testCase.basis.plot);
            testCase.log(3,FigureDiagnostic(f));
        end

    end
end