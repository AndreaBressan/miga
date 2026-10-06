classdef aSparsityTestInterface < aMigaTestCase
    properties
        sparsity
    end


    methods(Test,TestTags=["sparsity"])

        function testToFullAndBackColumnData(testCase)
            s=testCase.sparsity;
            data=rand(s.nnz,1);

            M=s.toFull(data);
            testCase.verifySize(M,s.shape);
            data2=s.toCompressed(M,1);
            testCase.verifyEqual(data,data2);
        end

        function testToFullAndBackRowData(testCase)
            s=testCase.sparsity;
            data=rand(1,s.nnz);

            M=s.toFull(data);
            testCase.verifySize(M,[1,s.shape]);
            data2=s.toCompressed(M,2);
            testCase.verifyEqual(data,data2);
        end

        function testToFullAndBackTensorData(testCase)
            s=testCase.sparsity;
            preSize=[3,5];
            postSize=5;
            data=rand([preSize,s.nnz,postSize]);

            M=s.toFull(data,numel(preSize)+1);
            testCase.verifySize(M,[preSize,s.shape,postSize]);
            data2=s.toCompressed(M,numel(preSize)+1);
            testCase.verifyEqual(data,data2);
        end

        function testToSparseAndBackColumnData(testCase)
            s=testCase.sparsity;
            data=rand(s.nnz,1);

            M=s.toSparse(data);
            testCase.verifySize(M,s.shape);
            data2=s.toCompressed(full(M),1);
            testCase.verifyEqual(data,data2);
        end

        function testToSparseAndBackRowData(testCase)
            s=testCase.sparsity;
            data=rand(1,s.nnz);
            M=s.toSparse(data);
            testCase.verifySize(M,s.shape);
            data2=s.toCompressed(full(M),1);
            testCase.verifyEqual(data,data2');
        end

        function testToSparseAndBackTensorData(testCase)
            s=testCase.sparsity;
            preSize=[3,5];
            data=rand([preSize,s.nnz]);

            M=s.toSparse(data,numel(preSize)+1);

            testCase.verifySize(M,[prod(preSize)*s.shape(1),s.shape(2)]);
            data2=s.toCompressed(reshape(full(M),[preSize,s.shape]),numel(preSize)+1);
            testCase.verifyEqual(data,data2);
        end

        function testNNZCol(testCase)
            s=testCase.sparsity;
            nnz=s.nnzCol;
            % compare partial with total and 
            odd=1:2:s.shape(2);
            even=2:2:s.shape(2);
            nnzO=s.nnzCol(odd);
            nnzE=s.nnzCol(even);
            testCase.verifyEqual(nnzO,nnz(odd));
            testCase.verifyEqual(nnzE,nnz(even));
            % compare total with matrix
            testCase.verifyEqual(nnz,sum(s.toFull(ones(s.nnz,1)),1)');
        end

       function testNNZRow(testCase)
            s=testCase.sparsity;
            nnz=s.nnzRow;
            % compare partial with total and 
            odd=1:2:s.shape(1);
            even=2:2:s.shape(1);
            nnzO=s.nnzRow(odd);
            nnzE=s.nnzRow(even);
            testCase.verifyEqual(nnzO,nnz(odd));
            testCase.verifyEqual(nnzE,nnz(even));
            % compare total with matrix
            testCase.verifyEqual(nnz,sum(s.toFull(ones(s.nnz,1)),2));
        end

        function testPosAndNNZPos(testCase)
            s=testCase.sparsity;
            M=s.toSparse(1:s.nnz);
            [r_p,c_p]=s.pos;
            % check partial results are coherent with all results
            odd=1:2:s.nnz;
            even=2:2:s.nnz;
            [r_po,c_po]=s.pos(odd);
            [r_pe,c_pe]=s.pos(even);
            testCase.verifyEqual(r_po,r_p(odd));
            testCase.verifyEqual(c_po,c_p(odd));
            testCase.verifyEqual(r_pe,r_p(even));
            testCase.verifyEqual(c_pe,c_p(even));
            % check coherence with constructed matrix
            [r_f,c_f,v_f]=find(M);
            [~,perm]=sort(v_f);
            testCase.verifyEqual(r_p,r_f(perm),'wrong row');
            testCase.verifyEqual(c_p,c_f(perm),'wrong col');
            testCase.verifyEqual(s.nzPos(r_f,c_f),v_f,'wrong nnz pos');
        end
    end
end