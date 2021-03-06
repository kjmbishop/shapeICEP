function tensors = C1Tensor
%TdTensor Prompts user for C1 tensor information and outputs Cijk and Dijk
%   outputs struct 'tensors' which contains the specified Cijk and Dijk
%   plus an example particle shape stored in shapeB

    % initialize tensors
    tensors.Cijk=zeros(3,3,3);
    tensors.Dijk=zeros(3,3,3);
    
    % example shape
    tensors.shapeB = sparse(zeros(5,5));
    tensors.shapeB(4,5) = 0.4;
    tensors.shapeB(8,9) = 0.4*1i;
    tensors.shapeB(3,4) = 0.4;
	tensors.shapeB(2,1) = -0.3;
    tensors.shapeB(3,1) = 0.3;
    
    %Display form of tensors
    for i = 1:3
        for j = 1:3
            for k = 1:3
                tensorDisp(i,j,k).C = '0';
                tensorDisp(i,j,k).D = '0';
            end
        end
    end
    
    
    
    tensorDisp(1,1,1).D = 'd1';
    tensorDisp(1,2,2).D = 'd2';
	tensorDisp(1,3,3).D = 'd3';
	tensorDisp(2,1,1).D = 'd4';
    tensorDisp(2,2,2).D = 'd5';
    tensorDisp(2,3,3).D = 'd6';
    tensorDisp(3,1,1).D = 'd7';
    tensorDisp(3,2,2).D = 'd8';
    tensorDisp(3,3,3).D = 'd9';
    
    tensorDisp(1,1,2).D = 'd10';
    tensorDisp(1,2,1).D = 'd10';
	
	tensorDisp(1,3,1).D = 'd11';
	tensorDisp(1,1,3).D = 'd11';
    
    tensorDisp(1,2,3).D = 'd12';
    tensorDisp(1,3,2).D = 'd12';
    
	tensorDisp(2,1,2).D = 'd13';
	tensorDisp(2,2,1).D = 'd13';
    
    tensorDisp(2,1,3).D = 'd14';
	tensorDisp(2,3,1).D = 'd14';
    
    tensorDisp(2,2,3).D = 'd15';
	tensorDisp(2,3,2).D = 'd15';
    
    tensorDisp(3,1,2).D = 'd16';
	tensorDisp(3,2,1).D = 'd16';
    
    tensorDisp(3,1,3).D = 'd17';
	tensorDisp(3,3,1).D = 'd17';
    
    tensorDisp(3,2,3).D = 'd18';
	tensorDisp(3,3,2).D = 'd18';
    
	tensorDisp(1,1,1).C = 'c1';
    tensorDisp(1,2,2).C = 'c2';
	tensorDisp(1,3,3).C = 'c3';
	tensorDisp(2,1,1).C = 'c4';
    tensorDisp(2,2,2).C = 'c5';
    tensorDisp(2,3,3).C = 'c6';
    tensorDisp(3,1,1).C = 'c7';
    tensorDisp(3,2,2).C = 'c8';
    tensorDisp(3,3,3).C = 'c9';
    
    tensorDisp(1,1,2).C = 'c10';
    tensorDisp(1,2,1).C = 'c10';
	
	tensorDisp(1,3,1).C = 'c11';
	tensorDisp(1,1,3).C = 'c11';
    
    tensorDisp(1,2,3).C = 'c12';
    tensorDisp(1,3,2).C = 'c12';
    
	tensorDisp(2,1,2).C = 'c13';
	tensorDisp(2,2,1).C = 'c13';
    
    tensorDisp(2,1,3).C = 'c14';
	tensorDisp(2,3,1).C = 'c14';
    
    tensorDisp(2,2,3).C = 'c15';
	tensorDisp(2,3,2).C = 'c15';
    
    tensorDisp(3,1,2).C = 'c16';
	tensorDisp(3,2,1).C = 'c16';
    
    tensorDisp(3,1,3).C = 'c17';
	tensorDisp(3,3,1).C = 'c17';
    
    tensorDisp(3,2,3).C = 'c18';
	tensorDisp(3,3,2).C = 'c18';
            
    %% Ask for Dijk values
    DString = ['Dijk tensor: \n','((',tensorDisp(1,1,1).D,' ',tensorDisp(1,1,2).D,' ',tensorDisp(1,1,3).D,') \t(',...
        tensorDisp(1,2,1).D,' ',tensorDisp(1,2,2).D,' ',tensorDisp(1,2,3).D,') \t(',...
        tensorDisp(1,3,1).D,' ',tensorDisp(1,3,2).D,' ',tensorDisp(1,3,3).D,'))\n',...
        '((',tensorDisp(2,1,1).D,' ',tensorDisp(2,1,2).D,' ',tensorDisp(2,1,3).D,') \t(',...
        tensorDisp(2,2,1).D,' ',tensorDisp(2,2,2).D,' ',tensorDisp(2,2,3).D,') \t(',...
        tensorDisp(2,3,1).D,' ',tensorDisp(2,3,2).D,' ',tensorDisp(2,3,3).D,'))\n',...
        '((',tensorDisp(3,1,1).D,' ',tensorDisp(3,1,2).D,' ',tensorDisp(3,1,3).D,') \t(',...
        tensorDisp(3,2,1).D,' ',tensorDisp(3,2,2).D,' ',tensorDisp(3,2,3).D,') \t(',...
        tensorDisp(3,3,1).D,' ',tensorDisp(3,3,2).D,' ',tensorDisp(3,3,3).D,'))\n'];
    fprintf(DString);
    
    DPrompt = 'Please enter values for [d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18]. \n';
    DArray = input(DPrompt);
    
    tensors.Dijk(1,1,1) = DArray(1);
    tensors.Dijk(1,2,2) = DArray(2);
    tensors.Dijk(1,3,3) = DArray(3);
    tensors.Dijk(2,1,1) = DArray(4);
    tensors.Dijk(2,2,2) = DArray(5);
    tensors.Dijk(2,3,3) = DArray(6);
    tensors.Dijk(3,1,1) = DArray(7);
    tensors.Dijk(3,2,2) = DArray(8);
    tensors.Dijk(3,3,3) = DArray(9);
    
    tensors.Dijk(1,1,2) = DArray(10);
    tensors.Dijk(1,2,1) = DArray(10);
    
    tensors.Dijk(1,1,3) = DArray(11);
    tensors.Dijk(1,3,1) = DArray(11);
    
    tensors.Dijk(1,2,3) = DArray(12);
    tensors.Dijk(1,3,2) = DArray(12);
    
    tensors.Dijk(2,1,2) = DArray(13);
    tensors.Dijk(2,2,1) = DArray(13);
    
    tensors.Dijk(2,1,3) = DArray(14);
    tensors.Dijk(2,3,1) = DArray(14);
    
    tensors.Dijk(2,2,3) = DArray(15);
    tensors.Dijk(2,3,2) = DArray(15);
    
    tensors.Dijk(3,1,2) = DArray(16);
    tensors.Dijk(3,2,1) = DArray(16);
    
    tensors.Dijk(3,1,3) = DArray(17);
    tensors.Dijk(3,3,1) = DArray(17);
    
    tensors.Dijk(3,2,3) = DArray(18);
    tensors.Dijk(3,3,2) = DArray(18);
    
    %% Ask for Cijk values
    CString = ['Cijk tensor: \n','((',tensorDisp(1,1,1).C,' ',tensorDisp(1,1,2).C,' ',tensorDisp(1,1,3).C,') \t(',...
        tensorDisp(1,2,1).C,' ',tensorDisp(1,2,2).C,' ',tensorDisp(1,2,3).C,') \t(',...
        tensorDisp(1,3,1).C,' ',tensorDisp(1,3,2).C,' ',tensorDisp(1,3,3).C,'))\n',...
        '((',tensorDisp(2,1,1).C,' ',tensorDisp(2,1,2).C,' ',tensorDisp(2,1,3).C,') \t(',...
        tensorDisp(2,2,1).C,' ',tensorDisp(2,2,2).C,' ',tensorDisp(2,2,3).C,') \t(',...
        tensorDisp(2,3,1).C,' ',tensorDisp(2,3,2).C,' ',tensorDisp(2,3,3).C,'))\n',...
        '((',tensorDisp(3,1,1).C,' ',tensorDisp(3,1,2).C,' ',tensorDisp(3,1,3).C,') \t(',...
        tensorDisp(3,2,1).C,' ',tensorDisp(3,2,2).C,' ',tensorDisp(3,2,3).C,') \t(',...
        tensorDisp(3,3,1).C,' ',tensorDisp(3,3,2).C,' ',tensorDisp(3,3,3).C,'))\n'];
    fprintf(CString);
    
    CPrompt = 'Please enter values for [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18]. \n';
    CArray = input(CPrompt);
    
    tensors.Cijk(1,1,1) = CArray(1);
    tensors.Cijk(1,2,2) = CArray(2);
    tensors.Cijk(1,3,3) = CArray(3);
    tensors.Cijk(2,1,1) = CArray(4);
    tensors.Cijk(2,2,2) = CArray(5);
    tensors.Cijk(2,3,3) = CArray(6);
    tensors.Cijk(3,1,1) = CArray(7);
    tensors.Cijk(3,2,2) = CArray(8);
    tensors.Cijk(3,3,3) = CArray(9);
    
    tensors.Cijk(1,1,2) = CArray(10);
    tensors.Cijk(1,2,1) = CArray(10);
    
    tensors.Cijk(1,1,3) = CArray(11);
    tensors.Cijk(1,3,1) = CArray(11);
    
    tensors.Cijk(1,2,3) = CArray(12);
    tensors.Cijk(1,3,2) = CArray(12);
    
    tensors.Cijk(2,1,2) = CArray(13);
    tensors.Cijk(2,2,1) = CArray(13);
    
    tensors.Cijk(2,1,3) = CArray(14);
    tensors.Cijk(2,3,1) = CArray(14);
    
    tensors.Cijk(2,2,3) = CArray(15);
    tensors.Cijk(2,3,2) = CArray(15);
    
    tensors.Cijk(3,1,2) = CArray(16);
    tensors.Cijk(3,2,1) = CArray(16);
    
    tensors.Cijk(3,1,3) = CArray(17);
    tensors.Cijk(3,3,1) = CArray(17);
    
    tensors.Cijk(3,2,3) = CArray(18);
    tensors.Cijk(3,3,2) = CArray(18);
		
end

