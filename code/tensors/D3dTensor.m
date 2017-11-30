function tensors = D3dTensor
%TdTensor Prompts user for D3d tensor information and outputs Cijk and Dijk
%   outputs struct 'tensors' which contains the specified Cijk and Dijk
%   plus an example particle shape stored in shapeB

    % initialize tensors
    tensors.Cijk=zeros(3,3,3);
    tensors.Dijk=zeros(3,3,3);
    
    % example shape
    tensors.shapeB = sparse(zeros(5,5));
    tensors.shapeB(4,4) = 0.5*1i;
    tensors.shapeB(2,1) = -0.1;
    
    %Display form of tensors
    for i = 1:3
        for j = 1:3
            for k = 1:3
                tensorDisp(i,j,k).C = '0';
                tensorDisp(i,j,k).D = '0';
            end
        end
    end
    
    tensorDisp(1,2,3).D = '-d1';
    tensorDisp(1,3,2).D = '-d1';
    tensorDisp(2,1,3).D = 'd1';
    tensorDisp(2,3,1).D = 'd1';
    
    tensorDisp(1,1,1).D = '-d2';
    tensorDisp(1,2,2).D = 'd2';
    tensorDisp(2,1,2).D = 'd2';
    tensorDisp(2,2,1).D = 'd2';

        
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
    
    DPrompt = 'Please enter values for [d1, d2]. \n';
    DArray = input(DPrompt);
    
    tensors.Dijk(1,2,3) = -DArray(1);
    tensors.Dijk(1,3,2) = -DArray(1);
    tensors.Dijk(2,1,3) = DArray(1);
    tensors.Dijk(2,3,1) = DArray(1);
    
    tensors.Dijk(1,1,1) = -DArray(2);
    tensors.Dijk(1,2,2) = DArray(2);
    tensors.Dijk(2,1,2) = DArray(2);
    tensors.Dijk(2,2,1) = DArray(2);

    
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
    
    fprintf('No independent terms in Cijk!');

end

