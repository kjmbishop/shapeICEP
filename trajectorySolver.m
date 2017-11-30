%% trajectorySolver.m
%
% Prompts user for particle symmetry and nonzero tensor terms. Then,
% generates a still image of the trajectory of the particle.


%% Paths
addpath('code');
addpath('code/tensors');

%% Prompt user for particle symmetry

inputValid = 0;
while ~inputValid
    symmetryPrompt = ['Please specify particle symmetry.\n',...
        'Choices: Td, Dnh, Dn, Cnv, D2d, D3h, D2h, D2, C2v, D3d, D3, C3v, Cnh,\n',...
        'Cn, C3h, S4, S6, C3, C2h, C2, Cs, Ci, C1\n'];
    symmetryChoice = input(symmetryPrompt,'s');

    %% Prompt user for tensor input based on particle symmetry

    if strcmp(symmetryChoice,'Td')
        tensors = TdTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'Dnh')
        tensors = DnhTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'Dn')
        tensors = DnTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'Cnv')
        tensors = CnvTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'D2d')
        tensors = D2dTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'D3h')
        tensors = D3hTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'D2h')
        tensors = D2hTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'D2')
        tensors = D2Tensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C2v')
        tensors = C2vTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'D3d')
        tensors = D3dTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'D3')
        tensors = D3Tensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C3v')
        tensors = C3vTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'Cnh')
        tensors = CnhTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'Cn')
        tensors = CnTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C3h')
        tensors = C3hTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'S4')
        tensors = S4Tensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'S6')
        tensors = S6Tensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C3')
        tensors = C3Tensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C2h')
        tensors = C2hTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C2')
        tensors = C2Tensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'Cs')
        tensors = CsTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,' Ci')
        tensors = CiTensor();
        inputValid = 1;
    elseif strcmp(symmetryChoice,'C1')
        tensors = C1Tensor();
        inputValid = 1;
    else
        fprintf('Not a valid choice! Please try again! \n\n');
    end
end
%% Solve trajectory
isFinished = 0;
while ~isFinished
    %% Request input parameters
    timePrompt = 'How long would you like to run the simulation?\n';
    eulerPrompt = 'Please give the starting orientation in format [phi, theta, psi]\n';

    inputValid = 0;

    timeInput = input(timePrompt);
    timeProduce = timeInput(1);

    while ~inputValid
        eulerInput = input(eulerPrompt);
        sizeEuler = size(eulerInput);
        if (sizeEuler(1) == 1) && (sizeEuler(2) == 3)
            inputValid = 1;
        else
            fprintf('Not a valid set of Euler angles! Please try again! \n\n');
        end
    end

    eulerAngle = eulerInput;

    %% Compute trajectory
    sol = trajectory(tensors.Cijk, tensors.Dijk, eulerAngle, timeProduce);

    %% Visualize
    tensors.rootPath = pwd();
    geometry = shapeGeometry(tensors);
    trajectoryVisualize(sol, geometry);
    
    %% Ask if user wishes to change angle and runtime
    finishPrompt = 'Would you like to change the run time and starting angle? Y/N [Y]\n';
    finishInput = input(finishPrompt,'s');
    if strcmp(finishInput,'N')
        isFinished = 1;
    end
    
end

