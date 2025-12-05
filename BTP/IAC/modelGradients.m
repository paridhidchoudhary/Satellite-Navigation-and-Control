function [loss, gradients] = modelGradients(dlnet, dlT)
    % Enable automatic differentiation with respect to dlnet parameters
    dlX = forward(dlnet, dlT);    % predicted states: [62, batchSize]

    % Compute dX/dt via automatic differentiation w.r.t. time
    dlDX = dlgradient(sum(dlX, 'all'), dlT);

    % Evaluate true dynamics for each time point
    numSamples = size(dlT, 2);
    residuals = zeros(size(dlX), 'like', dlX);
    
    for i = 1:numSamples
        t_val = extractdata(dlT(1,i));
        x_val = extractdata(dlX(:,i));
        
        % Evaluate dynamics
        [x_dot_true, ~, ~] = coupledEOM_smc_only(t_val, x_val);
        residuals(:,i) = dlDX(:,i) - x_dot_true;
    end

    % Compute physics loss (MSE over residuals)
    loss = mean(sum(residuals.^2, 1), 'all');

    % Get gradients w.r.t. the learnable parameters of the network
    gradients = dlgradient(loss, dlnet.Learnables);
end
