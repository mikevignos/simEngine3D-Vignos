%% ME 751 Homework 7 -- Vignos
clear; close all; clc;

%% Problem 3
if 0
    % Run forward Euler for a second of simulation time. Plot the errors in the
    % solution from the analytical solution at each timestep.
    func = @(y,lambda)lambda*y;
    yActual = @(t,lambda)(exp(lambda.*t));
    y0 = 1;
    
    % lambda = -10. This should go unstable at h = 0.2
    lambda = -10;
    timeStart = 0;
    hVec = [0.01 0.1 0.15 0.19 0.2];
    timeEnd = 10;
    hVecLabel = [{'0p01'} {'0p1'} {'0p15'} {'0p19'} {'0p2'}];
    hVecLegend = [{'0.01'} {'0.1'} {'0.15'} {'0.19'} {'0.2'}];
    figure
    hold on
    for iH = 1:length(hVec)
        h = hVec(iH);
        time = timeStart:h:timeEnd;
        
        % Compute actual value for this time step.
        yAct = yActual(time,lambda);
        prob3Results.minus10.(['h' hVecLabel{iH}]).yAct = yAct;
        
        % Perform forward euler for 1 second with these values of lambda and
        % this time step.
        y = zeros(1,length(time));
        for iT = 1:length(time)
            t = time(iT);
            
            % Perform forward euler step
            if (iT == 1)
                y(iT) = y0 + h*func(y0,lambda);
            else
                y(iT) = y(iT-1) + h*func(y(iT-1),lambda);
            end
        end
        
        yFE = [y0 y(1:end-1)];
        
        % Compute error
        prob3Results.minus10.(['h' hVecLabel{iH}]).forwardEulerSol = yFE;
        error = yFE - yAct;
        prob3Results.minus10.(['h' hVecLabel{iH}]).error = error;
        
        % Plot results
        subplot(2,1,1)
        if (iH == 1)
            plot(time,yAct)
        end
        hold on
        plot(time,yFE)
        if (iH == length(hVec))
            xlabel('Time (sec)')
            ylabel('Solution')
            legendEnt = [{'Actual'} hVecLegend];
            legend(legendEnt)
        end
        hold off
        
        subplot(2,1,2)
        hold on
        plot(time,error)
        if (iH == length(hVec))
            xlabel('Time (sec)')
            ylabel('Error in Forward Euler Solution')
            legend(hVecLegend)
            title('Error in Forward Euler Solution for lambda = -10')
        end
        hold off
        
    end
    
    % lambda = -100. This should go unstable at h = 0.02
    lambda = -100;
    hVec = [0.001 0.01 0.015 0.019 0.02];
    timeEnd = 1;
    hVecLabel = [{'0p001'} {'0p01'} {'0p015'} {'0p019'} {'0p02'}];
    hVecLegend = [{'0.001'} {'0.01'} {'0.015'} {'0.019'} {'0.02'}];
    figure
    hold on
    for iH = 1:length(hVec)
        h = hVec(iH);
        time = timeStart:h:timeEnd;
        
        % Compute actual value for this time step.
        yAct = yActual(time,lambda);
        prob3Results.minus10.(['h' hVecLabel{iH}]).yAct = yAct;
        
        % Perform forward euler for 1 second with these values of lambda and
        % this time step.
        y = zeros(1,length(time));
        for iT = 1:length(time)
            t = time(iT);
            
            % Perform forward euler step
            if (iT == 1)
                y(iT) = y0 + h*func(y0,lambda);
            else
                y(iT) = y(iT-1) + h*func(y(iT-1),lambda);
            end
        end
        
        yFE = [y0 y(1:end-1)];
        
        % Compute error
        prob3Results.minus10.(['h' hVecLabel{iH}]).forwardEulerSol = yFE;
        error = yFE - yAct;
        prob3Results.minus10.(['h' hVecLabel{iH}]).error = error;
        
        % Plot results
        subplot(2,1,1)
        if (iH == 1)
            plot(time,yAct)
        end
        hold on
        plot(time,yFE)
        if (iH == length(hVec))
            xlabel('Time (sec)')
            ylabel('Solution')
            legendEnt = [{'Actual'} hVecLegend];
            legend(legendEnt)
        end
        hold off
        
        subplot(2,1,2)
        hold on
        plot(time,error)
        if (iH == length(hVec))
            xlabel('Time (sec)')
            ylabel('Error in Forward Euler Solution')
            legend(hVecLegend)
            title('Error in Forward Euler Solution for lambda = -10')
        end
        hold off
        
    end
    
    % lambda = -1000. This should go unstable at h = 0.002
    lambda = -1000;
    hVec = [0.0001 0.001 0.0015 0.0019 0.002];
    timeEnd = 0.1;
    hVecLabel = [{'0p0001'} {'0p001'} {'0p0015'} {'0p0019'} {'0p002'}];
    hVecLegend = [{'0.0001'} {'0.001'} {'0.0015'} {'0.0019'} {'0.002'}];
    figure
    hold on
    for iH = 1:length(hVec)
        h = hVec(iH);
        time = timeStart:h:timeEnd;
        
        % Compute actual value for this time step.
        yAct = yActual(time,lambda);
        prob3Results.minus10.(['h' hVecLabel{iH}]).yAct = yAct;
        
        % Perform forward euler for 1 second with these values of lambda and
        % this time step.
        y = zeros(1,length(time));
        for iT = 1:length(time)
            t = time(iT);
            
            % Perform forward euler step
            if (iT == 1)
                y(iT) = y0 + h*func(y0,lambda);
            else
                y(iT) = y(iT-1) + h*func(y(iT-1),lambda);
            end
        end
        
        yFE = [y0 y(1:end-1)];
        
        % Compute error
        prob3Results.minus10.(['h' hVecLabel{iH}]).forwardEulerSol = yFE;
        error = yFE - yAct;
        prob3Results.minus10.(['h' hVecLabel{iH}]).error = error;
        
        % Plot results
        subplot(2,1,1)
        if (iH == 1)
            plot(time,yAct)
        end
        hold on
        plot(time,yFE)
        if (iH == length(hVec))
            xlabel('Time (sec)')
            ylabel('Solution')
            legendEnt = [{'Actual'} hVecLegend];
            legend(legendEnt)
        end
        hold off
        
        subplot(2,1,2)
        hold on
        plot(time,error)
        if (iH == length(hVec))
            xlabel('Time (sec)')
            ylabel('Error in Forward Euler Solution')
            legend(hVecLegend)
            title('Error in Forward Euler Solution for lambda = -10')
        end
        hold off
        
    end
    
end

%% Problem 4
if 0
    % Define parameters
    stepSize = 0.01;
    time = 0:stepSize:20;
    tolerance = 10^-8;
    maxIter = 50;
    alphaVec = -2:0.5:2;
    betaVec = -1:0.5:3;
    alphaConst = 1;
    betaConst = 1;
    
    % Define x and y vector and the intial values
    x = zeros(1,length(time));
    y = zeros(1,length(time));
    y(1) = 2;
    x(1) = 0;
    
    % Loop through each time step. At each time step perform Newton-Raphson to
    % compute an approximation of xn and yn.
    % Sensitivity analysis for alpha.
    beta = 1;
    alpha = 0;
    
    for iT = 2:length(time)
        t = time(iT);
        
        % Initialize Newton-Raphson method
        yIter = y(iT - 1);
        xIter = x(iT - 1);
        ynMinus = y(iT - 1);
        xnMinus = x(iT - 1);
        iter = 1;
        while iter < maxIter
            % Compute Jacobian and g for this iteration
            [ g, J ] = prob4Function( xIter, yIter, xnMinus, ynMinus, stepSize, alpha, beta );
            
            % Compute correction factor
            correction = J\-g;
            
            % Apply correction factor
            xIter = xIter + correction(1);
            yIter = yIter + correction(2);
            
            % Check if norm of residual is below the specified tolerance
            if norm(correction) < tolerance
                break;
            end
            iter = iter + 1;
        end
        x(iT) = xIter;
        y(iT) = yIter;
    end
    % Plot x and y
    figure
    hold on
    plot(time,x)
    plot(time,y)
    ylabel('X or Y')
    xlabel('Time (sec)')
    title('alpha = 0 and beta = 1')
    legend('X','Y')
    
    % Sensitivity analysis for alpha.
    beta = betaConst;
    
    for iA = 1:length(alphaVec);
        alpha = alphaVec(iA);
        y(1) = 2;
        x(1) = 0;
        for iT = 2:length(time)
            t = time(iT);
            
            % Initialize Newton-Raphson method
            yIter = y(iT - 1);
            xIter = x(iT - 1);
            ynMinus = y(iT - 1);
            xnMinus = x(iT - 1);
            iter = 1;
            while iter < maxIter
                % Compute Jacobian and g for this iteration
                [ g, J ] = prob4Function( xIter, yIter, xnMinus, ynMinus, stepSize, alpha, beta );
                
                % Compute correction factor
                correction = J\-g;
                
                % Apply correction factor
                xIter = xIter + correction(1);
                yIter = yIter + correction(2);
                
                % Check if norm of residual is below the specified tolerance
                if norm(correction) < tolerance
                    
                    break;
                end
                iter = iter + 1;
            end
            x(iT) = xIter;
            y(iT) = yIter;
        end
        % Plot x and y for changing alpha
        if (iA == 1)
            figure
        end
        subplot(2,1,1)
        hold on
        plot(time,x)
        legendEntry(iA) = cellstr(num2str(alphaVec(iA)));
        if (iA == length(alphaVec))
            xlabel('Time (sec)')
            ylabel('X')
            legend(legendEntry)
            title(['Sensitivity of X and Y to Alpha when Beta = ' num2str(betaConst)])
        end
        hold off
        
        subplot(2,1,2)
        hold on
        plot(time,y)
        if (iA == length(alphaVec))
            xlabel('Time (sec)')
            ylabel('Y')
        end
        hold off
    end
    
    % Sensitivity analysis for beta.
    alpha = alphaConst;
    
    for iB = 1:length(betaVec);
        beta = betaVec(iB);
        y(1) = 2;
        x(1) = 0;
        for iT = 2:length(time)
            t = time(iT);
            
            % Initialize Newton-Raphson method
            yIter = y(iT - 1);
            xIter = x(iT - 1);
            ynMinus = y(iT - 1);
            xnMinus = x(iT - 1);
            iter = 1;
            while iter < maxIter
                % Compute Jacobian and g for this iteration
                [ g, J ] = prob4Function( xIter, yIter, xnMinus, ynMinus, stepSize, alpha, beta );
                
                % Compute correction factor
                correction = J\-g;
                
                % Apply correction factor
                xIter = xIter + correction(1);
                yIter = yIter + correction(2);
                
                % Check if norm of residual is below the specified tolerance
                if norm(correction) < tolerance
                    break;
                end
                iter = iter + 1;
            end
            x(iT) = xIter;
            y(iT) = yIter;
        end
        % Plot x and y for changing alpha
        if (iB == 1)
            figure
        end
        subplot(2,1,1)
        hold on
        plot(time,x)
        legendEntry(iB) = cellstr(num2str(betaVec(iB)));
        if (iB == length(betaVec))
            xlabel('Time (sec)')
            ylabel('X')
            legend(legendEntry)
            title(['Sensitivity of X and Y to Beta, when Alpha = ' num2str(alphaConst)])
        end
        hold off
        
        subplot(2,1,2)
        hold on
        plot(time,y)
        if (iB == length(betaVec))
            xlabel('Time (sec)')
            ylabel('Y')
        end
        hold off
    end
    
end

%% Problem 5
if 1
    % Set various parameters
    y0 = 1;
    t0 = 1;
    tEnd = 10;
    
    % Set necessary functions
    yFunc = @(t)(1./t + 1./t.^2.*tan(1./t + pi - 1));
    gFunc = @(yn, ynMinus, h, t)(yn - ynMinus + h*yn^2 + h/t.^4);
    gJacobianFunc = @(yn, h, t)(1 + 2*h*yn);
    
    % Perform backward Euler with varying step size
    hVec = [0.1 0.05 0.01 0.001 0.0001 0.00001];
    error = zeros(1,length(hVec));
    for iH = 1:length(hVec)
        h = hVec(iH);
        
        % Compute actual solution
        time = t0:h:tEnd;
        yActual = yFunc(time);
        yBE = backwardEuler( gFunc, gJacobianFunc, y0, t0, tEnd, h);
        
        % Compute error at the final time step for this value of h
        error(iH) = abs(yActual(end) - yBE(end));
    end
    
    % Create convergence plot
    figure
    plot(log(hVec),log(error),'bs-')
    xlabel('log(h)')
    ylabel('log(error) @ t = 10')
    title('Convergence Plot for Backward Euler Method')
    
    
end