clc% Clear screen
clear% Clear all variables in memory

% Particle Swarm Topology Optimization
%         Matt Ireland 2019

%% Define Population 

% Generate or Load a population
gen = 0;% 0 to load, 1 to gen

% Load a population
fid = '20200124-1510_pop';% name of the saved .mat population

if gen == 1
    % Generate a population
    count = 1000;% the population size
    
    % Input Data
    % General parameters; Design space
    nndx = 7; % Number of nodes in x
    nndy = 3; % Number of nodes in y
    scale = 1; % Scaling factor for element size
    
    % Boundary conditions
    % Choose BC type
    type = 1;% 1 for MBB
    % Set point load
    pload = 1000;% unitless
    
    % Element mechanical properties
    % Modulus, section area, and MOI
    E = 100000; A = 1; I = pi/4;
    
    % Process interpretation
    % Scan direction
    scanDir = [0 1];% [x y]
    scan = 0;% logical yes/no
    maxRed = 0.1;% maximum percent reduction
    
    % File saving
    fis = 'pop';% filename
    fis = [char(datetime('now','Format','yyyyMMdd''-''HHmm')),'_',fis];% append date and time
end

if gen == 0
    % Process interpretation
    % Scan direction
    scanDir = [0 1];% [x y]
    scan = 0;% logical yes/no
    maxRed = 0.1;% maximum percent reduction
    
    % Load file
    load(fid)
    count = p.count;
   % Change filename
   fis = fid;
   % Check for scanning
   if (scan)
       % Append suffix to filename
       if scanDir(1) == 1
           fis = [fis '_scanX']; 
       else
           if scanDir(2) == 1
           fis = [fis '_scanY']; 
           end
       end
       % Apply process interpreter flag and direction
       for ii = 1:count
           p.popMember(ii).scan = 1;
           p.popMember(ii).scanDir = scanDir;
           p.popMember(ii).maxRed = maxRed;
       end
   end
end

%% Optimization Parameters

% Convergence Criteria
ctol = 1e-3;

% Number of iterations to keep for convergence calc
rollKeep = 3;

% Craziness interval
crazyIter = 25;

% Iteration limit
iterLimit = 2000;

% Velocity Update Tuning
omega = -0.2134;
pPhi = -0.3344;
gPhi = 2.3259;

%% Rendering
% Render best particle, fitness, and velocity values
pRender = 1;% set to zero for none

% Initialize structures for animations
fitMOV(iterLimit) = struct('cdata',[],'colormap',[]);
parMOV(iterLimit) = struct('cdata',[],'colormap',[]);
velMOV(iterLimit) = struct('cdata',[],'colormap',[]);
vebMOV(iterLimit) = struct('cdata',[],'colormap',[]);

% Change figure background color to white
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]);

%% Initialization

if gen == 1
    % Generate Population
    p =  population(count,nndx,nndy,scale,E,A,I,scanDir,scan,maxRed,type,pload);
    % Save population to file
    save(fis,'-v7.3','p','iterLimit','omega','pPhi','gPhi');
end


% Evaluate population fitness
for ii = 1:count
    % Call the fitnessEval method
    p.popMember(ii) = fitnessEval(p.popMember(ii));
end

% Set pBest to current
for ii = 1:count
    % Set pBestPos
    p.popMember(ii).pBestPos = p.popMember(ii).dVar;
    % Set pBestFit
    p.popMember(ii).pBestFit = p.popMember(ii).fitnessVal;
end

% Find gBest particle
bestPar = find([p.popMember.fitnessVal] == min([p.popMember.fitnessVal]));
% Find worst particle
worstPar = find([p.popMember.fitnessVal] == max([p.popMember.fitnessVal]));

% Set gBest of population
for ii = 1:count
    % Set gBestFit
    p.popMember(ii).gBestFit = p.popMember(bestPar).fitnessVal;
    % Set gBestPos
    p.popMember(ii).gBestPos = p.popMember(bestPar).dVar;
end

% Initialize velocities
for ii = 1:count
    p.popMember(ii).vel = datasample(p.popMember(ii).velPot,size(p.popMember(ii).dVar,1))';
end

% Initialize a fitness history vector
popFitHist = ones(iterLimit+1,4);
% Initialize a rolling fitness history vector
popFitRolling = ones(rollKeep,1);
% Find population fitness mean
popFitMean = mean([p.popMember.fitnessVal]);
% Append population fitness to history vector
popFitHist(1,:) = [p.popMember(bestPar).fitnessVal p.popMember(bestPar).fitnessVal...
    popFitMean p.popMember(worstPar).fitnessVal];
% Initialize iteration counter
iter = 1;
% Initialize fitness change checksum
popFitChange = 1;
% Initialize a scattering count
nScatter = 0;

%% Iteration

while iter <= iterLimit
    
    % Check if current interval needs scattering
    if mod(iter,crazyIter) == 0
        % Update the number of scatterings
        nScatter = nScatter + 1;
        % Inspect a rolling fitness vector for convergence
        if nScatter <= rollKeep
            % Write into rolling vector
            popFitRolling(nScatter) = popFitCurrentBest;
        else
            % Trim vector
            popFitRolling = [popFitRolling(2:end); popFitCurrentBest];
            % Find fitness change
            popFitChange = range(popFitRolling);
            % Check for convergence
            if popFitChange < ctol
                % Trim the fitness history vector
                popFitHist = popFitHist(1:iter+1,:);
                % Exit iteration
                break
            end
        end
        % Graph the rolling pre-scatter fitness
        if (pRender)
            fig3 = figure(3);
            cla
            hold on
            plot(popFitRolling);
            drawnow limitrate
        end
        % Random particle positions
        for ii = 1:count
            p.popMember(ii) = randomPosition(p.popMember(ii));
            % Generate new element distribution
            p.popMember(ii) = changeVar(p.popMember(ii),p.popMember(ii).dVarProposed);
            % Only evaluate valid movements
            if p.popMember(ii).validConnec == 1 && p.popMember(ii).validIsland == 1
                % Assign proposed values
                p.popMember(ii).dVar(:,2) = p.popMember(ii).dVarProposed;
                p.popMember(ii).elDist = p.popMember(ii).elDistProposed;
                p.popMember(ii).nodalCoords = p.popMember(ii).nodalCoordsProposed;
                % Update position dependent properties
                p.popMember(ii) = newPositionProps(p.popMember(ii));
                % Call the fitnessEval method
                p.popMember(ii) = fitnessEval(p.popMember(ii));
            else
                % Clear proposed variables
                clear p.popMember(ii).dVarProposed
                clear p.popMember(ii).elDistProposed
                clear p.popMember(ii).nodalCoordsProposed
                % Update the fitness components for this iteration
                p.popMember(ii).fitnessValComponents = ...
                    [p.popMember(ii).fitnessValComponents; p.popMember(ii).fitnessValComponents(end,:)];
            end
        end
        
        % Update pBest and gBest
        % Find gBestFit and gWorstFit
        bestPar = find([p.popMember.fitnessVal] == min([p.popMember.fitnessVal]));
        worstPar = find([p.popMember.fitnessVal] == max([p.popMember.fitnessVal]));
        
        % Update particle gBest and pBest
        for ii = 1:count
            % Check for new gBest
            if p.popMember(bestPar(1)).fitnessVal < p.popMember(ii).gBestFit
                % Set gBestFit
                p.popMember(ii).gBestFit = p.popMember(bestPar(1)).fitnessVal;
                % Set gBestPos
                p.popMember(ii).gBestPos = p.popMember(bestPar(1)).dVar;
            end
            % Check for new pBest
            if p.popMember(ii).fitnessVal < p.popMember(ii).pBestFit
                % Set pBestPos
                p.popMember(ii).pBestPos = p.popMember(ii).dVar;
                % Set pBestFit
                p.popMember(ii).pBestFit = p.popMember(ii).fitnessVal;
            end
        end
        
    else % Normal movement
        
        % Choose random personal and group coeff
        rp = rand; rg = rand;
        % Move particle and evaluate fitness
        for ii = 1:count
            % Update particle velocity according to Pedersen 2010:
            p.popMember(ii).vel = omega.*p.popMember(ii).vel + ...
                pPhi.*rp.*(p.popMember(ii).pBestPos(:,2) - p.popMember(ii).dVar(:,2)) + ...
                gPhi.*rg.*(p.popMember(ii).gBestPos(:,2) - p.popMember(ii).dVar(:,2));
            % Bound velocity
            if p.popMember(ii).vel(p.popMember(ii).vel > 139)
                p.popMember(ii).vel(p.popMember(ii).vel > 139) = 139;
            else
                if p.popMember(ii).vel(p.popMember(ii).vel < -139)
                    p.popMember(ii).vel(p.popMember(ii).vel < -139) = -139;
                end
            end
            % Move particle to proposed position
            p.popMember(ii).dVarProposed = p.popMember(ii).dVar(:,2) +  p.popMember(ii).vel;
            % Bound position
            if any(p.popMember(ii).dVarProposed > 140)
                p.popMember(ii).dVarProposed(p.popMember(ii).dVarProposed > 140,1) = 140;
            end
            
            if any(p.popMember(ii).dVarProposed < 1)
                p.popMember(ii).dVarProposed(p.popMember(ii).dVarProposed < 1,1) = 1;
            end
            % Generate new element distribution
            p.popMember(ii) = changeVar(p.popMember(ii),p.popMember(ii).dVarProposed);
            % Only evaluate valid movements
            if p.popMember(ii).validConnec == 1 && p.popMember(ii).validIsland == 1
                % Assign proposed values
                p.popMember(ii).dVar(:,2) = p.popMember(ii).dVarProposed;
                p.popMember(ii).elDist = p.popMember(ii).elDistProposed;
                p.popMember(ii).nodalCoords = p.popMember(ii).nodalCoordsProposed;
                % Update position dependent properties
                p.popMember(ii) = newPositionProps(p.popMember(ii));
                % Call the fitnessEval method
                p.popMember(ii) = fitnessEval(p.popMember(ii));
            else
                % Clear proposed variables
                clear p.popMember(ii).dVarProposed
                clear p.popMember(ii).elDistProposed
                clear p.popMember(ii).nodalCoordsProposed
                % Update the fitness components for this iteration
                p.popMember(ii).fitnessValComponents = ...
                    [p.popMember(ii).fitnessValComponents; p.popMember(ii).fitnessValComponents(end,:)];
            end
        end
    end
    
    % Update pBest and gBest
    % Find gBestFit and gWorstFit
    bestPar = find([p.popMember.fitnessVal] == min([p.popMember.fitnessVal]));
    worstPar = find([p.popMember.fitnessVal] == max([p.popMember.fitnessVal]));
    
    % Update particle gBest and pBest
    for ii = 1:count
        % Check for new gBest
        if p.popMember(bestPar(1)).fitnessVal < p.popMember(ii).gBestFit
            % Set gBestFit
            p.popMember(ii).gBestFit = p.popMember(bestPar(1)).fitnessVal;
            % Set gBestPos
            p.popMember(ii).gBestPos = p.popMember(bestPar(1)).dVar;
        end
        % Check for new pBest
        if p.popMember(ii).fitnessVal < p.popMember(ii).pBestFit
            % Set pBestPos
            p.popMember(ii).pBestPos = p.popMember(ii).dVar;
            % Set pBestFit
            p.popMember(ii).pBestFit = p.popMember(ii).fitnessVal;
        end
    end
    
    % Find population fitness mean, best, and worst
    popFitMean = mean([p.popMember.fitnessVal]);
    popFitCurrentBest = p.popMember(bestPar(1)).fitnessVal;
    popFitGlobalBest = p.popMember(1).gBestFit;
    popFitWorst = p.popMember(worstPar(1)).fitnessVal;
    % Append population fitness statistics to history vector
    popFitHist(iter+1,:) = [popFitGlobalBest popFitCurrentBest popFitMean popFitWorst];
    
    % Render a selected particle
    if (pRender)
        % Render the shape of the best particle
        fig1 = figure(1);
        cla
        Plot2DGeometryUndeformed(p.popMember(bestPar(1)).elDist,p.popMember(bestPar(1)).nodalCoords,0,['Particle ' num2str(bestPar)]);
        drawnow limitrate
        
        % Graph the population fitness
        fig2 = figure(2);
        cla
        plot(1:1:iter,popFitHist(1:iter,1),1:1:iter,popFitHist(1:iter,2),1:1:iter,popFitHist(1:iter,3));
        ylim([0 10]);
        title(['Population Fitness, Iteration ' num2str(iter)]);
        legend('Global Best','Current Best','Current Mean');
        drawnow limitrate
        
        % Show the velocity term for the first particle
        fig4 = figure(4);
        cla
        hold on
        bar(1:1:size(p.popMember(1).vel,1),p.popMember(1).vel);
        ylim([-139 139]); xlim([1 size(p.popMember(1).vel,1)]);
        title(['Particle 1 Velocity, Iteration ' num2str(iter)]);
        drawnow limitrate
        
        % Show the velocity term for the first particle
        fig5 = figure(5);
        cla
        hold on
        bar(1:1:size(p.popMember(bestPar(1)).vel,1),p.popMember(bestPar(1)).vel);
        ylim([-139 139]); xlim([1 size(p.popMember(bestPar(1)).vel,1)]);
        title(['Particle ' num2str(bestPar(1)) ' Velocity, Iteration ' num2str(iter)]);
        drawnow limitrate
        
        % Write into movie structure
        fitMOV(iter) = getframe(fig2);
        parMOV(iter) = getframe(fig1);
        velMOV(iter) = getframe(fig4);
        vebMOV(iter) = getframe(fig5);
    end
    
    % Terminal write outs
    % Write out iteration count
    fprintf('Iteration %6u of %6u ~~~~~~~~~~~~~~~~~~~~~~~~\n',iter,iterLimit);
    % Write out population fitness statistics
    fprintf(['\tGroup Best fitness:\t\t\t\t' num2str(popFitGlobalBest) '\n']);
    fprintf(['\tCurrent Best fitness:\t\t\t' num2str(popFitCurrentBest) '\n']);
    fprintf(['\tPopulation fitness mean:\t\t' num2str(popFitMean) '\n']);
    fprintf(['\tWorst fitness:\t\t\t\t\t' num2str(popFitWorst) '\n']);
    
    % Update iteration count
    iter = iter + 1;
    
end

%% Set to best positions

for ii = 1:count
    % Set current position to personal best
    p.popMember(ii).dVarProposed = p.popMember(ii).pBestPos(:,2);
    p.popMember(ii) = changeVar(p.popMember(ii),p.popMember(ii).dVarProposed);
    % Assign best position values
    p.popMember(ii).dVar(:,2) = p.popMember(ii).dVarProposed;
    p.popMember(ii).elDist = p.popMember(ii).elDistProposed;
    p.popMember(ii).nodalCoords = p.popMember(ii).nodalCoordsProposed;
    % Calculate new position-dependent properties
    p.popMember(ii) = newPositionProps(p.popMember(ii));
    % Evaluate fitness
    p.popMember(ii) = fitnessEval(p.popMember(ii));
end

%% Save solution set

% Save the final population
save([fis,'_solution'],'-v7.3','p','popFitHist','iter','omega','pPhi','gPhi');

% Write movis to gif
fitGifName = [fis '_fitnessAnimated.gif']; % Specify the output file name
parGifName = [fis '_bestParAnimated.gif'];
velGifName = [fis '_velocityAnimated.gif'];
vebGifName = [fis '_velocityBestParticleAnimated.gif'];
for idx = 1:(iter-1)
    disp(['Saving GIFs Frame ' num2str(idx)]);
    [aFIT,mapFIT] = rgb2ind(fitMOV(idx).cdata,256);
    [aPAR,mapPAR] = rgb2ind(parMOV(idx).cdata,256);
    [aVEL,mapVEL] = rgb2ind(velMOV(idx).cdata,256);
    [aVEB,mapVEB] = rgb2ind(vebMOV(idx).cdata,256);
    if idx == 1
        imwrite(aFIT,mapFIT,fitGifName,'gif','LoopCount',Inf,'DelayTime',1);
        imwrite(aPAR,mapPAR,parGifName,'gif','LoopCount',Inf,'DelayTime',1);
        imwrite(aVEL,mapVEL,velGifName,'gif','LoopCount',Inf,'DelayTime',1);
        imwrite(aVEB,mapVEB,vebGifName,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(aFIT,mapFIT,fitGifName,'gif','WriteMode','append','DelayTime',0.2);
        imwrite(aPAR,mapPAR,parGifName,'gif','WriteMode','append','DelayTime',0.2);
        imwrite(aVEL,mapVEL,velGifName,'gif','WriteMode','append','DelayTime',0.2);
        imwrite(aVEB,mapVEB,vebGifName,'gif','WriteMode','append','DelayTime',0.2);
    end
end
