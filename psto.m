clc% Clear screen
clear% Clear all variables in memory

% Particle Swarm Topology Optimization
%         Matt Ireland 2019

%% Define Population

% Generate or Load a population
gen = 1;% 0 to load, 1 to gen

% Load a population
fid = '20200124-1510_pop';% name of the saved .mat population

if gen == 1
    % Generate a population
    count = 160;% the population size
    
    % Input Data
    % General parameters; Design space
    nndx = 3; % Number of nodes in x
    nndy = 2; % Number of nodes in y
    scale = 5; % Scaling factor for element size
    
    % Boundary conditions
    % Choose BC type
    type = 1;   % 1 for MBB
                % 2 for Cantilever midplane tip-load
    
    % Set point load
    pload = 1000;% unitless
    
    % Element mechanical properties
    % Modulus, section area, and MOI
    E = 100000; A = 1; I = pi/4;
    
    % Process interpretation
    % Scan direction
    scanDir = [0 1];% [x y]
    scan = 0;% logical yes/no
    maxRed = 0.2;% maximum percent reduction
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

%% File Naming
% Automatically build filenames based on input parameters

% Date and time
fis = char(datetime('now','Format','yyyyMMdd''-''HHmm'));
% Grid Size
fis = [fis,'_',num2str(nndx),'x',num2str(nndy),'grid_'];
% Load type
if type == 1
    fis = [fis,'mbb'];
else 
    if type == 2
        fis = [fis,'canti'];
    end
end
% Scan
if scan
   if scanDir(1) == 1
       fis = [fis,'ScanX_'];
   else
       if scanDir(2) == 1
           fis = [fis,'ScanY_'];
       end
   end
end

%% Optimization Parameters

% Convergence Criteria
ctol = 1e-3;

% Number of pre-scatter iterations to keep for convergence calc
rollKeep = 7;

% Scattering trigger length
sTrig = 10;

% Iteration limit
iterLimit = 10000;

% Velocity maximum
vmax = 10;

% Velocity Update Tuning
omega = 0.4;
pPhi = 0.2;
gPhi = 0.6;

%% Rendering
% Render best particle, fitness, and velocity values
pRender = 1;% set to zero for none

if (pRender)
    % Initialize structures for animations
    fitMOV(iterLimit) = struct('cdata',[],'colormap',[]);
    parMOV(iterLimit) = struct('cdata',[],'colormap',[]);
    posMOV(iterLimit) = struct('cdata',[],'colormap',[]);
    vebMOV(iterLimit) = struct('cdata',[],'colormap',[]);
    
    % Change figure background color to white
    get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]);
end

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
    p.popMember(ii).gBestFit = p.popMember(bestPar(1)).fitnessVal;
    % Set gBestPos
    p.popMember(ii).gBestPos = p.popMember(bestPar(1)).dVar;
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
popFitHist(1,:) = [p.popMember(bestPar(1)).fitnessVal p.popMember(bestPar(1)).fitnessVal...
    popFitMean p.popMember(worstPar(1)).fitnessVal];
% Initialize iteration counter
iter = 1;
% Initialize fitness change checksum
popFitChange = 1;
% Initialize a scattering count
nScatter = 0;
% Initialize a position plotting array
popPos = ones(size(p.popMember(1).dVar,1),count);
% Initialize a scatter boolean
sctBool = 0;

%% Iteration

while iter <= iterLimit
    
    % Check if current interval needs scattering
    if (sctBool)
        % Update the number of scatterings
        nScatter = nScatter + 1;
        % Inspect a rolling fitness vector for convergence
        if nScatter <= rollKeep
            % Write into rolling vector
            popFitRolling(nScatter) = popFitCurrentBest(1);
        else
            % Trim vector
            popFitRolling = [popFitRolling(2:end); popFitCurrentBest(1)];
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
            % Data structure for population position plotting
            popPos(:,ii) = p.popMember(ii).dVar(:,2);
        end
        
        % Find gBestFit and gWorstFit
        bestPar = find([p.popMember.fitnessVal] == min([p.popMember.fitnessVal]));
        worstPar = find([p.popMember.fitnessVal] == max([p.popMember.fitnessVal])); %#ok<NASGU>
        
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
        
        % Reset scatter boolean
        sctBool = 0;
        
    else % Normal movement
        
        % Move particle and evaluate fitness
        for ii = 1:count
            % Choose random personal and group coeff
            rp = rand; rg = rand;
            % Update particle velocity according to Pedersen 2010:
            p.popMember(ii).vel = omega.*p.popMember(ii).vel + ...
                pPhi.*rp.*(p.popMember(ii).pBestPos(:,2) - p.popMember(ii).dVar(:,2)) + ...
                gPhi.*rg.*(p.popMember(ii).gBestPos(:,2) - p.popMember(ii).dVar(:,2));
            % Round to nearest integer
            p.popMember(ii).vel = round(p.popMember(ii).vel);
            % Bound velocity
            if p.popMember(ii).vel(p.popMember(ii).vel > vmax)
                p.popMember(ii).vel(p.popMember(ii).vel > vmax) = vmax;
            else
                if p.popMember(ii).vel(p.popMember(ii).vel < -vmax)
                    p.popMember(ii).vel(p.popMember(ii).vel < -vmax) = -vmax;
                end
            end
            % Move particle to proposed position
            p.popMember(ii).dVarProposed = p.popMember(ii).dVar(:,2) +  p.popMember(ii).vel;
            % Bound position
            if any(p.popMember(ii).dVarProposed > 56)
                p.popMember(ii).dVarProposed(p.popMember(ii).dVarProposed > 56,1) = 56;
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
            % Data structure for population position plotting
            popPos(:,ii) = p.popMember(ii).dVar(:,2);
            
            % Check if next iteration should be scattered
            if iter > sTrig
                % If the current and global best nad been equal too long
                if popFitHist(iter-sTrig:iter,1) == popFitHist(iter-sTrig:iter,2)
                    sctBool = 1;
                else
                    % If the current and global best aren't converging
                    if range(popFitHist(iter-sTrig:iter,1)) < ctol && range(popFitHist(iter-sTrig:iter,2)) < ctol
                        sctBool = 1;
                    end
                end
            end
        end
    end
    
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
    
    % Record to terminal variables
    popFitMean = mean([p.popMember.fitnessVal]);
    popFitCurrentBest = p.popMember(bestPar(1)).fitnessValComponents(iter+1,:);
    popFitGlobalBest = p.popMember(1).gBestFit;
    popFitWorst = p.popMember(worstPar(1)).fitnessVal;
    % Append population fitness statistics to history vector
    popFitHist(iter+1,:) = [popFitGlobalBest popFitCurrentBest(1) popFitMean popFitWorst];
    
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
        ylim([0 2]);
        title(['Population: Objective Function Value, Iteration ' num2str(iter)]);
        legend('Global Best','Current Best','Current Mean');
        xlabel('Iteration'); ylabel('Objective Value');
        drawnow limitrate
        
        % Show the positions of the population
        fig4 = figure(4);
        cla
        hold on
        [nn,~] = meshgrid(1:1:size(p.popMember(1).dVar,1),1:1:count);
        nn = reshape(nn',[],1); popPos = reshape(popPos,[],1);
        if (sctBool)
            scatter(nn,popPos,'r.');
        else
            scatter(nn,popPos,'k.');
        end
        ylim([1 56]); xlim([0 size(p.popMember(1).vel,1)]);
        title(['Population Positions, Iteration ' num2str(iter)]);
        xlabel('Unit Cell Number'); ylabel('Unit Cell Value');
        drawnow limitrate; clear popPos;
        
        % Plotting data structures for patch function
        nUnCell = (nndx-1)*(nndy-1); elPot = elementPotential(nndx,nndy); grd = nodalGrid(nndx,nndy,scale);
        x = grd(elPot(:,2:5),2); x = reshape(x,nUnCell,4); x = x'; x(3:4,:) = circshift(x(3:4,:),1,1);
        y = grd(elPot(:,2:5),3); y = reshape(y,nUnCell,4); y = y';
        % Remap unit cell values to "density" values
        dens = p.popMember(bestPar(1)).dVar(:,2);
        dens(dens == 1) = 0; dens(dens > 1 & dens <= 7) = 1; dens(dens > 7 & dens <= 22) = 2;
        dens(dens > 22 & dens <= 42) = 3; dens(dens > 42 & dens <= 49) = 4;
        dens(dens > 49 & dens <= 55) = 5; dens(dens == 56) = 6;
        
        % Show the unit cell value of the best particle
        fig5 = figure(5);
        cla
        hold on
        patch(x,y,dens);
        colormap(flipud(gray(7)));
        xlim([0 (nndx-1)*scale]); ylim([0 (nndy-1)*scale]); pbaspect([nndx nndy 1]);
        set(gca,'CLim',[0 7]); c = colorbar('Ticks',[0 1 2 3 4 5 6]); c.Label.String = 'Number of Elements';
        title(['Particle ' num2str(bestPar(1)) ', Element Count, Iteration ' num2str(iter)]);
        drawnow limitrate
        
        % Write into movie structure
        fitMOV(iter) = getframe(fig2);
        parMOV(iter) = getframe(fig1);
        posMOV(iter) = getframe(fig4);
        vebMOV(iter) = getframe(fig5);
    end
    
    % Terminal write outs
    % Write out iteration count
    fprintf('Iteration %6u of %6u\n',iter,iterLimit);
    % Write out population fitness statistics
    fprintf('\tGlobal Best fitness:\t\t\t%.5f\tSE-----\tVolPen-\tVolTar-\tVol----\n',popFitGlobalBest);
    fprintf('\tCurrent Best fitness:\t\t\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n',popFitCurrentBest(1:5));
    fprintf('\tPopulation fitness mean:\t\t%.5f\t-------\t-------\t-------\t-------\n',popFitMean);
    fprintf('\tWorst fitness:\t\t\t\t\t%.5f\t-------\t-------\t-------\t-------\n',popFitWorst);
    
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
if (pRender)
    fitGifName = [fis '_fitnessAnimated.gif']; % Specify the output file name
    parGifName = [fis '_bestParAnimated.gif'];
    posGifName = [fis '_populationPositionAnimated.gif'];
    vebGifName = [fis '_velocityBestParticleAnimated.gif'];
    for idx = 1:(iter-1)
        disp(['Saving GIFs Frame ' num2str(idx)]);
        [aFIT,mapFIT] = rgb2ind(fitMOV(idx).cdata,256);
        [aPAR,mapPAR] = rgb2ind(parMOV(idx).cdata,256);
        [aPOS,mapPOS] = rgb2ind(posMOV(idx).cdata,256);
        [aVEB,mapVEB] = rgb2ind(vebMOV(idx).cdata,256);
        if idx == 1
            imwrite(aFIT,mapFIT,fitGifName,'gif','LoopCount',Inf,'DelayTime',1);
            imwrite(aPAR,mapPAR,parGifName,'gif','LoopCount',Inf,'DelayTime',1);
            imwrite(aPOS,mapPOS,posGifName,'gif','LoopCount',Inf,'DelayTime',1);
            imwrite(aVEB,mapVEB,vebGifName,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(aFIT,mapFIT,fitGifName,'gif','WriteMode','append','DelayTime',0.2);
            imwrite(aPAR,mapPAR,parGifName,'gif','WriteMode','append','DelayTime',0.2);
            imwrite(aPOS,mapPOS,posGifName,'gif','WriteMode','append','DelayTime',0.2);
            imwrite(aVEB,mapVEB,vebGifName,'gif','WriteMode','append','DelayTime',0.2);
        end
    end
end
