classdef particle
    % Initialization, manipulation, and evaluation of a particle object
    %
	%	written by Matt Ireland 2019
	%
    %   The constructor method finds a valid random arrangement of design
    %   variables, builds the corresponding element distribution, and finds
    %   the position-dependent properties by calling:
    %       1. elementPotential
    %       2. dvarPlacement
    %       3. buildProps
    %       4. pInterp
    %
    %   Property and method blocks are listed generally in order of
    %   dependence; if the particle's position is changed via the
    %   optimization block, methods in the positional block will need to be
    %   called before methods from the simulation block reflect the new
    %   position, etc.
    %
    %   The optimization block contains meta-methods that call the methods
    %   from other blocks to simplify the upstream syntax
    
    
    % Fundamental Defining Properties
    % These properties are set when the particle is initialized
    properties (SetAccess = immutable)
        nndx
        nndy
        scale
        type
        pload
        nodof
        eldof
        elPot
        nomProp
        initVar
        velPot
        gsStrainEnergy
    end
    
    % Positional Defining Properties
    % These properties are dependent upon the position of the particle in
    % the design space
    properties
        scan
        maxRed
        scanDir
        elDist
        randVar
        graphRep
        nodalCoords
        nel
        prop
        nodalFreeMat
        numFreeDOF
        numNodes
        imposedLoad
    end
    
    % Simulated Properties
    % These properties are i/o for FE simulation
    properties
        KK
        F
        delta
        node_disp
    end
    
    % Optimization-Modified Properties
    % These properties are i/o for external interaction with the particle
    properties
        dVar
        dVarProposed
        elDistProposed
        nodalCoordsProposed
        validConnec
        validIsland
        fitnessVal
        fitnessValComponents
        pBestFit
        gBestFit
        pBestPos
        gBestPos
        vel
        memory = struct('pos',[],'fit',[],'pbf',[],'pbp',[])
    end
    
    % Particle initialization
    methods
        % Constructor method
        function obj = particle(nndx,nndy,scale,E,A,I,scanDir,scan,maxRed,type,pload)
            % Generate an element distribution and evaluate the FE problem
            
            if nargin ~=0 % allow for the no-input argument requirement
                % Assign Properties
                    % This section assigns the immutable properties of the object
                    % that can't be modified
                    obj.nndx = nndx;
                    obj.nndy = nndy;
                    obj.scale = scale;
                    obj.scan = scan;
                    obj.maxRed = maxRed;
                    obj.type = type;
                    obj.pload = pload;
                    % Write the section properties into a class property structure
                    obj.nomProp = [E A I];
                    % Generate the element potential matrix property
                    obj.elPot = elementPotential(nndx,nndy);
                    % Generate the particle velocity potential array
                    obj.velPot = -139:1:139;
                
                % Frame Element Definition
                    % This section provides element definitions that control data
                    % structure sizing and certain indexing actions elsewhere.
                    % Modifying these numbers will break the simulation, and so
                    % this code should be cleaned up
                    nne = 2 ; % Number of nodes per element
                    obj.nodof = 3 ; % Number of degrees of freedom per node
                    obj.eldof = nne * obj.nodof; % Number of degrees of freedom per element                    
                    
                % Ground Structure Calculation
                    % This section calculates the strain energy of the
                    % "fully dense", or ground structure condition for
                    % normalization of strain energy during optimization
                    obj.gsStrainEnergy = groundStructure(obj);
                    
                    
                % Initialize Particle Position
                    % This section generates position-dependent properties that are allowed to
                    % change after initialization
                    % Write the scan direction into a class property
                    obj.scanDir = scanDir;
                    % Generate an element distribution from a random selection
                    % of design values for each unit cell
                    [obj.elDist,obj.randVar,~,obj.nodalCoords] = dvarPlacement(obj.nndx,obj.nndy,obj.scale,obj.type);
                    % Assign the random dVar distribution to the property
                    % that will be changed during optimization
                    obj.dVar = obj.randVar;
                    % Record the initial position
                    obj.initVar = obj.randVar;
                    % Count number of elements
                    obj.nel = size(obj.elDist,1);
                    % Build the element property matrix
                    obj.prop = buildProps(obj.nel,E,A,I);
                    % Modify element stiffness
                    if obj.scan == 1
                        obj.prop = pInterp(obj.prop,obj.elDist,obj.nodalCoords,scanDir,maxRed);
                    end
            end
        end
    end
    
    % Methods for generating positional properties
    methods
        % Build the element property matrix
        function obj = buildProps(obj)
            % Call the element property assignment function
            obj.prop = buildProps(obj.nel,obj.nomProp(1),obj.nomProp(2),obj.nomProp(3));
        end
        
        % Impose boundary conditions
        function obj = boundaryCond(obj)
            % Call the boundary condition function
            [obj.nodalFreeMat,obj.numFreeDOF,obj.numNodes,obj.imposedLoad] = ...
                boundaryCond(obj.nndx,obj.nndy,obj.scale,obj.nodof,obj.nodalCoords,obj.type,obj.pload);
        end
        
        % Modify element stiffness with process interpretation
        function obj = pInterp(obj)
            % Call the process interpreter function
            obj.prop = pInterp(obj.prop,obj.elDist,obj.nodalCoords,obj.scanDir,obj.maxRed);
        end
    end
    
    % Methods for generating simulated properties
    methods
        % Build the global stiffness matrix
        function obj = buildGlobalStiffness(obj)
            % Call the global stiffness function
            obj.KK = buildGlobalStiffness(obj.numFreeDOF,obj.nel,obj.nodalCoords,obj.elDist,obj.prop,obj.nodalFreeMat,obj.eldof);
        end
        
        % Build the global force vector
        function obj = formGlobalF(obj)
            % Call the global force vector function
            obj.F = formGlobalF(obj.numFreeDOF,obj.numNodes,obj.nodof,obj.nodalFreeMat,obj.imposedLoad);
        end
        
        % Solve the FE system of equations
        function obj = feSolve(obj)
            % Call the finite element solver
            [obj.delta,obj.node_disp] = feSolve(obj.numNodes,obj.nodof,obj.nodalFreeMat,obj.KK,obj.F);
        end
    end
    
    % Methods for optimization interaction with the particle
    methods
        % Generate the strain energy of the ground structure
        function strainEnergy = groundStructure(obj)
            % Find the number of unit cells
            ncell = (obj.nndx-1)*(obj.nndy-1);
            % Initialize the ground structure value matrix
            groundVar = ones(ncell,1);% [n var]
%             groundVar(:,1) = 1:1:ncell;% unit cell numbering
            groundVar(:,1) = groundVar(:,1).*140;% set to fully dense
            % Build dependent properties
            obj = changeVar(obj,groundVar);
            obj.dVar = obj.dVarProposed;
            obj.elDist = obj.elDistProposed;
            obj.nodalCoords = obj.nodalCoordsProposed;
            % Count number of elements
            obj.nel = size(obj.elDist,1);
            % Build the element property matrix
            obj.prop = buildProps(obj.nel,obj.nomProp(1),obj.nomProp(2),obj.nomProp(3));
            % Simulate performance
            obj = simulate(obj);
            % Calculate ground structure strain energy
            strainEnergy = (obj.delta'*obj.KK*obj.delta);
        end
        
        % Generate a random distribution of design variables
        function obj = randomPosition(obj)
            % Call the design variable placement function
            [obj.elDistProposed,obj.randVar,~,obj.nodalCoordsProposed] = dvarPlacement(obj.nndx,obj.nndy,obj.scale,obj.type);
            % Set the valid config flags
            obj.validConnec = 1; obj.validIsland = 1;
            % Write the random variables to the propsed position
            obj.dVarProposed = obj.randVar(:,2);
        end
        
        % Generate a new element distribution from design values
        function obj = changeVar(obj,dVarProposed)
            % A means of updating the element distribution according to a
            % new position vector of design variables. Sets flags for
            % element distributions that don't meet the connectivity and
            % continuity conditions.
            
            % Initialize
            elDistNew = [];
            
            % Update the element distribution
                % for all of the unit cells
                for ii = 1:size(dVarProposed,1)
                    % Parse new design variable value into line elements for the unit cell
                    unitElements = dvarParse(obj.elPot(ii,:),dVarProposed(ii));
                    % For the nonzero case
                    if unitElements(1,1) ~= 0
                        % Change the local element numbering to the global values
                        if ii > 1 && size(elDistNew,1) > 0
                            unitElements(:,1) = unitElements(:,1) + elDistNew(end,1);
                        end
                        % Concatenate to the element distribution matrix
                        elDistNew = [elDistNew; unitElements];
                    end
                end
                % Write to the proposed elDist property
                obj.elDistProposed = elDistNew;
                
            % Check for the zero position
                if isempty(elDistNew)
                    % Trigger both flags
                    obj.validConnec = 0;
                    obj.validIsland = 0;
                    % Exit the function
                    return
                end
                
                % Check for BC connectivity and flag (connectivity condition)
                % Represent connectivity as a network object
                obj.graphRep = graph(elDistNew(:,2),elDistNew(:,3));
                
                % Create node variables to check for connectivity
                if obj.type == 1% (MBB loading)
                    n1 = 1; n2 = obj.nndx; n3 = obj.nndx *(obj.nndy-1) + 1;
                else
                    if obj.type == 2% (cantilever midplane tip loading)
                        n1 = 1; n2 = obj.nndx * ((obj.nndy-1)/2 + 1); n3 = obj.nndx *(obj.nndy-1) + 1;
                    end
                end
                % Put network components into bins
                [gbins,gbinsizes] = conncomp(obj.graphRep);

                % Check if vital nodes are in the same component (connected)
                if size(gbins,2) >= n3 && gbins(n1) == gbins(n2) == gbins(n3)
                    % No flag
%                     disp('Valid Manipulation');
                    obj.validConnec = 1;
                else
                    % Flag the change
%                     disp('Invalid Manipulation');
                    obj.validConnec = 0;
                end

            
            % Check for islands and flag (continuity condition)
                % Initialize flag
                obj.validIsland = 1;
                
                % Build logical array of invalid components
                h = gbinsizes < max(gbinsizes) & gbinsizes > 1;

                % Search through nodes
                for ii = 1:length(gbins)
                    % If node belongs to invalid component
                    if h(gbins(ii)) == 1
                        % Flag the change
                        obj.validIsland = 0;
                    end
                end

            
            % Make a new nodal coordinate matrix
                % Build a nodal location grid
                grid = nodalGrid(obj.nndx,obj.nndy,obj.scale);

                % Build vector of nodes connected by edges
                nodeList = cat(1,elDistNew(:,2),elDistNew(:,3));

                % Remove duplicate node entries
                nodeList = unique(nodeList);

                % Sort the nodes
                nodeList = sort(nodeList);

                % Initialize a location matrix
                obj.nodalCoordsProposed = zeros(size(nodeList,1),3);
                obj.nodalCoordsProposed(:,1) = nodeList;

                % Search through the location grid and apply coordinates
                for ii = 1:size(nodeList,1)
                    obj.nodalCoordsProposed(ii,2:3) = grid(nodeList(ii),2:3);
                end
        end
        
        % Re-calculate the position-dependent properties
        function obj = newPositionProps(obj)
            % Count number of elements
            obj.nel = size(obj.elDist,1);
            % Build the element property matrix
            obj = buildProps(obj);
            % Impose boundary conditions
            obj = boundaryCond(obj);
            % Modify element stiffness with process interpretation
            if obj.scan == 1
                obj = pInterp(obj);
            end
        end
        
        % Perform the FE simulation
        function obj = simulate(obj)
            % Build the boundary conditions
            obj = boundaryCond(obj);
            % Build the global stiffness matrix
            obj = buildGlobalStiffness(obj);
            % Build the global force vector
            obj = formGlobalF(obj);
            % Solve the FE system of equations
            obj = feSolve(obj);
        end
        
        % Evaluate the particle fitness
        function obj = fitnessEval(obj)
           % Perform the FE simulation
           obj = simulate(obj);
           
           % Note: This formulation of particle fitness is subject to
           % change and not necessarily the best formulation
               
               % Use the normalized strain energy and add the volume
               % fraction
               strainEnergy = (1 - obj.gsStrainEnergy/(obj.delta'*obj.KK*obj.delta));
               volumeFraction = obj.nel/(size(obj.elPot,1)*6);
               
               a = 5; b = 5; c = 5; d = 5;
               
               % Penalize a density above a threshold
               if volumeFraction > 0.4
                   obj.fitnessVal = a*strainEnergy + b*volumeFraction;
               else
                   obj.fitnessVal = c*strainEnergy + d*volumeFraction;
               end
               
               % Append to a history matrix
               obj.fitnessValComponents = [obj.fitnessValComponents; obj.fitnessVal strainEnergy volumeFraction];
        end
        
        % Memorize current positional state
        function obj = memorize(obj)
            
            % Append positional information
            obj.memory.pos = [obj.memory.pos obj.dVar(:,2)];
            % Append fitness value
            obj.memory.fit = [obj.memory.fit obj.fitnessVal];
        end
    end
    
    % Methods for data extraction
    methods
        % Render the non-deformed particle
        function [] = renderParticle(obj,id)
            % Call the Plot2DGeometry function and pass an index value for
            % a title
            Plot2DGeometryUndeformed(obj.elDist,obj.nodalCoords,0,['Particle ' num2str(id)]);
        end
    end
    
end

