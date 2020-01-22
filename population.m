classdef population
    % Generate an object array of particles
	%
	%	written by Matt Ireland 2019
    
	% Population properties
    properties
        count
        popMember
        time
    end
    
	% Initialize population
    methods
        function obj = population(count,nndx,nndy,scale,E,A,I,scanDir,scan,maxRed,type,pload)
			% Allow for the zero-input condition
            if nargin == 0
                obj.count = 0;
                obj.popMember = [];
                obj.time = [];
            end
			% Record the population size into the objects
            obj.count = count;
			% This loop can be changed to a parfor if needed
            for ii = 1:obj.count
                tic;
				% Initialize particle
                p(ii) = particle(nndx,nndy,scale,E,A,I,scanDir,scan,maxRed,type,pload);
				% Record particle initialization time
                time(ii) = toc;
            end
			% Record into population structure
            obj.popMember = p;
            obj.time = time;
        end
    end
    
end

