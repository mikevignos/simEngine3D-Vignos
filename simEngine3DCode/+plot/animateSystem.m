function animateSystem(sys,viewAngle)
% simple animation of system3D.
% INPUTS:
%   sys   : system3D object. Contains output from simulation
%   viewAngle: [Az, El] to set viewing angle of animation figure, optional

if ~exist('viewAngle','var') || isempty(viewAngle)
    viewAngle = [98,12];
end

% determine time step
timeStep = sys.myBodies{1}.myTimeTotal(2) - sys.myBodies{1}.myTimeTotal(1);

% determine number of free bodies in simulation
nFreeBodies = sys.myNumBodiesMinusGround;

% determine number of points on bodies
nFreePoints = 0;
for i = 1:nFreeBodies
    nFreePoints = nFreePoints + sys.myBodies{i+1}.myNumPoints;
end


% allocate for speed
bodyHandles = gobjects(nFreeBodies,1); % empty graphics handles array
pointHandles = gobjects(nFreePoints,1);
bodyColors = zeros(nFreeBodies,3);
% rBodies = sys.myRTotal; % body positions
time = sys.myBodies{1}.myTimeTotal;
rBodies = cell(length(time),1); % body positions
rPoints = cell(length(time),1); % point positions (as lines)

%% PULL DATA
% iterate over the time grid  
for i = 1:length(time)
    pointNum = 1;
    rBody = cell(nFreeBodies,3);
    rPoint = cell(nFreePoints,3);
    for j = 1:nFreeBodies % plot free bodies
        
        % save color of body
%         bodyID = sys.bodyIDs(j); % pull current free-body ID
%         bodyColors(j,:) = sys.body{bodyID}.color;
        
        % save position of body
        rBody{j,1} = sys.myBodies{j+1}.myRTotal(1,i)';
        rBody{j,2} = sys.myBodies{j+1}.myRTotal(2,i)';
        rBody{j,3} = sys.myBodies{j+1}.myRTotal(3,i)';
        r = cell2mat(rBody(j,:));
        
        % save points as lines
        p = sys.myBodies{j+1}.myPTotal(:,i);
        A = simEngine3DUtilities.p2A(p);
        if sys.myBodies{j+1}.myNumPoints ~= 0
            for k = 1:sys.myBodies{j+1}.myNumPoints
                sbar = sys.myBodies{j+1}.myPoints{k}.sBar; % local position of point
                Asbar = A*sbar; % rotated to global rf
                rAsbar = r + Asbar'; % global position of point
                rPoint(pointNum,:) = {[r(1);rAsbar(1)],[r(2);rAsbar(2)],[r(3);rAsbar(3)]};
                pointNum = pointNum + 1;
            end
        end
    end
    rBodies{i} = rBody;  % save body  positions for this timestep
    rPoints{i} = rPoint; % save point positions for this timestep
end

%% ANIMATE DATA

% set frame
figure();
fig = gcf;
fig.Color = [1 1 1]; % set background color to white

% plot ground bodies 
hold on
for j = 1:sys.myNumBodies 
    if sys.myBodies{j}.myIsGround
        % plot grounded bodies as frames
        plot.drawframe(sys.myBodies{j}.myRTotal(:,1),sys.myBodies{j}.myPTotal(:,1),[],2) 
        
        % plot points on ground bodies
        A = simEngine3DUtilities.p2A(sys.myBodies{j}.myP);
        r = sys.myBodies{j}.myR;
        if sys.myBodies{j}.myNumPoints ~= 0
            for k = 1:sys.myBodies{j}.myNumPoints 
                sbar = sys.myBodies{j}.myPoints{k}.sBar; % local position of point
                Asbar = A*sbar; % rotated to global rf
                rAsbar = r + Asbar; % global position of point
                plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-','MarkerSize',12); % line from BODY RF to point
            end
        end        
    end
end

% plot initial bodies
for j = 1:nFreeBodies 
%     bodyHandles(j) = scatter3(sys.myBodies{j+1}.myRTotal(1,1),sys.myBodies{j+1}.myRTotal(2,1),sys.myBodies{j+1}.myRTotal(2,1),100,'k','filled');
    bodyHandles(j) = scatter3(rBodies{1}{j,1},rBodies{1}{j,2},rBodies{1}{j,3},100,'k','filled');
end

% % plot initial points
for j = 1:nFreePoints
    pointHandles(j) = plot3(rPoints{1}{j,1},rPoints{1}{j,2},rPoints{1}{j,3},'ks-','MarkerSize',12); % line from BODY RF to point
end
hold off

% determine size of axes
maxs = [0,0,0]; mins = [0,0,0];
for i = 1:length(time) % find max and min positions over time
    rTemp = sys.myRTotal(:,i);
    R = reshape(rTemp,[3,sys.myNumBodies])';
    maxs = max([maxs; R]);
    mins = min([mins; R]);
end
border = 0.5*max((maxs - mins));
% border = 1;
axisWindow = [mins(1)-border maxs(1)+border mins(2)-border maxs(2)+border mins(3)-border maxs(3)+border];

% set frame properties
axis equal
axis(axisWindow); % set axes size
view(viewAngle(1),viewAngle(2)); 
pause;

% animate through time
propertyCell = {'XData','YData','ZData'};
for i = 1:length(time)
    % update bodies
    set(bodyHandles,propertyCell,rBodies{i});
    
    %     % update points on bodies
    set(pointHandles,propertyCell,rPoints{i});
    
    drawnow; % redraw updated points on the plot
    %pause(timeStep); % wait an actual time step
end


end