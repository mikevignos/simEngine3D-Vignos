function plotSystem(sys,frames)
% plotSystem  plots the bodies in the multibodySystem
%
%	plotSystem(sys)
%	plotSystem(sys, frames)
%
% sys        : created with the system3D class.
% frames = 1 : Show body reference frames.
%        = 0 : Don't show body reference frames.
%

if nargin == 1
    frames = 0;
end



figure()

%% PLOT GLOBAL REFERNCE FRAME
plot.drawframe([0 0 0]',[1 0 0 0]',1,1) % plot GLOBAL RF
xlabel('X');
ylabel('Y');
zlabel('Z');

%% PLOT BODIES
if sys.myNumBodies == 0; disp('No bodies in the system.'); return; end;


hold on;
for i = 1:sys.myNumBodies % plot bodies in system
    
    % plot marker for every body
    r = [sys.myBodies{i}.myR(1);sys.myBodies{i}.myR(2);sys.myBodies{i}.myR(3)];
    scatter3(r(1),r(2),r(3),500);
    
    %optionally plot body reference frames
    if frames 
        if sys.myBodies{i}.myIsGround
            plot.drawframe(sys.myBodies{i}.myR,sys.myBodies{i}.myP,[],2)
        else
            plot.drawframe(sys.myBodies{i}.myR,sys.myBodies{i}.myP)
        end
    end
    
    % add body labels
    r_text = r + 0.1; %adjust so label is not right over body origin
    text(r_text(1),r_text(2),r_text(3),['Body ' num2str(sys.myBodies{i}.myBodyNumber)],'FontWeight','bold');
    
    % IF POINTS ON BODY, PLOT THEM
    if ~isempty(sys.myBodies{i}.myNumPoints)
        for j = 1:sys.myBodies{i}.myNumPoints % plot bodies in system
            sbar = sys.myBodies{i}.myPoints{j}.sBar; % local position of point
            Asbar = simEngine3DUtilities.p2A(sys.myBodies{i}.myP)*sbar; % rotated to global rf
            rAsbar = r + Asbar; % global position of point
            
            plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-') % line from BODY RF to point
            
            % add point labels
            rAsbar_text = rAsbar + 0.05; %adjust so label is not right over point
            text(rAsbar_text(1),rAsbar_text(2),rAsbar_text(3),sys.myBodies{i}.myPoints{j}.name);
        end
    end
    
    % IF VECTORS ON BODY, PLOT THEM
    if ~isempty(sys.myBodies{i}.myNumVectors)
        for j = 1:sys.myBodies{i}.myNumVectors % plot bodies in system
            abar = sys.myBodies{i}.myVectors{j}.aBar; % local position of point
            aGlobal = simEngine3DUtilities.p2A(sys.myBodies{i}.myP)*abar; % rotated to global rf
%             rAsbar = r + Asbar; % global position of point
            
            quiver3(r(1), r(2), r(3), aGlobal(1), aGlobal(2), aGlobal(3),'r');
%             [r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-') % line from BODY RF to point
            
            % add point labels
            rAsbar_text = r + 0.05; %adjust so label is not right over point
            text(rAsbar_text(1),rAsbar_text(2),rAsbar_text(3),sys.myBodies{i}.myVectors{j}.name);
        end
    end
end
hold off;
axis equal
end