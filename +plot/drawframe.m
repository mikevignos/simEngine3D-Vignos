function drawframe(r, ORIENTATION, scale, style)
%DRAWFRAME  plots a graphical description of a coordinate frame
%
%	DRAWFRAME(r,ORIENTATION)
%	DRAWFRAME(r,ORIENTATION, SCALE)
%	DRAWFRAME(r,ORIENTATION, SCALE, style)
%
% ORIENTATION : a rotation matrix [3x3] OR euler parameters [4x1]
%           r : a point in GLOBAL RF.
%       scale : scale frame size by this constant (default = 1)
%        style : specify style of reference frame
%
% adapted from:
% $Id: drawframe.m,v 1.1 2009-03-17 16:40:18 bradleyk Exp $
% Copyright (C) 2005, by Brad Kratochvil

if ~exist('style','var') || isempty(style)
    style = 0;
end
if ~exist('scale','var') || isempty(scale)
    scale = 1;
end

% put orientation into orientation matrix
if numel(ORIENTATION) == 9 % orientation matrix
    A = ORIENTATION;
elseif numel(ORIENTATION) == 4 % euler parameters
    A = simEngine3DUtilities.p2A(ORIENTATION); % convert to orientation matrix
else
    error('Orientation parameter is the wrong size')
end



hchek = ishold;
hold on

% scaled, orthogonal unit vectors for frame
unitX = scale*[1;0;0]; 
unitY = scale*[0;1;0];
unitZ = scale*[0;0;1];

% rotated unit vectors
x = A*unitX;
y = A*unitY;
z = A*unitZ;

% choose style and plot
switch style
    case 0 
        % default,rgb frames
        quiver3(r(1),r(2),r(3),x(1),x(2),x(3),'Color',[1,0,0]);
        quiver3(r(1),r(2),r(3),y(1),y(2),y(3),'Color',[0,1,0]);
        quiver3(r(1),r(2),r(3),z(1),z(2),z(3),'Color',[0,0,1]);
    case 1
        % use black for frame color
        color = [0,0,0];
        lineWidth = 2;
        quiver3(r(1),r(2),r(3),x(1),x(2),x(3),'Color',color,'LineWidth',lineWidth);
        quiver3(r(1),r(2),r(3),y(1),y(2),y(3),'Color',color,'LineWidth',lineWidth);
        quiver3(r(1),r(2),r(3),z(1),z(2),z(3),'Color',color,'LineWidth',lineWidth);
        
        % label the frame axes
        Atext = A+scale*0.05;
        text(Atext(1,1),Atext(2,1),Atext(3,1), 'X');
        text(Atext(1,2),Atext(2,2),Atext(3,2), 'Y');
        text(Atext(1,3),Atext(2,3),Atext(3,3), 'Z');
    case 2
        % use DarkSlateGray for frame color
        color = [0.1836,0.3086,0.3086];
        quiver3(r(1),r(2),r(3),x(1),x(2),x(3),'Color',color);
        quiver3(r(1),r(2),r(3),y(1),y(2),y(3),'Color',color);
        quiver3(r(1),r(2),r(3),z(1),z(2),z(3),'Color',color);
end

if hchek == 0
    hold off
end
end