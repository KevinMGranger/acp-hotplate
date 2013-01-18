function [temp_array]  =  hot_plate(nx, ny, ta, tb, tc, td, tol)
%
%DESCRIPTION
%    Determine the steady-state temperature of a two-dimensional flat
%    metal plate, subject to a set of given boundary conditions.
%
%ARGUMENTS
%
%     nx                  is the number of sections in the
%                             X-direction 
% 
%     ny                  is the number of sections in the
%                             Y-direction 
% 
%     ta                  is the temperature of side A (the "bottom") 
%                                 (Kelvin) 
% 
%     tb                  is the temperature of side B (the "top")
%                                 (Kelvin) 
% 
%     tc                  is the temperature of side C (the "left")
%                                 (Kelvin) 
% 
%     td                  is the temperature of side D (the "right")
%                                 (Kelvin) 
% 
%     tol                 is the "tolerance": when divided by (nx * ny),
%                              the fractional error in the temperature
%                              of any element should be less than this
%
%RETURNS
%
%     temp_array          is an array of nx-by-ny elements,
%                             with the equilibrium temperature
%                             for each location on the plate
%
%AUTHOR
%    Kevin Granger <kmg2728@rit.edu>
%    2013-01-18

%{
Additional Documentation:

MATH

    This program attempts to solve the partial differential equation
    d^2 T / dx^2 = 0 = d^2 T / dy^2 iteratively. It does this by breaking 
    the plate into a given amount of small pieces in each dimension, and
    averaging the value of each piece between the pieces that surround it.
    The reasoning behind the average was determined by expressing the
    "second difference" of each of the temperature values in an equation,
    and rearranging to determine what the value of the current piece is
    when the second difference is 0.

    In other words, for a piece i compared to the pieces on the left (L),
    right (R), top (T), and bottom (B), the second difference between the
    values should be zero:

        (TR - Ti) - (Ti - TL) = 0 = (TT - Ti) - (Ti - TB) = 0

    or when rearranged:

        Ti = (TR + TL + TB + TT) / 4

    hence the average.

    New values are computed until the maximum fractional change for any
    piece is less than the determined convergence factor, which is given as
    an argument.


VARIABLE NAMING
    
    See the comments near where each variable is declared / invoked the
    first time.

%}




% Check starting values :
if min([ta tb tc td]) <= 0
    error('Temperatures are given in Kelvin, and as such must be above absolute zero.')
elseif (min([nx ny]) < 1) || (max(rem([nx ny],1)) ~= 0)
    error('You must give a positive, nonzero, integer number of pieces to break the rod into.')
end

% Populate variables :

% These are the temperature arrays.
% Since the averaging for each piece is done using the old values, two
% arrays are necessary.
old_array = ones(ny+2,nx+2) * mean([ta tb tc td]);
old_array(:,1) = tc; % left
old_array(:,nx+2) = td; % right
old_array(1,:) = tb; % top
old_array(ny+2,:) = ta; % bottom
old_array([1 ny+2],[1 nx+2]) = 0; % zero the corners for readability
temp_array = old_array;

% This is the fractional change in temperature.
% So the while loop runs at least once, make the fractional error larger
% than any possible convergence factor.
frac = Inf;

% Keep iterating until the fractional change is less than our determined
% convergence factor.
while frac > tol
    
    % For each piece of rod, calculate the average temperature between the
    % two other pieces.
    % Since the array also includes the two non-changing values at the end,
    % start at position 2 and go until we've reached the end of the rod
    % (num+1 since it's 1-indexed, and we're starting at 2)
    for i=2:ny+1
        for j=2:nx+1
        % The old temperatures are used so that values can't blow up, in
        % certain cases.
		temp_array(i,j) = mean([mean(old_array(i,[j-1 j+1])) mean(old_array([i-1 i+1],j))]);
        end
    end
    
    % We care about the maximum fractional change for any piece. It doesn't
    % matter if the one on the end isn't changing much, if we're still
    % calculating the center pieces, keep going!
	frac = max(max(abs(temp_array - old_array) ./ old_array));
    
    % What's new is old. Take our new values and get ready to use them for
    % next time, if there is a next time.
    old_array = temp_array;
end

% Shave off the parts of the matrix that aren't the plate itself, and
% rotate it so it fits with the desired coordinate system.
temp_array = rot90(temp_array(2:ny+1,2:nx+1),3);

surf(rot90(temp_array));
colorbar;
view(2);