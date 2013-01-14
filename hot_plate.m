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
%    2013-01-14

%{
Additional Documentation:

MATH

    This program attempts to solve the partial differential equation
    d^2 T / dx^2 = 0 iteratively. It does this by breaking the length of X
    into a given amount of small pieces, and averaging the value of each
    piece between the pieces that surround it. The reasoning behind the
    average was determined by expressing the "second difference" of each of
    the temperature values in an equation, and rearranging to determine
    what the value of the current piece is when the second difference is 0.

    In other words, for a piece i compared to the piece on the left (L) and
    a piece on the right (R), the second difference between the values
    should be zero:

        (TR - Ti) - (Ti - TL) = 0

    or when rearranged:

        Ti = (TR + TL) / 2

    hence the average.

    New values are computed until the maximum fractional change for any
    piece is less than the determined convergence factor, which is 0.01
    divided by the number of pieces the rod is broken into.


VARIABLE NAMING
    
    See the comments near where each variable is declared / invoked the
    first time.

%}




% Check starting values :
if min([ta tb tc td]) <= 0
    error('Temperatures are given in Kelvin, and as such must be above absolute zero.')
elseif min([nx ny]) < 1 || rem([nx ny],1) ~= 0
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

fprintf('Starting Temps: %u', old_array);

% This is the fractional change in temperature.
% So the while loop runs at least once, make the fractional error larger
% than any possible convergence factor.
frac = 1;

% This is the convergence factor.
% The data needs to be more accurate as the rod is broken into more pieces.
convergence = tol / (nx * ny);

% Keep iterating until the fractional change is less than our determined
% convergence factor.
while frac > convergence
    
    % For each piece of rod, calculate the average temperature between the
    % two other pieces.
    % Since the array also includes the two non-changing values at the end,
    % start at position 2 and go until we've reached the end of the rod
    % (num+1 since it's 1-indexed, and we're starting at 2)
    for i=2:num+1
        
        % The old temperatures are used so that values can't blow up, in
        % certain cases.
		temp_array(i) = (old_array(i-1) + old_array(i+1)) / 2;
    end
    
    % We care about the maximum fractional change for any piece. It doesn't
    % matter if the one on the end isn't changing much, if we're still
    % calculating the center pieces, keep going!
	frac = max(abs(temp_array - old_array) ./ old_array);
    
    % What's new is old. Take our new values and get ready to use them for
    % next time, if there is a next time.
    old_array = temp_array;
end

% Shave off the temperatures of the left and right ends, giving back only
% the temperatures of the pieces of the rod.
temp_array = temp_array(2:num+1);

% DON'T FORGET TO TRANSLATE COORDINATE SYSTEMS