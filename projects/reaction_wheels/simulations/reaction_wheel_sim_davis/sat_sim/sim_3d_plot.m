function sim_3d_plot( op, R, V, P )
%SIM_3D_PLOT  Display body and a vector in 3-D.
%   V is body vertices. Solid line is from origin to P in inertial space; 
%   dashed line is from origin to -P in inertial space.
%   Body face colors RGB correspond to +x, +y, +z directions, 
%   respectively, in body coordinates.

persistent hx hy hz hk hp hn
switch op
	case 'init'
		V = R*V;
		P = R*P;
		F = [5 6 7 8; 2 3 7 6; 3 4 8 7];
		figure('position', [600 357 575 725] )
			axis([-1 1 -1 1 -1 1]*20); hold on
			hx = patch( V(1,[1 2 6 5]), V(2,[1 2 6 5]), V(3,[1 2 6 5]), 'r' );
			hy = patch( V(1,[1 4 8 5]), V(2,[1 4 8 5]), V(3,[1 4 8 5]), 'g' );
			hz = patch( V(1,[1 2 3 4]), V(2,[1 2 3 4]), V(3,[1 2 3 4]), 'b' );
			hk = patch( 'Faces', F, 'Vertices', V', ...
				'Facecolor', 'white', 'FaceAlpha', .75 );
			hp = plot3( [P(1) 0], [P(2) 0], [P(3) 0], 'k', 'LineWidth', 2 );
			hn = plot3( [-P(1) 0], [-P(2) 0], [-P(3) 0], '--k', 'LineWidth', 2 );
		set(gca, 'ZDir', 'reverse', 'YDir', 'reverse' )
		xlabel('X')
		ylabel('Y')
		zlabel('Z')
		axis vis3d
		grid
		rotate3d on
%		view(-48,15)

	case 'update'
		V = R*V;
		P = R*P;
		set( hx, 'Vertices', V(:,[1 2 6 5])' )
		set( hy, 'Vertices', V(:,[1 4 8 5])' )
		set( hz, 'Vertices', V(:,[1 2 3 4])' )
		set( hk, 'Vertices', V' )
		set( hp, 'XData', [P(1) 0], 'YData', [P(2) 0], 'ZData', [P(3) 0] );
		set( hn, 'XData', [-P(1) 0], 'YData', [-P(2) 0], 'ZData', [-P(3) 0] );
		drawnow
	end
end
