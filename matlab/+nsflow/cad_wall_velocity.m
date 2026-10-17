function velocity = cad_wall_velocity(coordinates,time,component,settings)
%CAD_WALL_VELOCITY Prescribed heave, consistent with the moving immersed mask.
velocity=zeros(size(coordinates{1}));if component==2,velocity(:)=settings.amplitude*settings.omega*cos(settings.omega*time);end
end
