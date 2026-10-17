function inputFcn = pressure_history(time,pressure)
%PRESSURE_HISTORY Interpolated measured/CFD pressure difference, Pa.
% pressure: Nchannel x Ntime; positive lower-minus-upper, acting +w on plate.
% Time extrapolation is rejected. No flow/structure feedback is introduced.
time=time(:)';
assert(numel(time)>1 && all(isfinite(time)) && all(diff(time)>0) && ...
    size(pressure,2)==numel(time) && all(isfinite(pressure(:))), ...
    'ael:Pressure','Require increasing time and finite channel histories.');
inputFcn=@sample;
    function p=sample(t)
        assert(isscalar(t) && isfinite(t) && t>=time(1) && t<=time(end), ...
            'ael:PressureTime','Pressure history does not cover requested time.');
        p=interp1(time,pressure',t,'linear')';
    end
end
