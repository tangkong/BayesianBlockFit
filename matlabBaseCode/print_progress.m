function did_it = print_progress( id, num, skip, cpu_time )
% function print_progress( id, num, st_add, skip )
% id = runing index
% num = max value
% skip
% cputime

if nargin < 4
    st_left = '';
else
	cpu_now = cputime;
	time_so_far = ( cpu_now - cpu_time );
	rate = time_so_far / id;
	time_left = rate * ( num - id );
	st_left = sprintf(' Time left: %4.2f min  = %4.2f hours.', ...
        time_left / 60, time_left / 3600 );
end

did_it = 0;

if nargin > 2
    if rem( id, skip ) ~= 0
        return
    end
end

fprintf(1, ['%6.0f out of %6.0f ' st_left '\n'], id, num );
did_it = 1;
return