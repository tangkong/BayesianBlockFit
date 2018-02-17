function plotx( x, st_color, y1, y2 )

if nargin < 4
   y2 = 1;
end

if nargin < 3
   y1 = 0;
end

if nargin < 2
   st_color = '-r';
end

old_next = get( gca, 'NextPlot');

hold on
v = axis;

yy_range = v(4) - v(3);

yy_start = v(3) + y1*yy_range;
yy_end   = v(3) + y2*yy_range;


for ik = 1:length(x)
   plot( x(ik)*[1 1], [ yy_start yy_end ], st_color )
end

set( gca, 'NextPlot', old_next )

return

y_scale = get( gca, 'YScale' );

if strcmp( y_scale, 'log' ) & yy_start <= 0
'here'
   yy_start = eps;
else
   'not here'
end
