function  data_out = find_bumps( data_in )

num_pars_1 = 1 + 3;
num_pars_2 = 1 + 8*2*3;% Fudge factor for BIC

xx = data_in.xx;
xx = xx(:);
data_out.xx_raw = xx;% raw data truncated

fit_order = data_in.fit_order;
ncp_prior = data_in.ncp_prior;

id_length = length(xx);

% process the raw data and the polynomical fit 
id_start =  12; % shuck off a few points at both ends - was 12
id_end = id_length - 50;
%id_end   = 962;
% id_end   = 950;


xx = xx(id_start:id_end);

data_out.xx_trunc = xx;% raw data truncated
data_out.id_start = id_start;% raw data truncated
data_out.id_end   = id_end;% raw data truncated

num_use = length( xx );
tt = (1:num_use)' / num_use;% evenly spaced times
pp_fit = polyfit( tt, xx, fit_order );
data_out.pp_fit = pp_fit;

small_fac = 1.e-6;
xx_poly_fit = polyval( pp_fit, tt );
data_out.xx_poly_fit = xx_poly_fit;

xx = xx - xx_poly_fit;% remove polynomial trend
data_out.xx_use = xx;% data used in the BB algorithm

% guess at noise level
id_low = find( xx <= median( xx ) );% select low values
%sigma_guess = std( xx( id_low ) );  % take noise as the STD of these
sigma_guess = 1*std( xx( id_low ) );  % take noise as the STD of these

% data input to BB
cell_data = xx;
cell_data(:,2)    = sigma_guess * ones( size( xx ) );
data_in.cell_data = cell_data;
data_in.ncp_prior = ncp_prior;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_out_blocks = find_blocks( data_in );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate_vec      = data_out_blocks.rate_vec;
change_points = data_out_blocks.change_points;

num_blocks = length( rate_vec );
cp_use = [ 1 change_points length( xx ) ];
data_out.change_points = cp_use;

num_cp = length( cp_use ) - 1;
id_left_vec = zeros( num_cp, 1 );
ii_rite_vec = zeros( num_cp, 1 );

for ii_cp = 1: num_cp

    ii_1 = cp_use( ii_cp );
    ii_2 = cp_use( ii_cp + 1 );
    id_left_vec( ii_cp ) = ii_1;
    id_rite_vec( ii_cp ) = ii_2;

end
  
data_out.id_left_vec = id_left_vec;
data_out.id_rite_vec = id_rite_vec;

% find the maxima defining the watersheds ("bumps")
id_max = zeros( 1, num_blocks );

% Find the indices of the highest neighbors of each block
for ii = 1: num_blocks

    ii_left = ii - 1; % prevent out of bounds
    if ii_left <1;ii_left = 1;end
    
    ii_rite = ii + 1;
    if ii_rite > num_blocks;ii_rite = num_blocks;end

    rate_left  = rate_vec( ii_left );
    rate_cent  = rate_vec( ii );
    rate_right = rate_vec( ii_rite );

    [ ~, ii_this_max ] = max( [ rate_left rate_cent rate_right ] );
    id_max( ii ) = ii + ii_this_max - 2;

end

% adjust any outliers near the edges
id_max( find( id_max > num_blocks ) ) = num_blocks;
id_max( find( id_max < 1 ) ) = 1;

% Implement hill climbing (the HOP algorithm)
index_hop = 1: num_blocks;% Initial state: each block points to itself

while 1
       index_new = id_max( index_hop );% Each block hops to highest neighbor
       
       id_change = find( index_new ~= index_hop );% Locate index changes
       if isempty( id_change );break;end      % If no changes escape loop
       index_hop = index_new;                     % Implement the changes
end

id_max = unique( index_hop );
num_max = length( id_max );% Number of maxima in block representation
    
% collect data for each bump
clear data_struc_bumps
ii_count_max = 0;

for ii_max = 1: num_max
 
    id_this_max = id_max( ii_max );                % index of max block
    
    % find datum index of maximum
    id_datum_max_bump = fix( ( id_left_vec( id_this_max ) + ...
                               id_rite_vec( id_this_max ) )/2 );

    id_this_vec = find( index_hop == id_this_max );% indices of block(s)
    num_blocks_here = length( id_this_vec );

    % find dataum indices for the start and end of the bump
    id_first_block = id_this_vec(1);
    id_datum_left = id_left_vec( id_first_block );
    
    id_last_block = id_this_vec( num_blocks_here );
    id_datum_rite = id_rite_vec( id_last_block );

    % Now fit Gaussians 
  
    xx_all = xx( id_datum_left: id_datum_rite );
    %%
    num_all = length( xx_all );
    ii_fit = (1:num_all)';
    %small = 0.01 * max( xx_all );
    xx_offset = min( xx_all );% - small;
    xx_use = xx_all - xx_offset; 
    
    % These fit routines are part of the MatLab curve fitting toolbox
    options_1 = fitoptions('gauss1', 'Lower', [ 0 0 0 ]);
    options_2 = fitoptions('gauss2', 'Lower', [ 0 0 0 0 0 0 ]);
    [ fit_1, good_1 ] = fit( ii_fit, xx_use, 'gauss1', options_1 );
    
   % if length( xx_use ) >=  6 % if too few points, fit routine crashes 
     if length( xx_use ) >=  6 % if too few points, fit routine crashes
    
        ii_count_max = ii_count_max + 1;
        
        % These outputs are indices in the block array
        data_struc_bumps( ii_count_max ).id_block_max = id_this_max;% Block at the maximum
        data_struc_bumps( ii_count_max ).id_block_vec = id_this_vec;% Blocka in this bump

        % These outputs are indices in the adjusted data array
        % Datum at the maximum restoring the ofset due to id_start
        data_struc_bumps( ii_count_max ).id_datum_max  = id_datum_max_bump + id_start - 1;
        data_struc_bumps( ii_count_max ).id_datum_left = id_datum_left + id_start - 1;
        data_struc_bumps( ii_count_max ).id_datum_rite = id_datum_rite + id_start - 1;
    
        
        [ fit_2, good_2 ] = fit( ii_fit, xx_use, 'gauss2', options_2 );
        %[ fit_3, good_3 ] = fit( ii_fit, xx_use, 'gauss3');
        
        sigma_1 = good_1.rmse;
        sigma_2 = good_2.rmse;
        %sigma_3 = good_3.rmse;

        bic_vec(1) = num_all * log( sigma_1 .^ 2 ) + num_pars_1*log( num_all );
        bic_vec(2) = num_all * log( sigma_2 .^ 2 ) + num_pars_2*log( num_all );
        %bic_vec(3) = num_all * log( sigma_3 .^ 2 ) + num_pars_3*log( num_all );

        [ bic_min, num_gaussians ] = min( bic_vec );

        data_struc_bumps( ii_count_max ).num_gaussians = num_gaussians;

        if num_gaussians == 1
            
            fit_best = fit_1;

            a1 = fit_best.a1;
            b1 = fit_best.b1;
            c1 = fit_best.c1;

            xx_fit_1 = a1*exp(-((ii_fit-b1)/c1).^2);
            data_struc_bumps( ii_count_max ).xx_fit_1 = xx_fit_1;

        elseif num_gaussians == 2

            fit_best = fit_2;

            a1 = fit_best.a1;
            b1 = fit_best.b1;
            c1 = fit_best.c1;

            a2 = fit_best.a2;
            b2 = fit_best.b2;
            c2 = fit_best.c2;

            xx_fit_1 = a1*exp(-((ii_fit-b1)/c1).^2);
            xx_fit_2 = a2*exp(-((ii_fit-b2)/c2).^2);
            data_struc_bumps( ii_count_max ).xx_fit_1 = xx_fit_1;
            data_struc_bumps( ii_count_max ).xx_fit_2 = xx_fit_2;

        %elseif num_gaussians == 3
        %    fit_best = fit_2;
        end

        xx_best = feval( fit_best, ii_fit );
        data_struc_bumps( ii_count_max ).xx_best = xx_best + xx_offset;
        data_struc_bumps( ii_count_max ).fit_best = fit_best;
        
    else
        disp('this one is too short ')
    end

end
   
data_out.data_struc_bumps = data_struc_bumps;
