% main_bumps.m

% for id_case_run = 1:2
%     
%     % get list of the data files
%     if id_case_run == 1
%          data_files = dir('./mehta/Sample*B2_6*csv');
%     else
%          data_files = dir('./mehta/Sample*B2_14*csv');
%     end



%dir_root = ('C:\matlab_scripts\Bayesianbumps\mehta\');
dir_root =('C:\work\HiTp\July2017\Takeuchi\1D\');
Qmin = 1.5;
Qmax = 5.5;
minPkwidth = 3; %min pixel width for a peak

dir_name = [dir_root,  '*.csv'];
%mkdir(dir_root, 'output\');
mkdir(dir_root, 'PkFit\images\');
data_files = dir(dir_name);

    [ num_files, ~ ] = size( data_files );
    fprintf(1,'Found%5.0f data files.\n', num_files )
    
    fit_order = 2;% fit quadratic to full data interval
    ncp_prior = .5;%
    
    num_bumps_vec = zeros( num_files, 1 );
    cpu_0 = cputime;
    num_bumps_max = 100;% arbitrary maximum number of peaks
    
    id_block_max_mat = zeros( num_files, num_bumps_max );
    id_datum_max_mat = zeros( num_files, num_bumps_max );
    id_datum_left_mat = zeros( num_files, num_bumps_max );
    id_datum_rite_mat = zeros( num_files, num_bumps_max );
    
    
maxPk_mat(num_files, 6) = zeros;

    for ii_file = 1: num_files

         count_peaks = 0;
         
        coeff_mat = zeros( num_bumps_max, 3 );
         back_mat = zeros( num_bumps_max, 1 );

        file_name_this = data_files( ii_file ).name;
        id_cut = findstr('.csv', file_name_this );
        file_name_use = file_name_this(1:id_cut(1)-1 );
        id_bad = findstr('_', file_name_use );
        file_name_use( id_bad ) = '-';

        id_cut = findstr( '.csv', file_name_this );
        var_name = file_name_this(1: id_cut(1) - 1 );

%         Load the data
%         st_load = ['load ./mehta/' file_name_this ];
        filen_load = [dir_root file_name_this];
        load(filen_load);

        % Put data into array xx
        st_move = ['xx = ' var_name ';' ];
        eval( st_move );
                
        scalef = (xx(size(xx, 1), 1) - xx(1, 1))/size(xx, 1);
            startQ = xx(1, 1) - scalef;

        % call the bump finder
        clear data_bump_in
        data_bump_in.xx = xx(:,2);
        data_bump_in.fit_order = fit_order;
        data_bump_in.ncp_prior = ncp_prior;

        data_bump_out = find_bumps( data_bump_in );

             xx_poly_fit = data_bump_out.xx_poly_fit;
                id_start = data_bump_out.id_start;
        data_struc_bumps = data_bump_out.data_struc_bumps;
        
        num_bumps = length( data_struc_bumps );
        num_bumps_vec( ii_file ) = num_bumps;

        for ii_bump = 1: num_bumps

            id_block_max = data_struc_bumps( ii_bump ).id_block_max;

            id_datum_max  = data_struc_bumps( ii_bump ).id_datum_max;
            id_datum_left = data_struc_bumps( ii_bump ).id_datum_left;
            id_datum_rite = data_struc_bumps( ii_bump ).id_datum_rite;

            id_block_max_mat( ii_file, ii_bump ) = id_block_max;
            id_datum_max_mat( ii_file, ii_bump ) = id_datum_max;
            id_datum_left_mat( ii_file, ii_bump ) = id_datum_left;
            id_datum_rite_mat( ii_file, ii_bump ) = id_datum_rite;

        end
        
        figure(1);clf
        subplot(2,1,1)

        col_vec = 'bgrmc';
        col_count = 0;

        id_shift = data_bump_out.id_start - 1;
        xx_raw_plot = data_bump_out.xx_raw;
        vert_shift = 0.1*( max( xx_raw_plot ) - min( xx_raw_plot) );
        
        % Plot the fit
        set(plot(xx(:,1),  xx_raw_plot + vert_shift, '-k' ),'LineWidth',2)
           xlabel('Q (1/A)')
       %set(plot(xx_raw_plot + vert_shift, '-k' ),'LineWidth',2)
        hold on
    %plotx( id_datum_max_mat( ii_file, 1: num_bumps ), ':r')
     plotx((scalef*(id_datum_max_mat( ii_file, 1: num_bumps ))+startQ),  ':r')
   
        num_bumps = length( data_struc_bumps );
        for ii_bump = 1: num_bumps -1

            id_datum_left = data_struc_bumps( ii_bump ).id_datum_left;
            id_datum_rite = data_struc_bumps( ii_bump ).id_datum_rite;
            id_datum_max  = data_struc_bumps( ii_bump ).id_datum_max;
                 fit_best = data_struc_bumps( ii_bump ).fit_best;

            num_gaussians = data_struc_bumps( ii_bump ).num_gaussians;
            xx_best = data_struc_bumps( ii_bump ).xx_best;
            if id_datum_rite > length( xx_poly_fit )
                id_datum_rite = length( xx_poly_fit);% total bandaid!!
            end
            id_best = (id_datum_left: id_datum_rite);
            id_bump_offset = id_datum_left - 1;

            xx_back = xx_poly_fit( id_best );% background from poly fit

            col_count = col_count + 1;
            col_use = col_vec( col_count );
            if col_count >= length( col_vec);col_count = 0;end
            if length( xx_best ) > length( xx_back )
                xx_best = xx_best(1: length( xx_back ) );
            end

           set(plot((scalef* id_best) + startQ,  xx_best + xx_back, ['-' col_use ] ),'LineWidth',1)
           
       % set(plot( id_best,  xx_best + xx_back, ['-' col_use ] ),'LineWidth',1)
           
            set(gca,'TickDir','out')

            subplot(2,1,2)

            % Compute Gaussian fit (for either one or two components)
            
            if num_gaussians == 1

                count_peaks = count_peaks + 1;
                a1 = fit_best.a1;
                b1 = fit_best.b1;
                c1 = fit_best.c1;

                coeff_mat( count_peaks, : ) = [ a1, b1 + id_bump_offset, c1 ];
                back_mat( count_peaks ) = mean( xx_back );

                xx_fit_1 = data_struc_bumps( ii_bump ).xx_fit_1;
                %xx_fit_1 = xx_fit_1 / sum( abs( xx_fit_1 ) );
               % %normalization
                if length( xx_fit_1 ) > length( id_best )
                    xx_fit_1 = xx_fit_1( 1: length( id_best ) );
                end
                
               
                %set(plot( id_best,  xx_fit_1, ['-' col_use ] ),'LineWidth',1)
                set(plot( (scalef* id_best) + startQ,  xx_fit_1, ['-' col_use ] ),'LineWidth',1)
                hold on

            elseif num_gaussians == 2

                a1 = fit_best.a1;
                b1 = fit_best.b1;
                c1 = fit_best.c1;
                count_peaks = count_peaks + 1;
                coeff_mat( count_peaks, : ) = [ a1, b1 + id_bump_offset, c1 ];
                back_mat( count_peaks ) = mean( xx_back );

                a2 = fit_best.a2;
                b2 = fit_best.b2;
                c2 = fit_best.c2;
                count_peaks = count_peaks + 1;
                coeff_mat( count_peaks, : ) = [ a2, b2 + id_bump_offset ,c2 ];
                back_mat( count_peaks ) = mean( xx_back );

                xx_fit_1 = data_struc_bumps( ii_bump ).xx_fit_1;
                xx_fit_2 = data_struc_bumps( ii_bump ).xx_fit_2;
%                 xx_fit_1 = xx_fit_1 / sum( abs( xx_fit_1 ) );
%                 xx_fit_2 = xx_fit_2 / sum( abs( xx_fit_2 ) );

%                 set(plot( id_best,  xx_fit_1, ['-' col_use ] ),'LineWidth',1)
%                 set(plot( id_best,  xx_fit_2, ['-' col_use ] ),'LineWidth',1)
                 set(plot( (scalef* id_best) + startQ,  xx_fit_1, ['-' col_use ] ),'LineWidth',1)
                set(plot(  (scalef* id_best) + startQ,  xx_fit_2, ['-' col_use ] ),'LineWidth',1)
                
                hold on
                
            end
            set(gca,'TickDir','out')
            subplot(2,1,1)


        end

        v = axis;
        % Write number of components at the top of the plot
        for ii_bump = 1: num_bumps -1

            id_datum_left = data_struc_bumps( ii_bump ).id_datum_left;
            id_datum_rite = data_struc_bumps( ii_bump ).id_datum_rite;
            id_datum_max  = data_struc_bumps( ii_bump ).id_datum_max;
            num_gaussians = data_struc_bumps( ii_bump ).num_gaussians;
            hh=text( id_datum_max, v(4) - 0.05*(v(4)-v(3)), int2str( num_gaussians ) );
            hh=text((scalef*id_datum_max) +startQ, v(4) - 0.05*(v(4)-v(3)), int2str( num_gaussians ) );
           % (scalef* id_best) + startQ,
            set(hh,'FontWeight','bold','FontSize',12)

        end

        subplot(2,1,2)
        v=axis;
        v(3) = -.01;
        axis(v)

        % Shorten the arrays to the actual number of components
        [ num_act, ~ ] = size( coeff_mat );

        if count_peaks < num_act
            coeff_mat = coeff_mat(1: count_peaks, : );
             back_mat =  back_mat(1: count_peaks, : );
        end
        
        % Save Figure
        title(file_name_use);
                    
         filen_fig = [dir_root 'PkFit\images\' file_name_use '_fit.png' ];
         saveas(gcf, filen_fig); 


        % Save ascii data files 
        
%         if id_case_run == 1

%             st_save_coeff = ['save -ascii ./mehta_out/' file_name_use '_coeff.txt coeff_mat  ' ]
%             st_save_back  = ['save -ascii ./mehta_out/' file_name_use '_back.txt back_mat ' ];

% %         else
% % 
% %             st_save_coeff = ['save -ascii ./mehta_out_2/' file_name_use '_coeff.txt coeff_mat  ' ]
% %             st_save_back  = ['save -ascii ./mehta_out_2/' file_name_use '_back.txt back_mat ' ]

       % end
            filen_save_coeff = [dir_root 'PkFit\' file_name_use '_coeff.csv' ];
            %filen_save_back  = [dir_root 'output\' file_name_use '_back.txt'];
            
            
            coeff_mat(:, 3) = 2.366* (scalef*coeff_mat(:, 3)); %FWHM
            coeff_mat(:, 2) = (scalef*coeff_mat(:, 2)) + startQ;
            coeff_mat = [coeff_mat , back_mat];
           
              %  save(filen_save_coeff, 'coeff_mat', '-ascii');
              %save(filen_save_back, 'back_mat', '-ascii'); 
              
            fid = fopen(filen_save_coeff, 'wt');
            fprintf(fid, '%s,%s,%s,%s\n','intensity', 'position(Q)', 'FWHM(Q)', 'background');
            fprintf(fid, '%f,%f,%f,%f\n', coeff_mat');   %transpose is important!
            fclose(fid);
                
%         eval( st_save_coeff )
%         eval( st_save_back )

%Stop and Elaspsed time%

            %find the strongest peak
            range = find(coeff_mat(:, 2) > Qmin & coeff_mat(:, 2) < Qmax & coeff_mat(:, 3) > (minPkwidth *2.355*scalef) );
            maxind = find(coeff_mat(range,1) == max(coeff_mat(range, 1)));
       
                numpk = length(range);
               
%                 if numpk == 0
%                     maxPK_mat(ii_file, 3:6) = 0;
%                 else
              
                    maxPk_mat(ii_file, 1) = ii_file;
                    maxPk_mat(ii_file, 2) = numpk;
                    maxPk_mat(ii_file, 3) = coeff_mat(range(maxind), 1);
                    maxPk_mat(ii_file, 4) = coeff_mat(range(maxind), 2);
                    maxPk_mat(ii_file, 5) = coeff_mat(range(maxind), 3);
                    maxPk_mat(ii_file, 6) = coeff_mat(range(maxind), 4);
                   % maxPK_mat(ii_file, 3) = coeff_mat(range(maxind), 1)
%                end
                    
                 
               
    end
 filen_save_maxPk  = [dir_root 'output\' file_name_use '_maxPk.txt'];
  save(filen_save_maxPk, 'maxPk_mat', '-ascii');
 
%end % id_case_run
