function [p, d, field, polar_plot,fit] = pxd(spikes, times, p_boxcar, d_boxcar, max_iter, accuracy, tol, fig)
% [p, d, field, polar_plot] = pxd(spikes, times, p_boxcar, d_boxcar, max_iter, accuracy, tol);
% plots the fields and polar_plot resulting from a maximum likelihood factorial
% model of influences of place (p) and direction (d) on cell firing, assuming Poisson noise.
% Standard (spikes/dwell_time) plots are also shown for comparison.
% See Burgess, Cacucci, Lever, O'Keefe Hippocampus (2005) for details of procedure, although
% slightly different binning and smoothing of data is used in this reference.
%
% 'spikes' and 'times' are matrices with dimensions: d_bins y_bins x_bins.
% The field and polar plot are displayed after smoothing with p_boxcar*p_boxcar and d_boxcar
% boxcar kernels respectively.
%
% The algorithm is run until the fractional change in the likelihood of the data
% is less than 'accuracy' (e.g. 0.00001, i.e. until 'convergence') or for max_iter iterations
% (e.g. 30, i.e. 'non-convergence').
%
% The value tol (e.g. 0.1) replaces expected firing rate values less than tol (usually
% due to unvisited states) to avoid log(0) divergence in loglikelhood, and to avoid division
% by zero in p_estimate and d_estimate.
% cell_name is string printed on output titles.
%
% Output: p is a matrix with dimensions: y_bins x_bins (y increases downwards),
% d is a vector with dimension d_bins.
%
% Use e.g.:
% [p, d, field, polar_plot, fit] = pxd(spikes, times, 3, 9, 30, 0.000001, 0.1, 'pc 1');

converged = 0;
if ( sum( size(spikes) == size(times) ) ~= 3 )
    fprintf(1, ' matrices spikes and times have different dimensions - spikes: %d %d %d times: %d %d %d\n',...
        size(spikes), size(times));
else
    [nd ny nx] = size(spikes);
    fit = 1;
    for iter = 1:max_iter
        if( iter==1 )
            % Initial guess
            p = ones( ny, nx);
            d = ones( nd );
        else
            p = p_estimate(d, spikes, times, tol);
        end
        d = d_estimate(p, spikes, times, tol);
        prev_fit = fit;
        % NB log likelihhod is a negative (the smaller in magnitude the better)
        fit = loglikelihood(p, d, spikes, times, tol);
        if( abs(prev_fit - fit) < -accuracy*fit )
            converged = 1;
            break;
        end
    end
    
    if( converged == 0 )
        
        fprintf(1, ' Did not converge, %d spikes, nd=%d ny=%d nx=%d.\n', sum(sum(sum(spikes))), size(spikes));
        p = NaN; d = NaN; field = NaN; polar_plot = NaN; fit = NaN;
        
    else
        % for direct comparison with regular fields, normalise area under p and d separately to
        % equal the tot no spikes / tot dwell time.
        tot_spikes = sum(sum(sum(spikes)));
        pred_spikes = sum(sum(p.*(squeeze(sum(times,1)))));
        p = p.*(tot_spikes/pred_spikes);
        pred_spikes = sum(d.*(squeeze(sum(sum(times,2),3))));
        d = d.*(tot_spikes/pred_spikes);
        
        field = find_2Dratemap(squeeze(sum(spikes, 1)), squeeze(sum(times, 1)), p_boxcar);
        
        polar_plot = find_circ_ratemap(squeeze(sum(sum(spikes, 2),3)), squeeze(sum(sum(times, 2), 3)), d_boxcar);
        
        % plot p field
        % smooth p field if requested
        % smoothing before - NB don't assume zero for unvisited positions: just ignore.
        b = ones(p_boxcar, p_boxcar);
        c = ones(size(p));
        c( squeeze(sum(times,1))==0 ) = 0;
        denom = filter2(b, c);
        denom(denom==0) = NaN;
        
        fp = filter2(b, p);
        fp = fp./denom;
        
        if fig
            subplot(2,2,3);
            imAlpha=ones(size(fp));
            imAlpha(isnan(fp))=0;
            imagesc(fp,'AlphaData',imAlpha);
            axis xy; axis off; hold on; box off; axis image;
            colormap(viridis(255));
            colorbar;
            if( nx == ny )
                axis square;
            end
            title('Corrected');
            
            % plot d polar plot. imfilter uses image processing toolbox.
            % smooth, NB d is circular, again ignore unvisited directions
            b = ones(d_boxcar);
            c = ones(size(d));
            c( squeeze(sum(sum(times,2),3))==0 ) = 0;
            denom = imfilter(c, b, 'circular');
            denom(denom==0) = NaN;
            
            fd = imfilter(d, b, 'circular');
            fd = fd./denom;
            % set unnoccupied positions to NaN to be ignored
            fd(squeeze(sum(sum(times,2),3))==0 ) = NaN;
            
            subplot(2,2,4);
            theta = ([1:nd]' - 0.5).*(2*pi/length(spikes));
            polar([theta; theta(1)], [fd; fd(1)]);
            title('Corrected');
            
            % printf regular field and polar plot for comparison
            subplot(2,2,1);
            imAlpha=ones(size(field));
            imAlpha(isnan(field))=0;
            imagesc(field,'AlphaData',imAlpha);
            axis xy; axis off; hold on; box off; axis image;
            colormap([colormap(viridis(255)); 1 1 1]);
            colorbar;
            title('Observed');
            if( nx == ny )
                axis square;
            end
            
            subplot(2,2,2);
            theta = ([1:length(polar_plot)]' - 0.5).*(2*pi/length(polar_plot));
            polar([theta; theta(1)], [polar_plot; polar_plot(1)]);
            title('Observed');
        end
    end
end

