% The purpose here is to establish if there is any merit in applying
% tension at something other than T=S. We let T=1..S 
shouldLoadExistingTable = 1;

filename = 'OptimalTensionOrderTable.mat';

if exist(filename,'file')
    load(filename);    
else
    slopes = [-2; -3; -4];
%     slopes = -3;
    S_range = 1:5;
    strides = [5;20;80;200];
%     strides = 20;
    
    totalSlopes = length(slopes);
    totalStrides = length(strides);
    totalEnsembles = 101;
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 5; % for whatever reason, we chose this as the primary dt
    t_damp = 30*60;
    n = 250;
    
    optimal_lambda = nan(totalSlopes,length(strides),length(S_range),length(S_range),totalEnsembles);
    optimal_rmse = nan(totalSlopes,length(strides),length(S_range),length(S_range),totalEnsembles);
    optimal_neff = nan(totalSlopes,length(strides),length(S_range),length(S_range),totalEnsembles);
        
    for iSlope = 1:length(slopes)
        slope = slopes(iSlope);
        fprintf('slope %d, ',slope);        
        
        for iStride=1:length(strides)
            stride = strides(iStride);
            dt = stride*base_dt;
            fprintf('stride %d, ',stride);
            
            for iEnsemble = 1:totalEnsembles
                fprintf('..%d',iEnsemble);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generate the signal
                sampleFactor = 10; % 2500 points takes about 1/4 second to generate maternoise
                cv=maternoise(dt/sampleFactor,sampleFactor*n,sigma_u*sqrt(2),abs(slope),1/t_damp);
                cx = cumtrapz(cv)*dt/sampleFactor;
                t_all = (dt/sampleFactor)*(0:(sampleFactor*n-1))';
                x_all = real(cx);
                data = struct('t',dt*(0:n-1)','x',x_all(1:sampleFactor:end));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generate the noise
                sigma = 10;
                noiseDistribution = NormalDistribution(sigma);
                epsilon = noiseDistribution.rand(size(data.x));
                
                x_obs = data.x + epsilon;
                t_obs = data.t;
                
                compute_rms_error = @(spline) sqrt(mean(mean(  (x_all - spline(t_all)).^2,2 ),1));
                
                % For a given ensemble, we now try fitting to various
                % spline orders with varying tension order
                for iS = 1:length(S_range)
                    S = S_range(iS);
                    for T = 1:S
                        K = S+1;
                        
                        spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'T', T, 'knot_dof', 1,'lambda',Lambda.fullTensionExpected);
                        spline.minimizeMeanSquareError(data.t,data.x);
                        
                        optimal_lambda(iSlope,iStride,iS,T,iEnsemble) = spline.lambda;
                        optimal_rmse(iSlope,iStride,iS,T,iEnsemble) = compute_rms_error(spline);
                        optimal_neff(iSlope,iStride,iS,T,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                    end                   
                end
            end
            fprintf('\n');
        end
    end
    
    save(filename, 'S_range', 'totalEnsembles', 'slopes', 'strides','optimal_lambda','optimal_rmse','optimal_neff');

end

% minimum error across S,T, for each ensemble
min_rms_error = min(min(optimal_rmse,[],3),[],4);

% now we normalize them (min error should be 1.0 for each ensemble)
rel_rms_error = optimal_rmse./min_rms_error;

pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
rel_rms_error_SvsT = zeros(length(S_range),length(S_range));
rel_rms_error_std_SvsT = zeros(length(S_range),length(S_range));
rel_rms_error_SvsT_pct_low = zeros(length(S_range),length(S_range));
rel_rms_error_SvsT_pct_high = zeros(length(S_range),length(S_range));
for iS = 1:length(S_range)
    for iT=1:length(S_range)
        values = sort(reshape(rel_rms_error(:,:,iS,iT,:),1,[]));
        rel_rms_error_SvsT(iS,iT) = mean(values);
        rel_rms_error_std_SvsT(iS,iT) = std(values);
        rel_rms_error_SvsT_pct_low(iS,iT) = values(ceil( ((1-pct_range)/2)*length(values)));
        rel_rms_error_SvsT_pct_high(iS,iT) = values(floor( ((1+pct_range)/2)*length(values)));
    end
end

fprintf('\n\n');
fprintf('\\begin{tabular}{l *{%d}{l}}\n',length(S_range));
fprintf('\\toprule & \\multicolumn{%d}{c}{S} \\\\ \n',length(S_range));
fprintf('\\cmidrule(lr){2-%d} \n',length(S_range)+1);
fprintf('T ');
for iS = 1:length(S_range)
    fprintf('& %d ',iS);
end
fprintf('\\\\ \\midrule \n');
for iT = 1:length(S_range)
    fprintf('%d ',iT);
    for iS = 1:length(S_range)
        if isnan(rel_rms_error_SvsT_pct_low(iS,iT))
            fprintf('& ');
        else
            fprintf('& %.1f-%.1f\\%% ', 100*(rel_rms_error_SvsT_pct_low(iS,iT)-1),100*(rel_rms_error_SvsT_pct_high(iS,iT)-1) );
        end
    end
    fprintf('\\\\ \n');
end
fprintf(' \\bottomrule \n\\end{tabular} \n');

return

% now compute the ensemble average increase in error
% rel_rms_error = mean(rel_rms_error,5);

% now flatten across stride and slope
% rel_rms_error_SvsT =squeeze(mean(mean(rel_rms_error)));

% fprintf('\n\n');
% fprintf('\\begin{tabular}{r | llll} stride & full dof & reduced dof & blind initial & blind optimal \\\\ \\hline \\hline \n');
% for iSlope = 1:length(slopes)
%     
%     fprintf('%d slope &&&&  \\\\ \\hline \n',slopes(iSlope));
%     for iStride=1:length(result_stride)
% %         fprintf('%d & %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_blind_initial_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), dmse_reduced_dof_blind_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
%         fprintf('%d & %#.3g m$^2$ (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_true_optimal_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_blind_initial_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), 100*dmse_reduced_dof_blind_optimal_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
%     end
%     
% end
% fprintf('\\end{tabular} \n');

rms_error_SequalT = zeros(totalSlopes,length(strides),length(S_range),totalEnsembles);
for iS = 1:length(S_range)
    rms_error_SequalT(:,:,iS,:) = optimal_rmse(:,:,iS,iS,:);
end

% minimum error across S,T, for each ensemble
min_rms_error = min(rms_error_SequalT,[],3);

% now we normalize them (min error should be 1.0 for each ensemble)
rel_rms_error = rms_error_SequalT./min_rms_error;

% now compute the ensemble average increase in error
rel_rms_error = mean(rel_rms_error,4);

rel_rms_error_SlopeVsS = 100*(squeeze(mean(rel_rms_error,2)) - 1);

fprintf('\n\n');
fprintf('\\begin{tabular}{r | llll}  & S=1 & S=2 & S=3 & S=4 \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    
    fprintf('$\\omega^{%d}$ & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%%  \\\\ \n',slopes(iSlope),rel_rms_error_SlopeVsS(iSlope,1),rel_rms_error_SlopeVsS(iSlope,2),rel_rms_error_SlopeVsS(iSlope,3),rel_rms_error_SlopeVsS(iSlope,4));    
end
fprintf('\\end{tabular} \n');
